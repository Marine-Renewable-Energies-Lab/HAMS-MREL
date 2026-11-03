# Mean Wave Drift (Phases 0 to 2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add mean wave drift forces and moments to HAMS-MREL: far-field for the system, near-field per body, for single and multi-body runs and for unidirectional and bidirectional waves.

**Architecture:** New modules work on one symmetry-free, body-agnostic view of the mesh and the first-order solution, built by small adapters from the existing single-body and multi-body arrays. A driver is called inside the existing heading loops, after the motion solve. Everything is off unless `SecondOrder.in` exists.

**Tech Stack:** Fortran 90+ (ifort/ifx, MKL for LAPACK), OpenMP, fpm, test-drive. Python only for the mesh generator.

**Spec:** `docs/superpowers/specs/2026-11-04-second-order-waves-design.md`. Sections 5.2, 5.3 and 5.4 hold the conventions and formulas. Read them first.

## Global Constraints

- Without `SecondOrder.in`, all outputs must stay unchanged. `config1` and `config1_gmres` integration tests must keep passing after every task.
- Implement from the equations here and in the spec. Do not open, copy or translate Nemoh source (GPL-3.0-or-later; this project is Apache-2.0).
- Internal time convention is e^{−iωt}. Conjugate exactly once, in `SecondOrderOutput`.
- Use `Const_mod` values of `G` and `RHO` through arguments. Use `NEWUNIT` for new files and `outputdir` for paths, never a hard-coded `Output/`.
- New files get the Apache header used in `src/PotentWavForceMulti.f90`. Use free-form Fortran with `implicit none` and `real(8)`/`complex(8)`.
- `src/PatclVelct.f90` is not valid UTF-8, so preserve its encoding if it must ever be edited. No task needs to edit it.
- Do not change first-order behaviour. Existing issues are guarded against, not fixed.
- Acceptance thresholds are proposals from the spec. If a threshold fails, investigate and report. Do not loosen it silently.
- Work on branch `second-order-drift`, created from `origin/main` at `c40ad4d` (v2.3.2). Commit there. Do not push or open a PR unless asked. End commit messages with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- test-drive runs tests in parallel under OpenMP. Each test needs its own files and directory. A test that sets module globals, such as `WaveDyn_mod`, goes in its own suite run serially, by extending the `parallel=` condition in `test/unit_tester.f90`.
- Integration test commands must go through `native_path` from `testing_utilities`, because Windows CI runs them through `cmd.exe`. Follow `test/test_integration.f90`: a setup suite runs HAMS once, and a separate suite checks results with relative tolerances. Use `test_oamass_files_are_close`, since the absolute comparison no longer exists.
- CI now runs the unit and integration tests on Linux and Windows, so every test must pass on both, and within runner time limits.

**Toolchain.** Tests need Intel Fortran, MKL and fpm, per `docs/developer_documentation.md`:

```bash
export FPM_FC=ifx FPM_FFLAGS="-qmkl -qopenmp -heap-arrays 0" FPM_LDFLAGS="-qmkl -qopenmp"
fpm build && fpm install --prefix install
fpm test --target unit-tester
fpm test --target integration-tester -- install/bin/hamsmrel
```

On Windows, CI uses `set FPM_FFLAGS=/Qmkl /heap-arrays0 /fpp` and `install/bin/hamsmrel.exe`. See `.github/workflows/build_and_test_windows.yml`.

On this machine, oneAPI 2022.1 is at its default path, and `fpm` 0.13.0 is at `..\tools\fpm\fpm.exe`, outside the repository and not on PATH. In a `cmd` script started from Git Bash, clear `NoDefaultCurrentDirectoryInExePath` before `setvars.bat`. Task 1 starts by confirming a working baseline. If unavailable, stop and tell the user rather than skipping tests.

## Review Focus

1. Hull file contains panels lying in z = 0 (a lid or deck), or the waterline sits slightly off z = 0: expect a clear stop or correct tolerance, not silently missing waterline. Pinned in Task 3.
2. A gap in the waterline from a mesh defect: expect a warning naming the body and the gap, and the run continues. Pinned in Task 3.
3. `SYBO=1` adds zero and infinite frequencies: expect them skipped with no rows or NaNs in drift files. Pinned in Task 5.
4. A fully submerged body in an array with surface-piercing bodies: expect valid near-field results with term I zero for that body. Pinned in Task 8.
5. A multi-body run with a rotated `LCS` or `REFL≠1`: expect drift in global axes and scaling by `REFL`. Pinned in Tasks 5 and 8.

---

## File Structure

| File | Responsibility |
|---|---|
| `src/SecondOrderTypes.f90` | Options, `panel_mesh`, `first_order_state`, `surface_geom`, result types |
| `src/SecondOrderInput.f90` | Read and validate `SecondOrder.in` |
| `src/FirstOrderState.f90` | Incident wave helpers, image expansion, HAMS-to-physical conversion |
| `src/SurfaceGeometry.f90` | Neighbours, waterline, per-body waterplane and volume properties |
| `src/SurfaceVelocity.f90` | ∇φ by differencing, and by source solve |
| `src/FarFieldDrift.f90` | Kochin function and far-field F_x, F_y, M_z |
| `src/NearFieldDrift.f90` | Terms I to V |
| `src/WaterlineField.f90` | Field-point influence row, BIE waterline potential |
| `src/SecondOrderOutput.f90` | WAMIT-format writers |
| `src/SecondOrderDriver.f90` | Single and multi-body adapters, orchestration, calls from `HAMS_Prog.f90` |
| `test/test_second_order_*.f90` | Unit tests. `test/second_order_fixtures.f90` holds meshes and analytic fields. |
| `test/test_integration_drift.f90`, `test/data/drift_*` | Integration cases |
| `pre_postprocessing_python/SecondOrder/make_mesh.py` | Cylinder and box `.pnl` generator |

Each new `src/` file must also be added to `CompileHAMS.bat` line 19 (pure modules after `WavDynSubs.f90`, the driver immediately before `HAMS_Prog.f90`) and to `HAMS_MREL.vfproj`.

---

### Task 1: Input file and validation

**Files:**
- Create: `src/SecondOrderTypes.f90`, `src/SecondOrderInput.f90`, `test/test_second_order_input.f90`
- Modify: `src/HAMS_Prog.f90:137-146`, `CompileHAMS.bat:19`, `HAMS_MREL.vfproj`, `test/unit_tester.f90`, `docs/input_files.md`

**Interfaces:**
- Produces, in `SecondOrderTypes`:
  ```fortran
  type :: so_options
    logical :: active = .false.
    integer :: mean_drift = 0          ! 0 off, 1 unidirectional, 2 all heading pairs
    logical :: do_far = .false., do_near = .false.
    integer :: velocity_method = 1     ! 1 differencing, 2 source
    integer :: waterline_method = 1    ! 1 extrapolate, 2 field
    real(8) :: waterline_offset = 0.01d0
    logical :: motion_fixed = .false.
    integer :: n_azimuth = 360
    real(8) :: far_tolerance = 1d-4
    integer :: write_terms = 1         ! 0 none, 1 term file, 2 also panel diagnostics
    real(8) :: wl_z_tolerance = 1d-6   ! relative to body size
  end type
  ```
- Produces, in `SecondOrderInput`:
  - `subroutine read_second_order_options(inputdir, opts, ok, msg)`: absent file gives `opts%active=.false.`, `ok=.true.`.
  - `subroutine validate_second_order(opts, isol, h, isolv, irsp, has_gen_modes, ok, msg)`.

- [ ] **Step 1: Confirm the baseline.** Run the toolchain commands above before changing anything. Expected: all existing tests pass. Record any failure before continuing.

- [ ] **Step 2: Write failing tests** in `test/test_second_order_input.f90`, using a temp directory under the executable's path:
  - `absent_file_inactive`: `active` is false, `ok` is true.
  - `defaults_applied`: a file with only `mean_drift = 1` gives `n_azimuth == 360`, `velocity_method == 1`.
  - `all_keys_parsed`: file using the keys below gives the given values, with `!` comments and any key order.
  - `unknown_key_rejected`: `ok` is false and `msg` contains the key.
  - `validate_rejects`: each is `ok=.false.` with a specific message: `isol=2`; `h=600d0`; near with `has_gen_modes`; `velocity_method=2` with `isolv=2`; `velocity_method=2` with `irsp=1`.
  - `validate_accepts`: `h=-1d0`, `h=70d0`, far-only with `has_gen_modes`.

  Keys: `mean_drift`, `drift_methods` (words `far`, `near`), `velocity_method` (`differencing`, `source`), `waterline_method` (`extrapolate`, `field`), `waterline_offset`, `motion_state` (`free`, `fixed`), `n_azimuth`, `far_tolerance`, `write_terms`, `wl_z_tolerance`.

- [ ] **Step 3: Register and run.** Add `collect_second_order_input` to `test/unit_tester.f90`. Run `fpm test --target unit-tester`. Expected: fail to compile, because the module is missing.

- [ ] **Step 4: Implement** the types and both subroutines. `h > 500` is rejected with a message telling the user to set `H<0`, because `VINP` and `DINP` differ there (`src/PatclVelct.f90:64,106`) and first-order behaviour stays unchanged.

- [ ] **Step 5: Wire in.** In `HAMS_Prog.f90`, after `ReadControlFile` succeeds, read and validate, and stop with `msg` on failure. `has_gen_modes` is an `INQUIRE` over `gen_mod_{b}_{d}.txt`. Update the build lists. Document the file in `docs/input_files.md`.

- [ ] **Step 6: Run all tests.** Expected: all pass, and existing integration results are unchanged.

- [ ] **Step 7: Commit.**
```bash
git add src/SecondOrderTypes.f90 src/SecondOrderInput.f90 src/HAMS_Prog.f90 test CompileHAMS.bat HAMS_MREL.vfproj docs/input_files.md
git commit -m "feat: add optional SecondOrder.in input and validation"
```

### Task 2: First-order state

**Files:**
- Create: `src/FirstOrderState.f90`, `test/test_first_order_state.f90`, `test/second_order_fixtures.f90`
- Modify: `src/SecondOrderTypes.f90`, `test/unit_tester.f90`, build lists

**Interfaces:**
- Produces, in `SecondOrderTypes`:
  ```fortran
  type :: panel_mesh                ! always full: symmetry images expanded
    integer :: n = 0, nnode = 0, nbody = 1, nbase = 0, nsys = 1
    real(8), allocatable :: xc(:,:), nrm(:,:), area(:), psize(:)  ! (3,n),(3,n),(n),(n)
    real(8), allocatable :: nj(:,:)     ! (6,n) as used in the radiation BC, may be generalized
    real(8), allocatable :: xnode(:,:)  ! (3,nnode)
    integer, allocatable :: ncn(:), ncon(:,:), body(:), base(:), image(:)  ! (n),(4,n),(n),(n),(n)
    real(8), allocatable :: xr(:,:)     ! (3,nbody)
  end type
  type :: first_order_state
    real(8) :: omega, k, nu, h, amp, beta, g, rho   ! h<0: deep
    complex(8), allocatable :: phi_i(:), phi_p(:), dphi_p_dn(:), dphi_dn(:) ! (n)
    complex(8), allocatable :: grad_i(:,:), grad_p(:,:)  ! (3,n), grad_p filled by Task 6
    complex(8), allocatable :: xi(:,:)        ! (6,nbody), e^{-iwt}, translations then rotations
  end type
  ```
- Produces, in `FirstOrderState`:
  - `pure subroutine depth_factor(k, z, h, zf, dzf)`: Z and dZ/dz, as exponentials in deep water and overflow-safe for large kh. Finite depth: Z = cosh k(z+h)/cosh kh.
  - `pure complex(8) function incident_phi(x, omega, k, h, amp, beta, g)`: −i(g·amp/ω) Z e^{ik(x cosβ + y sinβ)}. No `XW` argument, since `ReadHydroStatic*` always sets `XW=0`.
  - `pure subroutine incident_grad(x, omega, k, h, amp, beta, g, grad)`.
  - `subroutine expand_images(mesh_base_arrays..., mesh)`: builds `panel_mesh`. Image rule from `EFORCE`: with `ISX=1` image 2 is x→−x, otherwise y→−y. Flip `xc`, `nrm`. `nj` flip factors for image 2: n₁,n₂,n₃ follow the coordinate; n₄ gets the y factor, n₅ the x factor, n₆ the product (as in `EFORCE_MULTI:75-80`).
  - `subroutine build_state(mesh, omega, k, nu, h, amp, beta, g, rho, mxpot, dspl, motion_fixed, state)`: `mxpot(nbase, 6*nbody+1, nsys)` and `dspl(6*nbody)`.

  Conversion, from spec section 5.2:
  ```
  xi(j)      = conj(dspl(j))                       (zero if motion_fixed)
  phi_p(i)   = i*mxpot(b,md,im) + sum_j (-i*omega) xi(j) mxpot(b,j,im),  md = 6*nbody+1
  dphi_dn(i)   = sum over the panel's body j=1..6 of (-i*omega) xi(j) nj(j,i)
  dphi_p_dn(i) = dphi_dn(i) - grad_i(:,i).nrm(:,i)
  ```
  Radiation sums run over all bodies, since `mxpot` for body b's modes is non-zero on other bodies' panels.

- [ ] **Step 1: Write failing tests:**
  - `incident_elevation`: (iω/g)·`incident_phi` at z = 0, x = 0 equals `amp`, finite and deep.
  - `incident_grad_fd`: `incident_grad` matches central differences of `incident_phi` to 1e-6 relative.
  - `incident_laplace_and_fs`: FD Laplacian near 0, and ∂φ/∂z = ν φ at z = 0 to 1e-6.
  - `incident_large_kh`: `k=5, h=400` gives finite values equal to deep-water values to 1e-10.
  - `matches_hams_vinp_dinp`, in a serial suite because it sets `WaveDyn_mod` globals: `i*VINP` and `i*DINP` match for `h=70` and for `h=-1`.
  - `build_state_conjugates`: 2-panel synthetic, `dspl=(1,2)` etc. gives `xi` as the conjugate, and `phi_p`, `dphi_dn` match hand values.
  - `build_state_fixed`: `motion_fixed` gives `xi=0`, `dphi_dn=0`.
  - `expand_images_y`: cube half mesh with `ISY=1`: `n=2*nbase`, image positions and normals mirrored, and `sum(area*nrm)` over a closed cube is zero.

- [ ] **Step 2: Run.** Expected: compile failure.

- [ ] **Step 3: Implement.** Add to fixtures `make_cylinder(radius, draft, ntheta, nz, nr, mesh, with_bottom)`, normals into the body.

- [ ] **Step 4: Run all tests.** Expected: pass.

- [ ] **Step 5: Commit** `feat: add first-order state conversion for second-order loads`.

### Task 3: Surface geometry

**Files:**
- Create: `src/SurfaceGeometry.f90`, `test/test_surface_geometry.f90`
- Modify: `src/SecondOrderTypes.f90`, `test/unit_tester.f90`, build lists

**Interfaces:**
- Produces, in `SecondOrderTypes`:
  ```fortran
  type :: surface_geom
    integer, allocatable :: nbr_start(:), nbr(:)       ! CSR, size n+1 and nnz
    integer :: nwl = 0
    real(8), allocatable :: wl_x1(:,:), wl_x2(:,:), wl_xm(:,:), wl_len(:) ! (3,nwl)...
    real(8), allocatable :: wl_nh(:,:), wl_nz(:)       ! unit horizontal normal into body; panel n_z
    integer, allocatable :: wl_panel(:), wl_body(:)
    real(8), allocatable :: awp(:), xf(:), yf(:), l11(:), l12(:), l22(:), vol(:), xb(:), yb(:), zb(:) ! per body, about xr
    integer, allocatable :: n_wl_gaps(:)
  end type
  ```
- Produces: `subroutine build_surface_geom(mesh, wl_z_tol, sharp_deg, geom, ok, msg)`, with `sharp_deg` default 45.

Rules:
- Merge nodes by coordinate with a tolerance of 1e-6 × body size (sort by x and sweep). Merged nodes join symmetry images.
- An edge is on the waterline if both nodes have |z| ≤ tol and it belongs to exactly one panel.
- Neighbours are panels sharing a merged node, plus panels of the same body with centroid within 1.5·`psize`, for non-conforming meshes. Exclude any with normal angle above `sharp_deg`.
- Waterplane properties come from Green's theorem on waterline segments, and volume from the divergence theorem on panels. Use signs so that `awp` and `vol` are positive given normals into the body.

- [ ] **Step 1: Write failing tests:**
  - `cylinder_waterline_length`: Σ`wl_len` equals the polygon perimeter to 1e-12 relative, and each node appears in exactly two segments.
  - `waterline_normal_into_body`: `wl_nh` is unit, horizontal and dot (axis − `wl_xm`) > 0.
  - `flared_normal_is_unit`: a cone-like hull with n_z ≠ 0 keeps |`wl_nh`| = 1 and `wl_nz` matches.
  - `multi_body_loops`: two cylinders give the right `wl_body` and length per body.
  - `submerged_body_no_waterline`: cylinder at z < 0 gives `nwl = 0` for it and `ok`.
  - `symmetry_images_join`: half cylinder after `expand_images` matches the full-mesh perimeter, with no gaps.
  - `nonconforming_neighbours`: a T-junction mesh gives neighbours across the junction.
  - `sharp_edge_excluded`: on a box, a side panel's neighbours exclude the bottom panels.
  - `wp_properties_box`: box a×b×d offset from `xr` gives exact `awp`, `xf`, `yf`, `l11`, `l12`, `l22`, `vol`, `xb`, `yb`, `zb`.
  - `lid_panels_rejected` (Review Focus 1): a panel in z = 0 gives `ok=.false.` and `msg` naming the panel.
  - `waterline_off_zero` (Review Focus 1): a hull with top nodes at z = −1e-8·size works, while top nodes at z = −0.01·size give zero segments for a surface-piercing-looking open hull, and `msg` warns.
  - `waterline_gap_reported` (Review Focus 2): one removed segment gives `n_wl_gaps(1) > 0` and `ok`.
  - `triangles_supported`: mixed triangles and quads pass the length test.

- [ ] **Step 2: Run.** Expected: compile failure.
- [ ] **Step 3: Implement** `build_surface_geom`.
- [ ] **Step 4: Run all tests.** Expected: pass. Also check the mesh in `examples/example5_single_body/Input/HullMesh.pnl` through a throwaway read: 30 waterline segments and no gaps.
- [ ] **Step 5: Commit** `feat: add waterline extraction and panel neighbours`.

### Task 4: Far-field drift

**Files:**
- Create: `src/FarFieldDrift.f90`, `test/test_far_field.f90`
- Modify: `src/SecondOrderTypes.f90`, `test/second_order_fixtures.f90`, `test/unit_tester.f90`, build lists

**Interfaces:**
- Consumes: `panel_mesh`, `first_order_state`, `depth_factor`.
- Produces:
  - `subroutine far_amplitude(mesh, s, theta, a, da)`: a(θ) and ∂a/∂θ, per unit `amp` when `s%amp=1`.
  - `pure real(8) function group_ratio(k, h)`: C_g/C_p.
  - `subroutine far_field_pair(mesh, s1, s2, opts, f)`: `f(3)` complex for (F_x, F_y, M_z), per unit amplitude pair, with M_z about the global origin. `s1` and `s2` share ω.
- Fixture: `maccamy_fuchs_state(mesh, a, h, k, omega, amp, beta, g, rho, s)`, setting `phi_i`, `phi_p`, `dphi_p_dn`, `grad_i`, and `grad_p` analytically on r = a, and `maccamy_fuchs_a(theta, ...)`.

Formulas (spec 5.3), with `C = k` for deep water, otherwise `C = k²/((k²−ν²)h + ν)`, and E = Z(ζ)e^{−ik(x cosθ + y sinθ)}:

```
H(θ) = Σ_i [dphi_p_dn(i) E − phi_p(i) (n·∇E)] area(i),   a(θ) = −(ωC/2g) √(2/(πk)) e^{−iπ/4} H(θ)
κ = (ρg/2) group_ratio
F_12 = κ[ −√(2π/k)( e_β1 conj(a2(β1)e^{iπ/4}) + e_β2 a1(β2)e^{iπ/4} ) − ∫ a1 a2* e_θ dθ ]
M_12 = −(κ/k)[ ∫(1/2i)(a2* ∂a1 − a1 ∂a2*)dθ + √(2π/k)( i conj(∂a2(β1)e^{iπ/4}) − i ∂a1(β2)e^{iπ/4} ) ]
```

Integrate with the trapezoid rule on [0, 2π), which converges fast for periodic integrands. Start at `n_azimuth` and double until the change is below `far_tolerance` (absolute after scaling by ρg, or relative). Stop doubling at 2¹⁴ points and warn. Cache a(θ) per state.

MacCamy-Fuchs on r = a, with β = 0 and ε₀ = 1, εₙ = 2:

```
φ_total = −i(gA/ω) Z(z) Σ_n ε_n i^n [2i/(π k a H'_n(ka))] cos nθ
a_series(θ) = −A √(2/(πk)) e^{−iπ/4} Σ_n ε_n [J'_n(ka)/H'_n(ka)] cos nθ
```

Use `bessel_jn`, `bessel_yn` with H'_n(x) = H_{n−1}(x) − (n/x)H_n(x), H₋₁ = −H₁, truncated at n ≈ ka + 20.

- [ ] **Step 1: Write failing tests,** using a bottom-mounted cylinder mesh (a = 5, h = 10, no bottom panels, ω = 1), no BEM:
  - `kochin_matches_series`: max|a − a_series|/max|a_series| < 1e-2.
  - `drift_matches_analytic`: F_x vs κ∫|a_series|²(1−cosθ)dθ within 2%.
  - `energy_identity`: X = ∫|a|²dθ within 1%, with X = −2√(2π/k) Re[a(β)e^{iπ/4}].
  - `lateral_zero`: |F_y|, |M_z/a| below 1e-6 of F_x.
  - `deep_limit`: for one prescribed panel field on a truncated cylinder, a(θ) with `h=-1` and with kh = 30 differs by less than 1e-8 relative.
  - `derivative_fd`: `da` matches FD of `a`.
  - `conjugate_symmetry`: F(β₁,β₂) = conj F(β₂,β₁) to 1e-10.
  - `diagonal_consistent`: `far_field_pair(s,s)` is real to 1e-12 and equals the single-heading form.
  - `reference_shift`: translate mesh by d with phases adjusted (state for shifted body, incident phase e^{ik d·e_β}): F unchanged, M_z changes by −(d×F)_z, to 1e-6 relative to F·|d|.
  - `group_ratio_limits`: ½ deep, 1 as kh → 0.

- [ ] **Step 2: Run.** Expected: compile failure.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run all tests.** Expected: pass. A failure of `drift_matches_analytic` or `energy_identity` means a normalization error. Fix the cause, not the threshold.
- [ ] **Step 5: Commit** `feat: add far-field mean drift`.

### Task 5: Driver, output and first end-to-end run

**Files:**
- Create: `src/SecondOrderDriver.f90`, `src/SecondOrderOutput.f90`, `test/test_integration_drift.f90`, `test/data/drift_cyl/`, `pre_postprocessing_python/SecondOrder/make_mesh.py`
- Modify: `src/HAMS_Prog.f90` (lines 237-253 and 300 single-body; 570-580 and 724 multi-body; final writes), `test/integration_tester.f90`, build lists, `docs/output_files.md`, `docs/python_scripts.md`

**Interfaces:**
- Consumes: Tasks 1 to 4.
- Produces, in `SecondOrderDriver`:
  - `subroutine so_setup_single(opts, outputdir)` and `so_setup_multi(opts, outputdir)`: copy module arrays into `panel_mesh` (after `CalNormals` and after `BuildCombinedArrays`), build geometry, open files.
  - `subroutine so_after_motion_single(opts, kk, ib)` and `so_after_motion_multi(opts, kk, ib)`: build the state, compute results. Skip when `ABS(TP)<1e-6` or `ABS(TP+1)<1e-6`.
  - `subroutine so_end_frequency(opts, kk)`: evaluate pairs when `mean_drift=2`, write rows.
- Produces, in `SecondOrderOutput`: `subroutine write_drift_row(unit, oufr, beta1_deg, beta2_deg, index, value, refl, moment)`: writes `PER BETA1 BETA2 I Mod Pha Re Im` format `2ES14.6` plus `ES14.6,I6,4ES14.6`. `value` is per unit amplitude pair. The routine conjugates it and divides by ρg·`refl` for forces, or ρg·`refl²` for moments.
- File: `Wamit_format/MeanDriftMomentum.8`, I = 1, 2, 6, system total.
- With `write_terms=2`: `SecondOrder/PanelDiagnostics.txt` with ω, β, panel, φ_D = φ_I + i·`mxpot(b,md,im)`, φ_total, plus first-order equation residual per frequency.

`make_mesh.py cylinder --radius --draft --ntheta --nz --nr --out`, and `box`, writes `HullMesh.pnl` in the format in `docs/input_files.md`, with a normal check on generated panels. Inputs for `test/data/drift_cyl`: cylinder a = 10, draft = 10, depth 100, 300 to 500 panels, 3 frequencies, headings 0° and 30°, no external loads, `motion_state=free`. Add `Hydrostatic.in` with consistent mass and restoring. Compute hydrostatics in the script, and document the values.

- [ ] **Step 1: Write failing integration tests** in `test/test_integration_drift.f90`, following `test/test_integration.f90`, parsing `.8`, `.1`, `.3`, `.6p`:
  - `run_drift_cyl`: exit code 0 and output files exist.
  - `convention_residual`: from diagnostics, |[−ω²(M+A) − iωB + C]ξ − F|/|F| < 1e-8, with F recovered as `i*EXFC`.
  - `diffraction_matches_6p`: with field points placed at panel centroids, conj(iωρφ_D/(ρgA)) equals the `PressureElevation.6p` diffraction value to 1e-4.
  - `radiated_energy`: ½ω²B_jj vs ½ρgC_g∫|a_j|²dθ within 1%, for j = 1, 3, 5 (add diagnostics for a_j).
  - `haskind`: excitation from the Kochin function at β+π vs `ExcForce.3`, within 1%. Derive the prefactor from X_j = iωρ∬(φ_I n_j − ψ_j ∂φ_I/∂n)dS with Green's second identity, and take magnitude first if the phase is in doubt.
  - `energy_balance`: free body, X vs ∫|a|²dθ within 1%.
  - `symmetric_lateral_zero`: at heading 0, F_y and M_z are below 1e-3 F_x in magnitude.
  - `sybo_skipped` (Review Focus 3): case with `SYBO=1` has no rows for limit frequencies and no `NaN` anywhere in `.8`.
  - `refl_scaling` (Review Focus 5): same case with `REFL=2`: forces halve and M_z quarters relative to `REFL=1`.
  - `no_input_no_output`: `config1` run creates no `MeanDrift*` file.
  - `gmres_agrees`: `ISOLV=2` vs 1, F_x within 1e-2 relative.

- [ ] **Step 2: Run** `fpm test --target integration-tester -- install/bin/hamsmrel`. Expected: fails because outputs are missing.

- [ ] **Step 3: Implement** the driver, output and calls in `HAMS_Prog.f90`. For `mean_drift=2`, keep a `first_order_state` per heading within the frequency. Pass `CRS`/`MATX`/`AMAS`/`BDMP` only for diagnostics.

- [ ] **Step 4: Run all tests.** Expected: pass, including unchanged `config1`.

- [ ] **Step 5: Document** `.8`, including the `OUFR` first column, in `docs/output_files.md`, and the script in `docs/python_scripts.md`.

- [ ] **Step 6: Commit** `feat: output far-field mean drift from single and multi-body runs`.

**Phase 1 is complete here.** Report results to the user before continuing.

### Task 6: Surface velocity by differencing

**Files:**
- Create: `src/SurfaceVelocity.f90`, `test/test_surface_velocity.f90`
- Modify: `test/unit_tester.f90`, build lists

**Interfaces:**
- Consumes: `surface_geom`, `first_order_state`.
- Produces: `subroutine velocity_differencing(mesh, geom, s, n_fallback)`: fills `s%grad_p`. Total velocity is `grad_i + grad_p`.

Algorithm: for each panel, build orthonormal tangents t₁, t₂. For neighbours j, Δ = (xc_j − xc_i) projected to the tangents. Solve the 2×2 weighted least-squares fit (weight 1/|Δ|²) for ∂φ_P/∂t₁ and ∂φ_P/∂t₂ from φ_j − φ_i. Then `grad_p = n·dphi_p_dn + t1·g1 + t2·g2`. With fewer than 2 non-collinear neighbours (test the 2×2 determinant), keep only the normal part and count `n_fallback`.

- [ ] **Step 1: Write failing tests:**
  - `linear_field_exact`: φ = c·x on a flat plate gives tangential gradient exact to 1e-10.
  - `known_normal_preserved`: grad_p·n equals `dphi_p_dn` to round-off.
  - `incident_field_converges`: treat φ_I on a cylinder as φ_P with exact ∂φ/∂n: error vs `incident_grad` falls by at least 1.5 per doubling of resolution.
  - `waterline_row_one_sided`: top-row panels have error below 3× the interior mean.
  - `fallback_counted`: an isolated panel gives `n_fallback=1` and no NaN.

- [ ] **Step 2: Run.** Expected: compile failure.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run all tests.** Expected: pass.
- [ ] **Step 5: Commit** `feat: add body velocity by tangential differencing`.

### Task 7: Near-field terms

**Files:**
- Create: `src/NearFieldDrift.f90`, `test/test_near_field.f90`
- Modify: `src/SecondOrderTypes.f90`, `test/unit_tester.f90`, build lists

**Interfaces:**
- Consumes: `first_order_state` with `grad_p`, `surface_geom`.
- Produces:
  - `pure complex(8) function pair_minus(f1, g1, f2, g2)`: `0.25*(f1*conjg(g2) + conjg(f2)*g1)`. This is the only product used, so later QTFs can reuse it (spec 5.4).
  - `subroutine waterline_extrapolate(mesh, geom, s, phi_wl)`: `phi_wl(nwl)` total at z = 0, φ_I exact at z = 0 plus φ_P from the adjacent panel with grad_p.
  - `subroutine near_field_pair(mesh, geom, s1, s2, phi_wl1, phi_wl2, terms)`: `terms(6, 6, nbody)` complex, index (component, term I, Ia, II, III, IV, V), about `xr` of each body. The total is the sum over terms.

Per panel of body b, with x relative to `xr(:,b)`, r = (n) or (x×n) for force or moment, φ_t amplitude −iωφ, ∇φ = grad_i + grad_p, using `pair_minus` (P) with components i and j from state 1 and state 2:

```
I:   ½ρg ∫ r' [P(ζ,ζ) − 2 P(ζ,X3)] dl          r' = n' or x_wl×n'
Ia:  ½ρg ∫ r'' [n_z/√(1−n_z²)] [P(ζ−X3, ζ−X3)] dl   r'' = k or x_wl×k
II:  −ρ ∬ r ½ Σ_i P(∂iφ, ∂iφ) dS
III: −ρ ∬ r Σ_i P(X_i, ∂iφ_t) dS
IV force:  Σ ε_ijk P(α_j, G_k),  G = −ρ ∬ n φ_t dS
IV moment: Σ ε_ijk P(α_j, H_k) + Σ ε_ijk P(ξ_j, G_k),  H = −ρ ∬ (x×n) φ_t dS
```

P(·,·) denotes the pair function over both components of each factor. For a symmetric bracket such as P(ζ,X₃), use f = ζ, g = X₃. V uses the same P on the products in eq. 15.56 and 15.57, with Z_o = `xr(3,b)`, and x_f, y_f, L_ij, ∀, x_b, y_b relative to `xr`:

```
F_S = −ρg A { α1α3 xf + α2α3 yf + ½(α1²+α2²) Zo } k
M_S,x = ρg { (ξ2α2 − ξ3α3) A xf − ξ2α1 A yf + α1α2 A xf Zo − (3/2 α1² + ½ α2²) A yf Zo − ξ2ξ3 A − ξ3α1 A Zo
           + α2α3 (L11−L22) − 2α1α3 L12 + α1α2 V xb − ½(α1²+α3²) V yb }
M_S,y = ρg { −ξ1α2 A xf + (ξ1α1 − ξ3α3) A yf + (½α1² + 3/2 α2²) A xf Zo − α1α2 A yf Zo + ξ1ξ3 A − ξ3α2 A Zo
           + α1α3 (L11−L22) + 2α2α3 L12 + ½(α2²+α3²) V xb }
M_S,z = 0
```

ξ and α in these formulas come from `xi(1:3,b)`, `xi(4:6,b)`. Rotations in HAMS are about `XR`, which matches O here.

- [ ] **Step 1: Write failing tests:**
  - `pair_mean_identity`: `pair_minus(f,g,f,g)` = ½Re(f g*), real to round-off.
  - `pair_time_average`: for random complex f₁, f₂, g₁, g₂ and A₁, A₂, a numerically time-averaged product of f(t) and g(t) over a period equals |A₁|²P₁₁ + |A₂|²P₂₂ + 2Re(A₁A₂*P₁₂), to 1e-10.
  - `fixed_cylinder_analytic`: `maccamy_fuchs_state`, exact `grad_p`, exact waterline φ: Σ terms F_x vs the Task 4 analytic drift, within 2%; III, IV, V are zero to round-off.
  - `term3_uniform`: ∇φ uniform vector V, X uniform: III_z = −ρ·½Re(X*·(−iωV))·A_wp on a closed-below cylinder, to 1e-3.
  - `term4_cross`: synthetic: G and α known give ½Re(conj(α)×G) to round-off.
  - `hydrostatic_exact_box`: wall-sided box offset from origin, still water (ζ = 0, φ = 0). Compare V + Ia with a direct integration of ρg × submerged volume under a finite rotation (roll-pitch-yaw order, eq. 15.49) and translation, using a 400×400 column grid. Second-order mean for amplitude ε equals ½[½(F(+ε) + F(−ε)) − F(0)], evaluated with ε = 1e-3 for each single mode and for pairs. Agreement within 1e-3 relative for F_z, M_x, M_y. A failure here may be a transcription error in eq. 15.56 and 15.57, so re-read the manual first.
  - `flare_only_vertical`: Ia has zero horizontal force. Ia is zero on a wall-sided hull.
  - `form_equivalence`: on a free cylinder, I + IV equals the alternative form I' + IV', with (ζ−X₃)² and rotation including −Cξ-type force from `A_wp`, within 1e-3 of F_x magnitude.
  - `pair_symmetry`: terms(s1,s2) = conj(terms(s2,s1)) to round-off.
  - `reference_point`: M_z about a shifted `xr` with rigid motion recomputed: M_z changes by −(d×F)_z to 1e-6 relative.

- [ ] **Step 2: Run.** Expected: compile failure.
- [ ] **Step 3: Implement.** Keep each term in its own helper.
- [ ] **Step 4: Run all tests.** Expected: pass.
- [ ] **Step 5: Commit** `feat: add near-field mean drift terms`.

### Task 8: Near-field end to end and cross-checks

**Files:**
- Modify: `src/SecondOrderDriver.f90`, `src/SecondOrderOutput.f90`, `test/test_integration_drift.f90`, `docs/output_files.md`
- Create: `test/data/drift_cyl_fixed/`, `test/data/drift_two_cyl/`, `test/data/drift_half/`, `test/data/drift_mixed/`, `test/data/drift_rotated/`

**Interfaces:**
- Consumes: `velocity_differencing`, `waterline_extrapolate`, `near_field_pair`.
- Produces:
  - `Wamit_format/MeanDriftPressure.9`: index 6(b−1)+I, about each body's `XR`, in global axes.
  - `Wamit_format/MeanDriftTerms.txt`: `PER BETA1 BETA2 body I term Re Im`.
  - When both methods run, the log prints the far-near relative difference per frequency and heading for F_x, F_y, M_z. M_z about global origin requires adding `xr × F` per body.

- [ ] **Step 1: Write failing integration tests:**
  - `fixed_near_far`: `motion_state=fixed` cylinder: F_x within 5% (see the mesh note below).
  - `free_near_far`: free cylinder, F_x within 5% away from resonance. Use only the frequencies in the fixture.
  - `near_mz_far_mz`: an asymmetric heading on a box or offset cylinder: near M_z about origin vs far, within 10%. This also checks the far-field M_z form.
  - `external_stiffness`: add surge and heave `KSTF`: near and far still agree within 5%.
  - `two_body_total`: two cylinders: Σ near F_x, F_y vs far total, within 5%.
  - `far_apart_isolated`: separation of 50 wavelengths: each body's F_x within 5% of the isolated cylinder's.
  - `half_vs_full`: `ISY=1` half mesh vs full mesh: F_x within 1e-3 relative.
  - `submerged_mix` (Review Focus 4): a surface-piercing and a submerged body: term I and Ia are zero for the submerged body, and its total force is finite.
  - `rotated_lcs` (Review Focus 5): a fixed box rotated by 30° via `LCS` with heading +30° vs unrotated at heading 0: forces are the rotation of each other, within 1e-3.
  - `bidirectional`: `mean_drift=2`: F(β₁,β₂) = conj F(β₂,β₁) in file, and diagonal rows equal the `mean_drift=1` run.
  - `irsp_stability`: `IRSP=1` vs 0 at a non-irregular frequency: F_x within 1%.

  Mesh note: use the 300 to 500 panel fixture only if run time is acceptable and thresholds hold. Otherwise make the CI case coarse with looser thresholds, and put the 5% checks in Task 11's study. Record which was chosen.

- [ ] **Step 2: Run.** Expected: fail on missing outputs.
- [ ] **Step 3: Implement** the driver branch and writers. Compute `grad_p` only when near-field is on.
- [ ] **Step 4: Run all tests.** Expected: pass. Investigate any failure with `MeanDriftTerms.txt`, starting from the fixed case.
- [ ] **Step 5: Document** `.9` and the terms file.
- [ ] **Step 6: Commit** `feat: output per-body near-field mean drift`.

### Task 9: Waterline potential by field evaluation

**Files:**
- Create: `src/WaterlineField.f90`, `test/test_waterline_field.f90`
- Modify: `src/SecondOrderDriver.f90`, `test/test_integration_drift.f90`, build lists

**Interfaces:**
- Consumes: `SGLINTBD_TRI2`, `SGLINTBD_QUAD2`, their `_MULTI_` versions, `INFGREEN3D`, `FINGREEN3D`.
- Produces:
  - `subroutine field_influence_row(xp, base_mesh_selector, row_g, row_dg)`: `row_g(nbase,nsys)`, `row_dg(nbase,nsys)`, the integrated G and ∂G/∂n_q, built like `DUM(1:2)` in `CalPotentialMulti`, with the field point mirrored per image. The driver wraps the single and multi Rankine routines behind this.
  - `subroutine waterline_field(mesh, geom, s, offset, rows_g, rows_dg, phi_wl)`: φ_P = (1/4π)Σ_images Σ_j [G ∂φ_P/∂n − (∂G/∂n) φ_P]. Rows are computed once per frequency, then applied to every state.

Evaluation point: `wl_xm − offset·wl_len·wl_nh`, with z = −offset·`wl_len`, which lies in the fluid. φ_I is added exactly at z = 0.

- [ ] **Step 1: Write failing tests** (run inside integration, since they need Green functions):
  - `field_matches_existing`: at an off-body point, `waterline_field` gives the same φ_P as `CalPotentialMulti` (scattered: subtract F0) to 1e-8.
  - `offset_sensitivity_recorded`: offsets 0.005, 0.01, 0.05, 0.1 on the fixed cylinder give term I within 5% of each other and of the MacCamy-Fuchs value. Print the table.

- [ ] **Step 2: Run.** Expected: fail, since the key is unsupported.
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run all tests.** Expected: pass, or a documented failure of the sensitivity bound for Task 11 to resolve.
- [ ] **Step 5: Commit** `feat: add field-evaluated waterline potential`.

### Task 10: Surface velocity by source solve

**Files:**
- Modify: `src/SurfaceVelocity.f90`, `src/SecondOrderDriver.f90`, `test/test_integration_drift.f90`

**Interfaces:**
- Produces: `subroutine velocity_source(mesh, s_list, cgrn, rkbn, ...)`: fills `grad_p` for all states of a frequency, using one factorization of A_s per symmetry block.

Formulation (spec D2), with n out of the fluid, φ_P = (1/4π)∬σG, G's derivatives as stored being with respect to the source q:

```
(2π I + D) σ = 4π ∂φ_P/∂n,   D_ij = n_i · ∫_j ∇_p G dS
∇_p φ_P(i) = (1/4π)[ 2π σ_i n_i + Σ_j σ_j ∫_j ∇_p G dS ]
```

∇_pG from stored arrays, for field image IS, in that image's coordinates:
- FLAG = 1, panel-integrated Rankine `DIJ` covers 1/R: ∇_p = −`DIJ` for all three components.
- Wave and image part `GRN·DS`: ∇_p,h = −GRN(2:3), ∂z_p = +GRN(4).
- FLAG = 0, where `GRN` also includes 1/R: remove the 1/R z-part by adding 2(ZF−ZP)/R³ to GRN(4) before use.

Reuse `CGRN*`/`RKBN*` for the current frequency. Assemble as a dense matrix for each symmetry block using the `RXY` pattern of `ASSB_LEFT`, and call `ZGETRF` and `ZGETRS`. `GRN_KIND=4` precision applies. `ISOLV=2` and `IRSP=1` stay rejected by Task 1's validation, unless the Task 11 study requires them.

- [ ] **Step 1: Write failing integration tests:**
  - `source_normal_consistent`: on each panel, grad_p·n recovers ∂φ_P/∂n within 1%.
  - `source_potential_consistent`: (1/4π)ΣσG vs direct-BIE φ_P, max relative difference below 5% on the fixture.
  - `source_near_far`: fixed and free cylinders, F_x within 5% of far-field.
  - `source_symmetry`: half vs full mesh within 1e-3.

- [ ] **Step 2: Run.** Expected: fail, since the method is unsupported.
- [ ] **Step 3: Implement.** Add a timing line to the log.
- [ ] **Step 4: Run all tests.** Expected: pass.
- [ ] **Step 5: Commit** `feat: add body velocity by auxiliary source solve`.

### Task 11: Convergence study, defaults and documentation

**Files:**
- Create: `pre_postprocessing_python/SecondOrder/drift_convergence.py`, `docs/mean_drift.md`
- Modify: `src/SecondOrderInput.f90` (defaults, only if the study supports a change), `docs/superpowers/specs/2026-11-04-second-order-waves-design.md`, `readme.md`, `docs/developer_documentation.md`, `test/test_integration_drift.f90`

- [ ] **Step 1: Run the study.** The script generates meshes and runs HAMS-MREL for these combinations, and writes a table:
  - Geometry: bottom-mounted cylinder (analytic reference), truncated cylinder, box.
  - Mesh: three resolutions, panels per side doubling.
  - Velocity: differencing, source. Waterline: extrapolate, field.
  - Settings: `GRN_KIND` 4 and 8 (rebuild required), `ISOLV` 1 and 2 (differencing only), `IRSP` 0 and 1, including a frequency near an irregular frequency.
  - Record F_x error against analytic or far-field, and run time.

- [ ] **Step 2: Apply the gate.** Near-field and far-field must agree and converge on at least one free body. Report this to the user with the table. Set default methods only from this evidence. If differencing fails the gate, make `source` the default, then implement `IRSP` support for it as a follow-up.

- [ ] **Step 3: External reference: Nemoh.** Run Nemoh v3 as a separate program, built outside this repository, and compare only through its documented input and output files. Do not add Nemoh source or binaries to the repository, and do not run it in CI. Compare stored numbers instead, after confirming with the user that committing them is acceptable.
  1. Setup already done on 2026-11-04. Nemoh v3.0.3 is built in `..\Nemoh-3.0.3\build` with ifort and static MKL, and `8b_QTF_Cylinder` has been run in `..\Nemoh-runs\8b_mesh634`. Rebuild recipe: from a `cmd` shell, clear `NoDefaultCurrentDirectoryInExePath` (Git Bash sets it, which breaks `setvars.bat`), run `setvars.bat intel64 vs2019`, then `cmake -S. -Bbuild -G "NMake Makefiles" -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_C_COMPILER=cl -DALLOW_MKL=ON -DBLA_VENDOR=Intel10_64lp_seq -DBLA_STATIC=ON` and `cmake --build build`. Run stages from a copy of the case, after creating `Mechanics`, `Motion`, `results\QTF`, `results\sources`, `QTFPreprocOut` and `mesh`, because Nemoh's own folder creation fails on Windows for at least `Mechanics`. Run `hydrosCal` only for cases that use it.
     First result: on the `8b` diagonal for ω ≥ 0.9 rad/s, Nemoh and HydroSTAR's near-field F_x differ by a median of 9.7% and M_y by 16.9%, on this 634-panel mesh and without Nemoh's hydrostatic term. The bundled `DriftFxFF` and `DriftFxMF` files contain only zeros, so they give no far-field reference. The bundled HydroSTAR near-field drift can also be compared with HAMS on the same geometry.
  2. Add `--format nemoh` to `make_mesh.py`, so both solvers use the same geometry and panel layout. Use a wall-sided waterline, because Nemoh does not renormalize its waterline normal. Use no lid and `IRSP=0`, no external loads, and frequencies away from irregular frequencies.
  3. Run `DUOK` only (contributing terms = 1) on a uniform grid, and read the diagonal ω₁ = ω₂, β₁ = β₂ of `results/QTF/OUT_QTFM_N.dat`, with values multiplied by ρg. Include the hydrostatic flag when comparing F_z, M_x and M_y.
  4. First compare first-order A, B, X and RAOs, because a first-order difference invalidates everything after it.
  5. Compare the fixed cylinder, then the free cylinder, then a bidirectional pair. Fixed needs `motion_state=fixed`, or for Nemoh a very large external stiffness. For each mesh, compare F_x, F_y and M_z, and F_z, M_x and M_y as a second tier. Refine the mesh until each solver has converged before applying a threshold. Thresholds: 5% on converged F_x. On an unconverged mesh, only require HAMS to be no further from Nemoh than the Nemoh-HydroSTAR spread, and report the rest.
  6. Treat Nemoh as a code comparison, not as truth: it shares low-order panels and near-field error sources. Multi-body results are informational only, because Nemoh's documentation states its QTFs are untested for multi-body problems.

  Record versions, inputs and results in `docs/mean_drift.md`. WAMIT remains an option if a licence becomes available, and it is still needed to settle the conjugation and axes questions for WAMIT-format output.

- [ ] **Step 4: Document.** `docs/mean_drift.md` covers usage, methods, mesh guidance (finer panels near the waterline and corners), limits (generalized modes, narrow gaps), outputs and test results. Update the spec's D2 and D3 with the decision, and add a feature line to `readme.md`.

- [ ] **Step 5: Keep CI fast.** Check integration runtime stays acceptable. Move long cases out of CI.

- [ ] **Step 6: Final verification.** Run the full toolchain commands. Expected: all unit and integration tests pass, and `config1` outputs are unchanged.

- [ ] **Step 7: Commit** `docs: add mean drift guide and convergence results`.
