# Control-surface mean drift: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add Chen's (2007) middle-field mean drift as a third method (`drift_methods = control`). Phase A gives F_x, F_y and M_z; Phase B gives F_z, M_x and M_y. Results are per body, from a control surface C around each body.

**Architecture:**
- `ControlSurfaceMesh` builds the evaluation points per body: C panels, Γc segments, F panels.
- `FieldVelocity` gives the total φ and ∇φ at those points for every heading of a frequency, by the Green-formula differences or by sources.
- `ControlSurfaceDrift` turns two headings' fields into the six mean loads, with every product through `pair_minus`.
- The driver runs all of this in `so_end_frequency`, over the heading states it has kept.

**Tech Stack:** Fortran 2008, OpenMP, `fpm` with test-drive, Intel ifort/ifx with MKL. Python 3 for the fixtures (`pre_postprocessing_python/SecondOrder`).

**Spec:** `docs/superpowers/specs/2026-11-06-control-surface-drift-design.md`. Section numbers below (§) refer to it.

## Global Constraints

- **Branch.** Work on `feature/control-surface-drift`. Do not push, and do not open a PR unless the user asks.
- **Commits.** End every commit message with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Sources.** Do not read or copy Nemoh source (GPL). Cite equations from Chen (2007) and the WAMIT v7.5 manual only.
- **Snellius.** Full example runs go to Snellius, not the laptop. Before any Snellius job, give the user its size, node hours and memory, and wait for a yes.
- **Writing.** Simple, crisp language in comments and docs. Every claim is cited or justified, and each result is labelled measured or derived.
- **New files.** Every new `.f90` carries the Apache header of `src/SurfaceVelocity.f90`, author Vaibhav Raghavan.
- **Build.** `FPM_FFLAGS="-qmkl -qopenmp -heap-arrays 0"`, `FPM_LDFLAGS="-qmkl -qopenmp"`, `fpm build`, `fpm install --prefix install`.
- **Tests.**
  - `fpm test --target unit-tester` and `fpm test --target integration-tester -- install/bin/hamsmrel` must pass before each commit, apart from tests already failing on the base commit; record those in Task 1.
  - New unit suites that set module globals run serially, the way `waterline_field` does (see `test/unit_tester.f90`).
- **Build lists.** Add every new source file to `CompileHAMS.bat` and `HAMS_MREL.vfproj` by hand.
- **Conventions.**
  - Time factor e^{−iωt}.
  - η = iωφ/g_w at z = 0, with g_w = `wave_gravity(omega, k, h)`. This was checked: `incident_phi` gives η = A for the incident wave. NearFieldDrift's private `elevation(s, phi)` already does this, so reuse it.
  - The prefactors ρg use `s%rho` and `s%g`, as NearFieldDrift does.
  - Normals of C and F points point into the fluid domain D (§3.1). On C that is towards the body; on F it is downward.
  - n̄ on Γc and on Γ points into F.
  - Moments are about `mesh%xr(:, b)`.
- **Outputs.** Write through `write_drift_row` and `write_terms_row` (`SecondOrderOutput`), so the conjugation and scaling stay in one place.

## Review Focus

1. **Symmetric hull (`ISX`/`ISY` = 1).** C is a full surface and the images come from the field rows. A half mesh must give the full-mesh loads. Test: Task 8, `control_half_matches_full`.
2. **Generated C that would contain or cut another body** (a gap smaller than the margin). The run must stop with a message naming both bodies, not give a silent wrong load. Test: Task 3, `generated_c_rejects_other_body`.
3. **Zero and infinite frequency limits** (`is_limit_frequency`). There is no field there; the control step must skip them as far and near do. Test: Task 7, `control_skips_limits`.
4. **User C mesh with the normals the wrong way.** The run must stop with a message, not flip signs silently. Test: Task 4, `file_c_wrong_orientation_stops`.
5. **`mean_drift = 2`.** The diagonal pairs must equal the `mean_drift = 1` values, and F(β₂, β₁) must be the conjugate of F(β₁, β₂). Test: Task 8, `control_pairs_match_unidirectional`.

---

### Task 1: Bring in the finite-depth small-R Green function fix

The green route differences field rows, so G must be smooth in the field point (§4.4).

**Files:**
- Merge: branch `fix/finite-depth-small-r-green` (tip 0f92fec at the time of writing; take its current tip)

- [ ] **Step 1:** Run both test targets on the current commit and note any test that already fails (the baseline).
- [ ] **Step 2:** `git merge --no-ff fix/finite-depth-small-r-green`, then resolve conflicts, keeping both sides' tests.
- [ ] **Step 3:** Rebuild and run both test targets. Expected: no new failures against the baseline. The fix branch's `test_fingreen3d` cases pass.
- [ ] **Step 4:** Commit the merge (with the attribution line).

### Task 2: Options and input keys

**Files:**
- Modify: `src/SecondOrderTypes.f90` (`so_options`)
- Modify: `src/SecondOrderInput.f90` (`read_second_order` key cases, `drift_methods` parser, `validate_second_order`)
- Modify: `docs/input_files.md` (`SecondOrder.in` table and rejection list)
- Test: `test/test_second_order_input.f90`

**Interfaces:**
- Produces, in `so_options`:
  - `LOGICAL :: do_control = .FALSE.`
  - `INTEGER :: cs_velocity = 1` (1 green, 2 source)
  - `INTEGER :: cs_shape = 1` (1 cylinder, 2 box, 3 file)
  - `REAL(8) :: cs_margin = 0D0, cs_bottom_margin = 0D0, cs_panel_size = 0D0`
  - `REAL(8) :: cs_fd_step = 0.02D0`
  - `INTEGER :: cs_fs_rows = 0`
- Produces: `validate_second_order(opts, isol, h, isolv, irsp, has_gen_modes, ok, msg)`, signature unchanged, with these rules:
  - `do_control` with `has_gen_modes` gives `'drift_methods = control cannot be used with generalized modes.'`
  - `cs_velocity = 2` with `isolv == 2` gives `'cs_velocity = source needs the direct solver (ISOLV=1) in ControlFile.in.'`
  - `cs_velocity = 2` with `irsp == 1` gives `'cs_velocity = source cannot be used with irregular frequency removal (IRSP=1) in ControlFile.in.'`
  - `cs_margin`, `cs_bottom_margin` or `cs_panel_size` below 0, `cs_fd_step` outside (0, 0.5], or `cs_fs_rows` below 0: reject, with a message naming the key.

- [ ] **Step 1: Write the failing tests** in `test_second_order_input.f90`:
  - `control_keys_parse`: a file with `drift_methods = far control`, `cs_velocity = source`, `cs_shape = box`, `cs_margin = 12.5`, `cs_fd_step = 0.05` gives `do_control .and. do_far .and. .not. do_near`, `cs_velocity == 2`, `cs_shape == 2`, `abs(cs_margin - 12.5d0) < 1d-12` and `abs(cs_fd_step - 0.05d0) < 1d-12`.
  - `control_defaults`: `drift_methods = control` alone gives `cs_velocity == 1`, `cs_shape == 1`, `cs_margin == 0d0` and `cs_fd_step == 0.02d0`.
  - `control_rejections`: each rule above gives `ok == .false.` and its exact message.
- [ ] **Step 2:** Run `fpm test --target unit-tester`. Expected: the three tests fail.
- [ ] **Step 3:** Add the fields, the key cases (with `choose(value, 'green', 'source', ...)` and a three-way choice for `cs_shape`), `control` in the `drift_methods` parser, and the rules. Document the keys in `docs/input_files.md`, using §7's table text and the defaults from §5.1.
- [ ] **Step 4:** Run the tests. Expected: they pass.
- [ ] **Step 5:** Commit: `feat: read the control-surface drift options`.

### Task 3: Generated control surfaces and Γc

**Files:**
- Create: `src/ControlSurfaceMesh.f90` (module `ControlSurfaceMesh`)
- Modify: `src/SecondOrderTypes.f90` (type `cs_points`)
- Test: create `test/test_control_surface_mesh.f90`, and register it in `test/unit_tester.f90`

**Interfaces:**
- Produces, in `SecondOrderTypes`:

```fortran
TYPE :: cs_points
   ! Evaluation points of the control-surface method, all bodies together. w is the panel area or segment length.
   INTEGER :: n = 0
   REAL(8), ALLOCATABLE :: x(:,:), nrm(:,:), w(:), psize(:)   ! (3,n), (3,n), (n), (n)
   INTEGER, ALLOCATABLE :: kind(:), body(:)                   ! (n), (n)
END TYPE cs_points
INTEGER, PARAMETER :: CS_WALL = 1, CS_LINE = 2, CS_FREE = 3   ! C panel, Gamma_c segment midpoint, F panel
```

  On `CS_LINE` points, `nrm` holds n̄; on `CS_FREE` points it holds (0, 0, −1).
- Produces, in `ControlSurfaceMesh`:
  - `SUBROUTINE cylinder_surface(centre, radius, depth, panel_size, h, cm)`, with `REAL(8) :: centre(2)` and `TYPE(panel_mesh), INTENT(OUT) :: cm`.
  - `SUBROUTINE box_surface(centre, half, angle, depth, panel_size, h, cm)`, with `half(2)` the half-lengths and `angle` the rotation about z in radians.
  - Both give a full C (`nsys = 1`, `nbase = n`) with normals towards the body and the wall rows graded towards z = 0, top row about `panel_size/2`. When `h > 0` and `depth >= h`, the wall reaches z = −h and there is no bottom (§5.1).
  - `SUBROUTINE add_surface_points(cm, b, pts)` appends the C panels as `CS_WALL` points of body `b`, and the Γc segments as `CS_LINE` points. Γc is found by `build_surface_geom(cm, wl_z_tol, DEFAULT_SHARP_DEG, g, ok, msg)`, from `wl_xm`, `wl_len` and `wl_nh`.
  - `SUBROUTINE build_control_surface(mesh, geom, opts, h, inputdir, pts, cms, ok, msg)`, with `TYPE(panel_mesh), ALLOCATABLE, INTENT(OUT) :: cms(:)` (one per body). It dispatches on `opts%cs_shape` and applies the defaults of §5.1:
    - margins of `4*median(psize of body b)`;
    - a panel size of `median(psize of body b)`;
    - the cylinder centred on the centre of the body's waterplane bounding box, with radius the largest horizontal distance from that centre to a hull node plus `cs_margin`;
    - the box aligned with the body's `LCS` angle, with half-lengths the hull extent plus `cs_margin`;
    - depth equal to the deepest hull node plus `cs_bottom_margin`.
  - `SUBROUTINE check_control_surface(mesh, pts, cms, ok, msg)`:
    - for a generated C, every node of another body must be outside body b's C, by an exact cylinder or box test;
    - for any C, no C point may lie within one local panel size of a hull panel centroid of any body.
  - `SUBROUTINE write_control_surface(outputdir, cms)` writes `CSMesh_{b}.pnl` in the `HullMesh.pnl` format, which the Task 4 reader can read back.

- [ ] **Step 1: Write the failing tests:**
  - `cylinder_area_and_normals`: `cylinder_surface([0d0,0d0], 20d0, 15d0, 2d0, 100d0, cm)`.
    - Σ area of the wall panels is within 0.5% of 2π·20·15, and Σ area of the bottom panels within 0.5% of π·20².
    - Every wall normal satisfies `dot(nrm, x_c) < 0` (towards the axis), and every bottom normal has `nz > 0.99`.
  - `cylinder_to_seabed`: with `h = 10d0` and depth 12, no panel has `nz > 0.5`, and the deepest node is at z = −10 within 1e-12.
  - `box_area_and_normals`: `box_surface([5d0,-3d0], [12d0,7d0], 0.3d0, 8d0, 1d0, 50d0, cm)`. The wall area is within 0.5% of 2·(24+14)·8, and the bottom area within 0.5% of 24·14. The normals point towards the centre.
  - `gamma_c_line`: after `add_surface_points` on the cylinder:
    - the `CS_LINE` lengths sum to within 0.5% of 2π·20;
    - every `CS_LINE` point has z = 0 within 1e-12;
    - every `nrm` is horizontal, of unit length, with `dot(nrm, x) < 0`.
    - If `wl_nh` comes out pointing away from the axis, negate it in `add_surface_points` and say so in a comment. Do not change `build_surface_geom`.
  - `generated_c_rejects_other_body` (Review Focus 2): two `make_cylinder` hulls of radius 5, 12 m apart, with `cs_margin = 10`. `check_control_surface` returns `ok = .false.`, and `msg` contains `'body 1'` and `'body 2'`.
- [ ] **Step 2:** Run the unit tests. Expected: the new suite fails to compile or fails.
- [ ] **Step 3:** Implement the routines and add `ControlSurfaceMesh.f90` to the build lists.
- [ ] **Step 4:** Run the tests. Expected: they pass.
- [ ] **Step 5:** Commit: `feat: generate control surfaces around each body`.

### Task 4: Control surface from a user file

**Files:**
- Modify: `src/ControlSurfaceMesh.f90`
- Modify: `docs/input_files.md` (`CSMesh_{b}.pnl`)
- Test: `test/test_control_surface_mesh.f90`

**Interfaces:**
- Produces: `SUBROUTINE read_control_mesh(path, cm, ok, msg)`. It reads the `HullMesh.pnl` format, and requires ISX = ISY = 0 in the file (C is always full, §2).
- `build_control_surface` with `cs_shape = 3` reads `inputdir//'CSMesh_'//b//'.pnl'` (`b` without leading zeros) for every body. A missing file stops the run with `'cs_shape = file needs CSMesh_<b>.pnl in the input directory.'`
- Orientation check: S = Σ dot(nrm_i, xc_i − x_b)·area_i, with x_b the centre of the body's waterplane bounding box at z = 0. S must be negative. Otherwise `ok = .false.`, with `msg = 'CSMesh_<b>.pnl: the panel normals must point towards the body.'`

- [ ] **Step 1: Write the failing tests:**
  - `file_c_round_trip`: write a generated cylinder with `write_control_surface`, then read it back. The node coordinates and the panel areas agree to 1e-9.
  - `file_c_wrong_orientation_stops` (Review Focus 4): the same file with every panel's node order reversed gives `ok = .false.` and the exact message.
- [ ] **Step 2:** Run them. Expected: they fail.
- [ ] **Step 3:** Implement the reader and the check. Document the file in `docs/input_files.md`.
- [ ] **Step 4:** Run them. Expected: they pass.
- [ ] **Step 5:** Commit: `feat: read a control surface from CSMesh files`.

### Task 5: Field φ and ∇φ by Green-formula differences

**Files:**
- Modify: `src/WaterlineField.f90` (`field_influence_row`: optional `x_decide`)
- Create: `src/FieldVelocity.f90` (module `FieldVelocity`)
- Modify: `test/second_order_fixtures.f90` (exact MacCamy-Fuchs field)
- Test: create `test/test_field_velocity.f90` (serial suite), and register it

**Interfaces:**
- Consumes: `cs_points` (Task 3); `field_influence_row`, `field_potential` (WaterlineField); `incident_phi`, `incident_grad`, `wave_gravity` (FirstOrderState).
- Produces:
  - `field_influence_row(mesh, xp, is_multi, row_g, row_dg, x_decide)`, with `REAL(8), INTENT(IN), OPTIONAL :: x_decide(3)`. When present, the near-or-far flag of each panel and image is computed from `image_signs(is)*x_decide` instead of `xf` (§4.1). Without it, behaviour is unchanged.
  - In `SecondOrderTypes`: `TYPE :: cs_fields` with `COMPLEX(8), ALLOCATABLE :: phi(:), grad(:,:)` of size (n) and (3,n), the total potential and its gradient at the points.
  - `SUBROUTINE fields_green(mesh, pts, states, step, is_multi, fields)`, with `TYPE(first_order_state), INTENT(IN) :: states(:)` and `TYPE(cs_fields), INTENT(OUT) :: fields(SIZE(states))`.
    - Each point uses δ = `step*pts%psize(i)`, with every stencil point's flag decided at `pts%x(:, i)`.
    - Points with z < −δ: 6 points, central differences.
    - Points with z ≥ −δ: the 4 horizontal points, plus the centre for φ, with ∂φ/∂z = `s%nu` × φ_total.
    - φ_total is computed at the centre for `CS_LINE` and `CS_FREE`; `CS_WALL` stores φ only if it is computed anyway.
    - The incident part comes from `incident_phi` and `incident_grad`, with g_w as in `build_state`.
    - The rows of one point are applied to every state before the next point. No rows are stored (§4.3).
- Produces, in fixtures: `subroutine maccamy_fuchs_field(x, a_rad, h, k, omega, amp, beta, g, phi, grad)`, the exact total potential and gradient at any r ≥ a of a fixed bottom-mounted cylinder. It uses φ = −(igA/ω)·Z(z)·Σ ε_n iⁿ [J_n(kr) − (J_n′(ka)/H_n′(ka)) H_n(kr)] cos n(θ−β) (MacCamy and Fuchs 1954), with the depth factor of `depth_factor`. It agrees with `maccamy_fuchs_phi` at r = a.

- [ ] **Step 1: Write the failing tests:**
  - `field_rows_match_analytic_gradient`: on the `make_cylinder(5d0, 5d0, 12, 2, 2, ...)` hull at h = 10, take points 1, 3 and 8 panel sizes off a panel. Compare the central difference of `field_influence_row`'s `row_g` with step 0.02·psize (decided at the centre) against `source_panel_integrals`' `grad_int` for the same panel and image. Every entry must agree to 1e-4 relative to |grad_int|, and also in deep water (h = −1).
  - `incident_stencil`: with `phi_p = dphi_p_dn = 0` (incident only), `fields_green` at 20 points (some at z = 0, some below) reproduces `incident_grad` to 1e-5 relative, and `incident_phi` to 1e-12.
  - `maccamy_fuchs_field_consistent`: `maccamy_fuchs_field` at r = a equals `maccamy_fuchs_phi`. Its gradient matches a central difference of itself with step 1e-5·a to 1e-6.
  - `green_field_matches_exact`: `maccamy_fuchs_state` on the 64 × 16 bottom-mounted mesh of `drift_field_*` (a = 5, h = 10), with `fields_green` at 40 points on r = 7.5 (some at z = 0). φ must agree with `maccamy_fuchs_field` to 1% of max|φ|, and ∇φ to 2% of max|∇φ|. Record the measured values in the test comment; the limits are those of the test.
- [ ] **Step 2:** Run the tests. Expected: they fail.
- [ ] **Step 3:** Implement `x_decide`, `cs_fields` and `fields_green`, and add `FieldVelocity.f90` to the build lists.
- [ ] **Step 4:** Run the tests. Expected: they pass. The existing `waterline_field` suite is unchanged.
- [ ] **Step 5:** Commit: `feat: field potential and velocity by Green-formula differences`.

### Task 6: Field φ and ∇φ by the source route

**Files:**
- Modify: `src/SurfaceVelocity.f90` (`velocity_source`: optional `sigma_full`)
- Modify: `src/FieldVelocity.f90`
- Test: `test/test_field_velocity.f90`

**Interfaces:**
- Produces: `velocity_source(mesh, op, s, phi_source, sigma_full)`, with `COMPLEX(8), INTENT(OUT), OPTIONAL :: sigma_full(:)` of size (n). It is σ on every panel, rebuilt from the blocks as `grad_p` is: `sigma_full((im-1)*nb+i) = SUM(RXY(im,1:ns)*sigma(i,1:ns))/ns`.
- Produces: `SUBROUTINE fields_source(mesh, pts, states, sigmas, is_multi, fields)`, with `COMPLEX(8), INTENT(IN) :: sigmas(:,:)` of size (n, nstate). For each point, base panel j and image is:
  - `xf = image_signs(is)*x`, then `source_panel_integrals(mesh, xf, j, flag, is_multi, g_int, grad_int)`;
  - φ_P += σ·g_int/(4π) and ∇φ_P += `image_signs(is)`·grad_int·σ/(4π), with σ = `sigmas((is-1)*nb+j, :)`;
  - flag = 1 within 50 panel sizes, as in `build_source_operator`;
  - the incident part as in Task 5.

- [ ] **Step 1: Write the failing tests:**
  - `sigma_full_reproduces_velocity`: on the source fixtures' hull with `ISY = 1`, Σ over panels of σ·∫∇G at a panel centroid, using `sigma_full`, gives `grad_p` from `velocity_source` to 1e-10.
  - `source_field_matches_exact`: as `green_field_matches_exact`, with `fields_source`, φ to 2% and ∇φ to 3%. Record the measured values.
  - The comparison of green against source on a BEM-solved floating hull needs the driver. It is `control_source_far` in Task 8.
- [ ] **Step 2:** Run them. Expected: they fail.
- [ ] **Step 3:** Implement `sigma_full` and `fields_source`.
- [ ] **Step 4:** Run them. Expected: they pass. The `surface_velocity` and `source_velocity` suites are unchanged.
- [ ] **Step 5:** Commit: `feat: field potential and velocity from the source distribution`.

### Task 7: Phase A loads, driver and output

**Files:**
- Create: `src/ControlSurfaceDrift.f90` (module `ControlSurfaceDrift`)
- Modify: `src/NearFieldDrift.f90` (make `elevation` and `force_and_moment` public)
- Modify: `src/SecondOrderDriver.f90`, `src/SecondOrderOutput.f90` (only if a new writer is needed; none is expected)
- Modify: `docs/output_files.md`
- Test: create `test/test_control_drift.f90`, and register it; plus `test/test_integration_drift.f90`

**Interfaces:**
- Consumes: `cs_points`, `cs_fields`, `pair_minus`, `elevation(s, phi)`, `force_and_moment(x, v) -> COMPLEX(8)(6)`.
- Produces:
  - `INTEGER, PARAMETER :: CS_T_LINE = 1, CS_T_FLUX = 2, CS_T_PRESSURE = 3, CS_T_WATERLINE = 4, CS_T_FREE = 5, CS_T_HYDROSTATIC = 6`
  - `SUBROUTINE control_surface_pair(pts, xr, s1, s2, f1, f2, terms)`, with `REAL(8) :: xr(:,:)` (3, nbody) and `COMPLEX(8), INTENT(OUT) :: terms(6, 6, nbody)` (component, term, body). Phase A fills three terms, as written by Chen eq. (37) and (38) for the C and Γc parts (§3.2):
    - `CS_T_LINE`: (ρg/2) Σ_{CS_LINE} pair_minus(η₁, η₁, η₂, η₂)·w·force_and_moment(x − xr, n̄);
    - `CS_T_FLUX`: ρ Σ_{CS_WALL} w·force_and_moment(x − xr, P(Φ_n, ∇Φ)), where P(Φ_n, ∇Φ)_i = pair_minus(n·g₁, g₁ᵢ, n·g₂, g₂ᵢ);
    - `CS_T_PRESSURE`: −(ρ/2) Σ_{CS_WALL} w·Σᵢ pair_minus(g₁ᵢ, g₁ᵢ, g₂ᵢ, g₂ᵢ)·force_and_moment(x − xr, n).
  - Driver behaviour:
    - When `do_control` is set, `states` holds every heading of the frequency (as for `mean_drift = 2`).
    - `build_control_surface` and `check_control_surface` run once in `common_setup`, and `write_control_surface` writes to `outputdir`.
    - In `so_end_frequency` (skipped at limit frequencies), `fields_green` or `fields_source` runs once for all headings. For `cs_velocity = 2`, the σ of each heading is kept from `after_motion`.
    - Then `control_surface_pair` runs for each pair: (ib, ib) for `mean_drift = 1`, all pairs for `mean_drift = 2`.
  - Output:
    - `Wamit_format/MeanDriftControlSurface.7`, rows I = 1, 2, 6 per body, index 6(b−1)+I. Moments go through `write_drift_row(..., i > 3)`.
    - With `write_terms >= 1`, rows in `MeanDriftTerms.txt` with term index `10 + CS_T_*`.
    - The log: the point counts, the field time per frequency, and with far-field the differences F_x, F_y, M_z against far, with M_z moved to the origin as `near_field_output` does.

- [ ] **Step 1: Write the failing unit tests** in `test_control_drift.f90`, with fields built directly, no BEM:
  - `incident_only_zero`: incident fields on a generated cylinder C (radius 30, depth 30, h = 100, ω = 0.8, 1.2, 2.0, headings 0 and 0.7). F_x, F_y and M_z must be within 1e-3 ρgA²·R of zero (R the C radius). The error must fall by at least 3× when `panel_size` is halved.
  - `incident_only_zero_box`: the same on a box C, at heading 0.7.
  - `maccamy_fuchs_exact`: fields from `maccamy_fuchs_field` on a cylinder C of radius 8 around a = 5, h = 10, reaching the seabed, heading 0, k from ω = 1.0 and 1.6. F_x must match `maccamy_fuchs_drift(a, k, kappa, 720)`, with kappa = (ρg/2)·`group_ratio(k, h)`, to 0.5%, and |F_y|, |M_z| must be below 1e-6 of it.
  - `maccamy_fuchs_independent_of_c`: the same with radius 6, 8 and 15, and with a box of half-lengths 7 × 9. F_x agree with each other to 0.5%.
  - `pair_symmetry` (Review Focus 5, unit level): `control_surface_pair(pts, xr, s1, s2, f1, f2)` is the conjugate of `(pts, xr, s2, s1, f2, f1)`.
- [ ] **Step 2:** Run them. Expected: they fail.
- [ ] **Step 3:** Implement `control_surface_pair` and make the two NearFieldDrift helpers public.
- [ ] **Step 4:** Run them. Expected: they pass.
- [ ] **Step 5: Write the failing integration tests.** Add the fixture `control_cyl`: the hull of `drift_cyl_fixed` (448 panels), with `drift_methods = far control`, `motion_state = fixed`, 0.70, 0.95 and 1.10 rad/s, headings 0 and 30. Generate it with `make_drift_fixtures.py`, and document it in `test/data/drift_cyl/README.md`. Tests:
  - `control_runs`: `MeanDriftControlSurface.7` exists, with 3 frequencies × 2 headings × 3 rows, and `CSMesh_1.pnl` is written.
  - `control_skips_limits` (Review Focus 3): `control_cyl_sybo`, as `drift_cyl_sybo` with `control`, has no rows at the limit frequencies.
- [ ] **Step 6:** Implement the driver and output wiring. Document the file in `docs/output_files.md`.
- [ ] **Step 7:** Run both test targets. Expected: they pass, and every older test is unchanged.
- [ ] **Step 8:** Commit: `feat: control-surface drift for F_x, F_y and M_z`.

### Task 8: Phase A physics tests on the laptop fixtures

**Files:**
- Modify: `pre_postprocessing_python/SecondOrder/make_drift_fixtures.py`, `test/data/drift_cyl/README.md`
- Create: fixtures under `test/data/control_*`
- Test: `test/test_integration_drift.f90`

Each fixture reuses an existing hull and `ControlFile.in`, with `drift_methods = far control` (and `near` where noted).

| Fixture | From | Test and limit |
| :--- | :--- | :--- |
| `control_cyl`, `control_cyl_free` | `drift_cyl_fixed`, `drift_cyl` (sharp) | `control_far_sharp`: F_x against far, at each frequency and heading. Record the measured value; the limit is 2× that, capped at 3%. If the measured value exceeds 3% at 448 panels, stop and report: that contradicts §1 |
| `control_cyl_wide` | `control_cyl` with `cs_margin = 30` | `control_independent_of_c`: F_x matches `control_cyl` to 1% |
| `control_cyl_box` | `control_cyl` with `cs_shape = box` | `control_cylinder_box`: F_x matches `control_cyl` to 1% |
| `control_bottom` | `drift_field_0p01` (bottom-mounted, a = 5, h = 10) | `control_maccamy_fuchs`: F_x against `maccamy_fuchs_drift`. Record the measured value; the limit is 2× that, capped at 2%. F_y and M_z below 1e-4 of F_x |
| `control_box` | `drift_box` (sharp, fixed, heading 30) | `control_far_box`: F_x, F_y and M_z against far. Record; 2× measured, capped at 3% |
| `control_two_round` | `drift_two_round`, plus `near` | `control_two_body_total`: the sum over bodies of F_x and F_y against far, 2%. Each body's F_x and F_y against near, 2% |
| `control_far_apart` | `drift_far_apart` | `control_far_apart_isolated`: each body against `control_isolated` (`drift_isolated` with `control`), 0.5% |
| `control_irsp`, `control_irsp_ref` | `drift_irsp`, `drift_irsp_ref` | `control_irsp_stable`: at 0.95 and 1.10, IRSP = 1 against IRSP = 0 to 0.5%. At 1.548, control (IRSP = 1) against far (IRSP = 1) to 3% |
| `control_half` | `drift_cyl_half` | `control_half_matches_full` (Review Focus 1): matches `control_cyl_free` to 1e-6 relative |
| `control_gmres` | `drift_cyl_gmres` | `control_gmres_agrees`: matches `control_cyl_free` to 1e-5 relative |
| `control_pairs` | `drift_cyl_sybo` hull, `mean_drift = 2` | `control_pairs_match_unidirectional` (Review Focus 5): diagonal pairs equal `control_cyl_free` to 1e-10; F(β₂, β₁) = conj F(β₁, β₂) |
| `control_source` | `control_cyl_free` with `cs_velocity = source` | `control_source_far`: against far. Record; 2× measured, capped at 3% |
| `control_step_*` | `control_cyl` with `cs_fd_step` 0.005, 0.01, 0.02, 0.05, 0.1 | `control_step_plateau`: the spread over 0.01 to 0.05 below 0.2% of F_x. Record all five in the README. If the plateau is not centred on 0.02, change the default to its centre, and change `docs/input_files.md` and the spec §4.1 to match |

- [ ] **Step 1:** Add the fixtures to `make_drift_fixtures.py`, generate them, and commit the inputs.
- [ ] **Step 2:** Write the tests in the table, with placeholders only for the measured values.
- [ ] **Step 3:** Run the integration tests once, and record every measured value in `test/data/drift_cyl/README.md`. Set each limit by the rule in its row.
- [ ] **Step 4:** Run them again. Expected: all pass. If a "stop and report" rule fires, stop and report to the user, with the numbers.
- [ ] **Step 5:** Commit: `test: control-surface drift against far-field, exact and itself`.

### Task 9: Phase A convergence study on Snellius

**Files:**
- Modify: `pre_postprocessing_python/SecondOrder/drift_convergence.py` (a `control` method column, green and source)
- Results: `../validation-runs/snellius-drift-control/` (outside the repo)

- [ ] **Step 1:** Extend the driver to add `control` to the existing cases:
  - sharp and rounded cylinders, fixed and free, m1 to m4, 11 frequencies;
  - the boxes at heading 30;
  - the bottom-mounted cylinder.
- [ ] **Step 2:** Estimate the cost from the run-time table of `../validation-runs/snellius-drift/tables.md`, scaled by the Green call count of spec §9. Report to the user the number of jobs, the node hours, the peak memory per job and the wall time. **Wait for a yes.**
- [ ] **Step 3:** Submit, collect, and write `tables.md` in the new results folder with the layout of the old one.
- [ ] **Step 4:** Check the Phase A gate (spec §8.2): sharp cylinder and sharp box at m2, within 1%, or within 0.001 ρgA²a where far-field is small. Report the result to the user.
- [ ] **Step 5:** Put a summary table and the gate result in `test/data/drift_cyl/README.md`. Commit: `docs: record the control-surface convergence study`.

**Phase B starts only after the Phase A gate passes, and after the user says to go on.**

### Task 10: Free-surface patch F

**Files:**
- Modify: `src/ControlSurfaceMesh.f90`, `docs/input_files.md`
- Test: `test/test_control_surface_mesh.f90`

**Interfaces:**
- Produces:
  - `SUBROUTINE free_surface_strips(mesh, geom, b, pts, nrows, fm, ok, msg)`, with `TYPE(panel_mesh), INTENT(OUT) :: fm`.
    - Quads go from each hull waterline segment of body b out to Γc of body b, along rays from the centre x_b, with `nrows` rows graded towards the hull (the inner row 1/4 of the outer).
    - The outer nodes are the intersections of the rays with the `CS_LINE` segments of body b.
    - When the waterline is not star-shaped about x_b (the polar angle of the segment end points is not monotonic over one turn), `ok = .false.`, with `msg = 'body <b>: the waterline is not star-shaped about its centre; give FSMesh_<b>.pnl.'`
  - `FSMesh_{b}.pnl` is read with `read_control_mesh` when present, and takes precedence over the strips.
  - `build_control_surface` appends the F panels as `CS_FREE` points with `nrm = (0, 0, -1)`, only when Phase B output is on (Task 11).
  - When C reaches the seabed, F is still built. The seabed patch is not built: Task 11 then writes no rows 3 to 5 for that body (see the note there).
  - The default for `cs_fs_rows` is the nearest integer to (Γc distance − waterline distance)/`cs_panel_size`, at least 2.

- [ ] **Step 1: Write the failing tests:**
  - `strips_area`: the `make_cylinder(5d0, ...)` hull with a C of radius 15. Σ F area is within 0.5% of π(15² − 5²), and every F centroid lies in 5 < r < 15 at z = 0.
  - `strips_box_hull`: a `drift_box`-like hull (20 × 10) with a cylinder C of radius 25. The area is within 0.5% of π·25² − 200.
  - `not_star_shaped_stops`: an L-shaped waterline gives `ok = .false.` and the exact message.
- [ ] **Step 2:** Run them. Expected: they fail.
- [ ] **Step 3:** Implement the strips and the file path.
- [ ] **Step 4:** Run them. Expected: they pass.
- [ ] **Step 5:** Commit: `feat: free-surface patch between the waterline and the control surface`.

### Task 11: Phase B loads

**Files:**
- Modify: `src/ControlSurfaceDrift.f90`, `src/SecondOrderDriver.f90`, `docs/output_files.md`
- Test: `test/test_control_drift.f90`, `test/test_integration_drift.f90`

**Interfaces:**
- Consumes: `elevation`, `vertical_motion` and `waterline_point` from NearFieldDrift. The hull waterline φ comes from `waterline_extrapolate` or `waterline_field` (as chosen by `waterline_method`, with the hull velocity from `velocity_method` when it is extrapolated). The hydrostatic term comes from `add_hydrostatic_term`, made public, but only after Step 1 below.
- Produces: `control_surface_pair` gains optional arguments `(mesh, geom, phi_wl1, phi_wl2)` and fills three more terms (§3.3):
  - `CS_T_WATERLINE`: −(ρg/2) Σ_Γ ⟨Ξ⟩·wl_len·force_and_moment(x − xr, k), with Ξ = (η² − 2ηζ₃)·n₃/√(1−n₃²) − 2η(X·n̄);
    - n₃ = −`wl_nz` and n̄ = −`wl_nh`, since Chen's normals point into the fluid and NearFieldDrift's out of it;
    - X is the body displacement at the waterline point, T + R∧r. Add `displacement(xi, x) -> COMPLEX(8)(3)` to NearFieldDrift next to `vertical_motion` if it is not there.
  - `CS_T_FREE`: ρ Σ_{CS_FREE} w·[pair_minus(η₁, −iων₁φ₁, η₂, −iων₂φ₂) + ½ Σᵢ pair_minus(g₁ᵢ, g₁ᵢ, g₂ᵢ, g₂ᵢ)]·force_and_moment(x − xr, k);
  - `CS_T_HYDROSTATIC`: from `add_hydrostatic_term`.
- Output: rows I = 3, 4, 5 in `MeanDriftControlSurface.7`, for bodies whose C has a bottom. For a body whose C reaches the seabed, the rows are not written, and the log says `'control surface of body <b> reaches the sea bed: F_z, M_x, M_y not computed.'` This narrows spec §5.1, which also planned a seabed patch: that needs the hull's seabed contact line, which no current routine gives, and no Phase B validation case needs it. Put the narrowing in the spec in the same commit.

- [ ] **Step 1: Check the hydrostatic term (spec §3.3 open item).** Read Molin and Hairault (1983) if the user can supply it. Otherwise compare in algebra: Chen's near-field (17)–(18) waterline term η(2ζ₃ − η)n′ plus (F_S, M_S), against NearFieldDrift's terms I, Ia and V (WAMIT 15.54–15.56).
  - Write the result in the spec §3.3 as derived, or as unresolved.
  - If unresolved, go on, with the integration gate below as the test, and say so in the commit message.
- [ ] **Step 2: Write the failing unit test** `incident_only_vertical_zero`. Incident fields, with F covering the whole disc inside Γc (strips from the centre: build the points directly in the test) and no hull. F_z, M_x and M_y must be within 1e-3 ρgA²R of zero, falling by at least 3× on halving.
- [ ] **Step 3: Write the failing integration tests:**
  - `control_vertical_round`: the fixtures `control_round` and `control_round_fixed` (`drift_round` hulls with `near control`). F_z, M_x and M_y from control against near. Record; the limit is 2× measured, capped at 3%.
  - `control_vertical_independent_of_c`: margins 4 and 8 panel sizes agree to 1%.
  - `control_seabed_rows_absent`: `control_bottom` has no rows 3 to 5, and has the log line.
- [ ] **Step 4:** Implement it. Run both test targets. Expected: they pass. If `control_vertical_round` is beyond 3%, stop and report, with the per-term breakdown from `MeanDriftTerms.txt`. Do not change formulas to fit.
- [ ] **Step 5:** Commit: `feat: control-surface drift for F_z, M_x and M_y`.

### Task 12: Phase B study on Snellius and wrap-up

- [ ] **Step 1:** Add rows 3 to 5 to the study collector. Estimate the cost as in Task 9 (only the rounded and sharp hulls, m1 to m4), report it to the user, and **wait for a yes.**
- [ ] **Step 2:** Run, collect, and check the Phase B gate (spec §8.3): within 1% of near-field on the rounded hulls at m3. Report it.
- [ ] **Step 3:**
  - Update `docs/input_files.md`, `docs/output_files.md` and `test/data/drift_cyl/README.md` with the final defaults and measured accuracy.
  - Add a short release note in the place the fix branches use.
  - Commit: `docs: record the vertical control-surface results`.
- [ ] **Step 4:** Use superpowers:finishing-a-development-branch. Do not push or open a PR without the user's word.
