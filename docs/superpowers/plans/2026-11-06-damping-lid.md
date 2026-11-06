# Damping Lid Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A lid of any 2D planform on z = 0 on which the total potential is zero, for single and multiple bodies, with no change to runs without a lid.

**Architecture:** The lid adds one unknown per panel, q = dphi/dz - nu phi, and one row per panel: the `IRSP` row (zero diagonal, field point on z = 0) plus a single-layer column for q. The new code lives in `src/DampingLid.f90` and `src/DampingLidMulti.f90`, with lid data in a new `Lid_mod`. The existing no-lid and `IRSP` paths gain only branches that are off unless `Damping_lid 1` is set.

**Tech Stack:** Fortran 90 (ifort 2022.1), MKL LAPACK (`ZGETRF`/`ZGETRS`, `ZGELS`), fpm with test-drive, Python 3.7 + numpy 1.21 + scipy 1.7 for fixtures and the independent reference.

**Spec:** `docs/superpowers/specs/2026-11-05-damping-lid-design.md`

## Global Constraints

- Branch `feature/damping-lid` from `origin/main` 6f4cc3b. Never push or open a PR unless the user asks.
- Runs without `Damping_lid 1` give unchanged results. Never regenerate stored expected outputs.
- Build and test only with ifort through `cmd //c "C:\Users\vraghavan1\Downloads\MREL_HAMS-MREL\validation-runs\tools\fpmenv_baseline.bat" C:\Users\vraghavan1\Downloads\MREL_HAMS-MREL\HAMS-MREL-damping-lid <fpm args>`. ifx cannot compile test-drive.
- No full example cases on the laptop. Larger runs go to Snellius by key login only (`ssh -o BatchMode=yes -o PasswordAuthentication=no vraghavan@snellius.surf.nl`), never the password file, with job size and cost shown to the user first.
- Do not read Nemoh source.
- Commit trailer: `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Commit only files the task changed.
- New lid arrays are `REAL*8` / `COMPLEX*16`, not `GRN_KIND`.
- Keep the stop in `CalNormals` for hull panels with centroid on z = 0.
- Every fixtures README records the measured values and the limits set from them. A claim in docs separates what was measured from what was derived.

## Review Focus

1. A lid touching the hull waterline (annulus starting at the hull radius) must give finite, convergent results. Test: `lid_touch` in Task 4.
2. A lid in two separate patches must work like any other mesh. Test: `lid_patches` in Task 4.
3. Triangle lid panels must agree with quads on the same planform. Tests: `lid_rankine_triangle` in Task 2 and `lid_ann_tri` in Task 4.
4. Zero and infinite frequency in a lid run must not give NaN, and at infinite frequency the result must equal the no-lid result. Test: `lid_limits` in Task 4.
5. Field points on a lid node, on an edge, or inside a panel off its centroid must give finite elevation near 0. Test: `lid_field_off_centroid` in Task 5.

---

## Commands used by every task

- Build: `<env> build`, then `<env> install --prefix install`
- Unit tests: `<env> test --target unit-tester`
- Integration tests: `<env> test --target integration-tester -- install/bin/hamsmrel`
- Fixtures: `python -I pre_postprocessing_python/DampingLid/make_lid_fixtures.py`, which writes `test/data/damping_lid/<case>/input/`. Each task adds its cases to this script.

`<env>` is the `fpmenv_baseline.bat` call from Global Constraints. Expected for every build: exit code 0, no warnings in new files. Expected after every task: the full unit and integration suites pass, including all pre-existing tests.

### Task 0: Baseline

**Files:** none.

- [ ] **Step 1:** Build, install, and run both suites on the unchanged branch. Record pass counts and wall times in the scratchpad. Expected: all pass. If any fails, stop and report: later tasks rely on this baseline.

### Task 1: Lid input, geometry and stops

**Files:**
- Modify: `src/WavDynMods.f90` (new `MODULE Lid_mod` after `Inerfs_mod`)
- Modify: `src/InputFiles.f90:177-201` (`ReadControlFile`), `src/InputFiles.f90:236-282` (`VerifyInputFilesExist`)
- Create: `src/DampingLid.f90`
- Modify: `src/HAMS_Prog.f90` (read and check after `CalNormals` / `CalNormalsMulti`; report counts)
- Modify: `CompileHAMS.bat`, `HAMS_MREL.vfproj` (add `DampingLid.f90`, and later `DampingLidMulti.f90`, after `PressureElevationMulti.f90` and before `HAMS_Prog.f90`)
- Create: `pre_postprocessing_python/DampingLid/make_lid_fixtures.py`, `test/data/damping_lid/README.md`, `test/test_integration_lid.f90`
- Modify: `test/integration_tester.f90`, `docs/input_files.md`

**Interfaces:**
- Produces, in `Lid_mod`: `INTEGER LIDON` (0 or 1), `LNELEM`, `LNTND`, `LNCN(:)`, `LNCON(:,4)`, `REAL*8 LXYZ(:,3)`, `LXYZ_P(:,3)`, `LDS(:)`, `LPNSZ(:)`, `COMPLEX*16 MXQ(:,:,:)`. `MXQ` is allocated `(LNELEM,7,NSYS)` for one body and `(LNELEM,6*NBODY+1,NSYS)` for several.
- Produces, in `DampingLid`:
  - `SUBROUTINE READ_LID_MESH(INPUTDIR)` reads `DampingLid.pnl` in the `WaterPlaneMesh.pnl` format, without touching `ISX`/`ISY`.
  - `SUBROUTINE LID_GEOMETRY()` fills `LXYZ_P`, `LDS` and `LPNSZ` with `CalPanelCentre`, `CalPanelArea` and `CalPanelChartLength`.
  - `SUBROUTINE CHECK_LID(NBODYIN)` runs every stop below.
  - `REAL*8 FUNCTION HULL_SOLID_ANGLE(XF)` returns the sum over hull panels, all symmetry images, of the double-layer Rankine integral plus its image at `XF` on z = 0. Single body: `SGLINTBD_TRI2`/`QUAD2` with `IRR=1` and `GEONRM_P`; multi: the `_MULTI` routines on the combined arrays.

`Damping_lid` sits on the line after `Solver_type`, or in its place when `Solver_type` is absent. It is read with the same backspace-on-miss pattern and the integer parsed after the keyword. Any value other than 0 or 1 warns and uses 0.

Every stop prints its message and then executes `STOP 1`:

| Condition | Message contains |
| --- | --- |
| `LIDON = 1`, no `DampingLid.pnl` | `DampingLid.pnl input file missing` |
| node with abs(z) > 1e-10 | `lid node is not on z = 0` |
| lid `ISX`/`ISY` differ from hull | `lid symmetry differs from the hull` |
| `ISOLV = 2` with lid | `damping lid needs Solver_type 1` |
| abs(`HULL_SOLID_ANGLE`) > 2 pi at a lid centroid | `lid panel inside a hull` |

- [ ] **Step 1: Fixtures.** Write `make_lid_fixtures.py`. It copies hull, hydrostatics and control file from `test/data/image_term/image_cos/input` (cylinder r = 5, draft 5, depth 50, 1.35 rad/s), sets `Damping_lid 1`, and writes the lid from a function `annulus(r0, r1, nt, nr, tri=False, reverse=False)`. Cases for this task:
  - `lid_stop_missing` (no file)
  - `lid_stop_offplane` (one node at z = 1e-3)
  - `lid_stop_symmetry` (lid header ISY = 1, hull 0)
  - `lid_stop_gmres` (`Solver_type 2`)
  - `lid_stop_inside` (disc r < 3)
  - `lid_flag_off` (lid file present, `Damping_lid 0`)

- [ ] **Step 2: Failing tests** in `test/test_integration_lid.f90`, module `test_integration_lid` with `setup_lid` and `collect_tests_lid`:
  - `stops_lid_missing`, `stops_lid_offplane`, `stops_lid_symmetry`, `stops_lid_gmres`, `stops_lid_inside`: exit code non-zero and the message from the table in the log. Copy the pattern of `run_case_expect_stop` from `test_integration_image.f90`, with the message as an argument.
  - `flag_off_matches_no_lid`: `AmssDamp.1` and `ExcForce.3` of `lid_flag_off` equal those of `image_cos` exactly. Reuse its output by running `image_cos` here as well.
  - Register both suites in `integration_tester.f90`. Run: the stop tests fail (the run does not stop), and `flag_off_matches_no_lid` fails to read the new line or passes.
- [ ] **Step 3:** Implement `Lid_mod`, the control-file read, the file check, `READ_LID_MESH`, `LID_GEOMETRY`, `HULL_SOLID_ANGLE` and `CHECK_LID`. Call them in `HAMS_Prog.f90` right after `CalNormals` / `CalNormalsMulti`, and print ` Number of panels on the damping lid:` with `LNELEM` next to the existing panel counts and in `ErrorCheck.txt` (`CreateErrorCheckFile`, only when `LIDON = 1`, so no-lid files stay identical).
- [ ] **Step 4:** Build, run both suites. Expected: all pass.
- [ ] **Step 5:** Document the `Damping_lid` line and `DampingLid.pnl` in `docs/input_files.md`:
  - the file format and global coordinates
  - phi = 0 means total potential, and the lid reflects, it does not absorb
  - the stops
  - lid panels below about 0.01 h in finite depth share the scatter of `FINGREEN3D` measured in `limes-measurement-result.md`

  Start `test/data/damping_lid/README.md` with the case table.
- [ ] **Step 6: Commit** `feat: read and check a damping lid mesh`.

### Task 2: Rankine integrals over lid panels

**Files:**
- Modify: `src/SingularIntgr.f90` (`SGLINTBD_TRI2`, `SGLINTBD_QUAD2`), `src/SingularIntgrMulti.f90` (`SGLINTBD_MULTI_TRI2`, `SGLINTBD_MULTI_QUAD2`)
- Modify: `test/test_panel_integrals.f90` (make `reference` public)
- Create: `test/test_damping_lid.f90`; modify `test/unit_tester.f90` (suite `damping_lid`, serial)

**Interfaces:**
- Consumes: `Lid_mod` arrays.
- Produces: `IRR = 5` in all four routines. Source vertices are `LXYZ(LNCON(JEL,k),:)` and `JEL` indexes the lid. `RK(1)` is the integral of 1/r and `RK(2:4)` its source gradient.

- [ ] **Step 1: Failing tests** in `test/test_damping_lid.f90`, module `test_damping_lid`, `collect_damping_lid`:
  - `lid_rankine_quad`: for a 2 x 1 quad on z = 0 and field points at its centroid, at a node, at (0.3, 0.2, 0) inside, and at (5, 0, -1), `RK(1)` equals `reference(...)` within `1d-8` relative.
  - `lid_rankine_triangle`: same for an obtuse-at-node-2 triangle.
  - `lid_node_order`: the same panel in every cyclic and reversed vertex order gives the same `RK(1)` within `1d-12` relative.
  - `lid_multi_matches_single`: the `_MULTI` routines with `IRR = 5` equal the single-body ones bit for bit.
  - Run: fail with `Please define where the field point locates.`
- [ ] **Step 2:** Add the `IRR.EQ.5` branches. If `lid_node_order` fails for reversed order, take the absolute value of the integral for `IRR = 5` only, and record why in a comment.
- [ ] **Step 3:** Build, run unit tests. Expected: pass.
- [ ] **Step 4: Commit** `feat: Rankine panel integrals for damping lid sources`.

### Task 3: Single-layer integral of G over a lid panel

**Files:**
- Modify: `src/DampingLid.f90`; test: `test/test_damping_lid.f90`

**Interfaces:**
- Consumes: `IRR = 5` integrals, `INFGREEN3D`, `FINGREEN3D`, `HavelockGF`.
- Produces: `SUBROUTINE LID_SINGLE_LAYER(XF, LEL, SL)`, with `REAL*8 XF(3)` (any field point, z <= 0), `INTEGER LEL`, and `COMPLEX*16 SL`, the integral of G over lid panel `LEL`. It reads `V`, `WVN`, `NK` and `H` from `WaveDyn_mod` as the existing routines do.

The parts of SL:

- **Rankine:** `RK(1)` with `IRR = 5`, plus the image. At z_F = 0 the image is `RK(1)` again; otherwise it is `SGLINTBD_*2` at the mirrored point with `IRR = 5`. Use the same signs as `IMAGE_TERM_EXACT`, including `SGN = -1` when `ABS(V+1) < 1e-8`.
- **Wave term, deep water:** `-V*GF` from `HavelockGF`, no Rankine parts (call `INFGREEN3D` with `FLAG = 1` and subtract `1/R1`). Choose the rule by distance:
  - dist(XF, centroid) > 3 `LPNSZ(LEL)`: one point at the centroid times `LDS(LEL)`.
  - Otherwise, field point inside the panel's projection: split into triangles at the projection of XF, then a Duffy transform with 6 x 6 Gauss-Legendre per triangle.
  - Otherwise: 6 x 6 Gauss on the panel.
- **Wave term, finite depth (H > 0):** the same rule applied to the deep-water wave term at the same `V`, plus `LDS(LEL)` times the remainder `R(rho) = G_finite - G_deep`. The remainder is evaluated at rho = max(dist, 0.03 H) on z = 0, both from `FINGREEN3D` and `INFGREEN3D` with `FLAG = 1`, wave parts only.
- **Limits:** at `ABS(V) < 1e-8` and at infinite frequency the wave term is zero and only the Rankine and image parts remain.

- [ ] **Step 1: Measure first.** In a throwaway program in the scratchpad, print `R(rho)` at z = 0 for rho/H in {1e-3, 3e-3, 0.01, 0.03, 0.05, 0.1, 0.3}, with H = 50 and nu in {0.05, 0.19, 1.0}. Compare against the direct series reference `h1/ref_green.f90` of the `measure-limes-scatter-in` session's scratchpad if it is reachable; otherwise compare against `LINTON` at rho/H >= 0.03.
  - Expected: R varies smoothly from rho/H = 0.03 up, and the jump in R between rho/H = 0.03 and 0.01 is below 1% of abs(G) there.
  - Record the numbers in `test/data/damping_lid/README.md`.
  - If R is not bounded as rho falls, the log coefficients differ. Stop and report to the user before going on.
- [ ] **Step 2: Failing tests:**
  - `duffy_converges`: on a unit square and a triangle, with the field point at the centroid, a node, and inside off the centroid, `nu = 0.19`, deep water: the 6 x 6 result equals a 24 x 24 Duffy result within `1d-6` relative.
  - `far_rule_consistent`: at dist = 3 `LPNSZ`, the one-point and Gauss rules agree within `1d-3` relative.
  - `electrified_disc`: `ABS(V) < 1e-8`, disc of radius 1 meshed with `annulus`-style polar quads plus a centre ring of triangles at 3 refinements (about 64, 256, 1024 panels). Solve `sum_l' SL(p_l, l') q_l' = 1`. Total charge `sum q LDS` approaches `1/pi` (kernel 2/r on z = 0, capacitance 2a/pi). The error falls at each refinement and is below 1% at the finest.
  - `finite_tends_to_deep`: `H = 5000`, `nu = 1`: finite equals deep within `1d-4` relative for a self panel of size 1.
- [ ] **Step 3:** Implement `LID_SINGLE_LAYER`.
- [ ] **Step 4:** Run unit tests. Expected: pass. Record the disc errors in the README.
- [ ] **Step 5: Commit** `feat: single-layer integral of the wave Green function over lid panels`.

### Task 4: Single-body solve

**Files:**
- Modify: `src/AssbMatx.f90` (move the fill loop of `ASSB_LEFT` into `SUBROUTINE FILL_HULL_LEFT(AMAT,NELEM,NSYS)`, which `ASSB_LEFT` calls before factorising)
- Modify: `src/DampingLid.f90`, `src/HAMS_Prog.f90` (single-body frequency loop)
- Test: `test/test_integration_lid.f90`, fixtures script, README

**Interfaces:**
- Consumes: `FILL_HULL_LEFT`; `ASSB_RBC` and `ASSB_DBC` for the hull rows of the right-hand side; `LID_SINGLE_LAYER`; `CALGREEN` (no IRSP) or `CALGREEN_IRR` (IRSP).
- Produces:
  - `SUBROUTINE LID_RADIATION(NELEM,NSYS)` and `SUBROUTINE LID_DIFFRACTION(NELEM,NSYS)` fill `MXPOT(1:NELEM,MD,:)` and `MXQ(:,MD,:)` the way `RADIATION_SOLVER` and `DIFFRACTION_SOLVER` fill `MXPOT`, including the `RXY` recombination and the `/NSYS`.
  - The matrix is built once per frequency by `SUBROUTINE LID_ASSEMBLE(NELEM,NSYS)`, kept in module arrays `ALID`, `IPIVL` (and `INELEM` interior rows when `IRSP = 1`).

System for each `IP`. Unknowns are `[phi(1:NELEM), q(1:LNELEM)]`; rows are hull, lid, then `IRSP` interior.

- **Hull block:** `FILL_HULL_LEFT`.
- **Hull rows, lid columns:** `-sum_IS RXY(IS,IP) SL(XF_IS(i), l)`. `XF_IS` is the hull centroid reflected as in `CALGREEN`.
- **Lid rows, hull columns:** the double-layer term at lid field points, built as the `IRR = 3` loop of `CALGREEN_IRR` with the lid arrays: `SGLINTBD_*2(XF,JEL,RK,1)`, the Green function, and `IMAGE_TERM_EXACT` when `FLAG = 1`. Diagonal 0.
- **Lid rows, lid columns:** `-sum_IS RXY(IS,IP) SL(XF_IS(l), l')`.
- **Interior rows (IRSP):** the `IRR = 3` rows of `ASSB_LEFT_IRR` plus lid columns, as for lid rows.
- **Right-hand side, radiation:** hull rows from `ASSB_RBC`. Lid and interior rows get `sum_j S_lj V_j`, with S from the same integrals as the lid-row double layer and the mode signs of `RBC_IRR`.
- **Right-hand side, diffraction:** hull rows from `ASSB_DBC`. Lid rows get the analogue of `DBC_IRR`, plus `sum_IS RXY(IP,IS) 4 pi VINP(XF_IS(l))`. Interior rows get no incident term.
- **Solve:** `ZGETRF`/`ZGETRS` when there are no interior rows, `ZGELS` with them.
- **Infinite frequency (`ABS(V+1) < 1e-8`):** call the existing no-lid or `IRSP` path and set `MXQ = 0`.

- [ ] **Step 1: Fixtures.** Annulus r = 6 to 10:
  - `lid_ann_40` (20 x 2), `lid_ann_160` (40 x 4), `lid_ann_640` (80 x 8)
  - `lid_ann_160_rev` (reversed vertex order) and `lid_ann_tri` (each quad of 40 x 4 split in two)
  - `lid_touch` (r = 5 to 9, 40 x 4) and `lid_patches` (sectors 0 to 60 and 180 to 240 degrees of the 40 x 4 annulus)
  - `lid_shrink_1` (r = 6 to 6.5, 40 x 1) and `lid_shrink_2` (r = 6 to 6.25, 40 x 1)
  - `lid_deep` (`lid_ann_160` with depth -1)
  - `lid_limits` (`lid_ann_40` with 0 and infinite frequency limits on; look up the control value in `docs/input_files.md` before writing it)
  - Field points of every case: lid centroids of its first 10 panels.
- [ ] **Step 2: Failing tests:**
  - `lid_runs`: all cases exit 0, and no NaN in `AmssDamp.1` or `ExcForce.3`.
  - `lid_damping_nonnegative`: B11, B22, B33 >= -1e-6 max(abs(B)) for every case.
  - `lid_converges`: for A11, A33, B11, B33 and abs(X1), the change from 160 to 640 is smaller than the change from 40 to 160. Report the observed rate.
  - `lid_symmetric_ab`: abs(A_ij - A_ji) and abs(B_ij - B_ji) below 1e-2 of max(abs(A)) on `lid_ann_640`, `lid_touch` and `lid_patches`. Set the final limit from the measurement.
  - `lid_haskind_heave`: as `haskind_heave` on `lid_ann_640`, limit 2e-2 at first.
  - `lid_order_free`: `lid_ann_160_rev` equals `lid_ann_160` within 1e-10 relative.
  - `lid_tri_close`: `lid_ann_tri` within 2% of `lid_ann_160` for A33, B33, abs(X1).
  - `lid_shrink_to_none`: relative distance from `image_cos` falls from `lid_shrink_1` to `lid_shrink_2`.
  - `lid_limits_finite`: no NaN, and the infinite-frequency A of `lid_limits` equals that of a no-lid run of the same hull within 1e-10.
- [ ] **Step 3:** Implement `FILL_HULL_LEFT` and run the full suites before going on. Expected: all pre-existing tests pass with unchanged stored outputs.
- [ ] **Step 4:** Implement `LID_ASSEMBLE`, `LID_RADIATION` and `LID_DIFFRACTION`, and branch the single-body loop: `IF (LIDON.EQ.1 .AND. ABS(V+1.D0).GE.1.D-8)`.
- [ ] **Step 5:** Run both suites. Expected: pass. Write the measured values and the limits set from them to the README, then tighten each limit to about twice its measured value.
- [ ] **Step 6: Commit** `feat: solve the single-body problem with a damping lid`.

### Task 5: Field output with a lid

**Files:**
- Modify: `src/PressureElevation.f90` (`CalPotential`); tests and fixtures

**Interfaces:**
- Consumes: `MXQ`, `LID_SINGLE_LAYER`.
- Produces: inside `CalPotential`, before the division by `SLD`, `POT = POT + sum_l sum_IS SL(XP_IS, l) MXQ(l, MD, IS)`, with the same reflection of `XET` as the hull loop. Only when `LIDON = 1`.

- [ ] **Step 1: Fixtures.** `lid_field_off_centroid`: `lid_ann_160` with field points at a lid node, an edge midpoint, a point inside a panel off its centroid, the 10 first centroids, and the points (0, 8, 0) and (20, 0, 0).
- [ ] **Step 2: Failing tests:**
  - `lid_elevation_zero_at_centroids`: diffraction elevation magnitude at lid centroids below 1e-6 of the incident amplitude, on `lid_ann_640`. Radiation elevation below 1e-6 of the largest field-point value.
  - `lid_field_off_centroid`: finite at node, edge and inside points. Report the values; a first limit is 0.1 of the incident amplitude at 640 panels, tightened from the measurement.
  - `lid_field_matches_collocation`: as `field_points_match_collocation` in `test_integration_image.f90`, on `lid_ann_160` with field points at hull centroids.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run, record, tighten, and commit `feat: field output includes the damping lid`.

### Task 6: Symmetry and IRSP with a lid

**Files:** fixtures, tests, and fixes in `src/DampingLid.f90` if the tests show them.

- [ ] **Step 1: Fixtures:**
  - `lid_ann_160_half`: `image_cos_half` hull, `ISY = 1`, half annulus.
  - `lid_ann_160_irsp`: `IRSP = 1` with the uniform waterplane mesh of `image_cos_lid`, plus the annulus.
- [ ] **Step 2: Tests:**
  - `lid_half_matches_full`: A, B, X of half equal full within 1e-4 relative, the limit of `half_matches_full`.
  - `lid_irsp_matches`: `lid_ann_160_irsp` equals `lid_ann_160` within 5e-3 relative, the limit of `lid_matches_no_lid`.
- [ ] **Step 3:** Run. If either fails, find the cause with superpowers:systematic-debugging before changing code.
- [ ] **Step 4:** Commit `test: damping lid with symmetry and with IRSP`.

### Task 7: Independent reference

**Files:**
- Create: `pre_postprocessing_python/DampingLid/annular_lid_reference.py`
- Fixtures: `mm_cyl_nolid`, `mm_cyl_lid_a`, `mm_cyl_lid_b`
- Tests in `test/test_integration_lid.f90`; README

Bottom-mounted vertical cylinder, radius a = 2, depth h = 10, omega such that kh = 1.5, heading 0. The hull runs from z = 0 to z = -h, 24 around and 10 down with cosine rows, and has no bottom. The lid is the annulus a to b = 4, at two refinements (24 x 2 and 48 x 4). Keep panels above 0.03 h if the mesh allows; if not, say so in the README.

The script solves the problem by eigenfunction matching, one azimuthal order m at a time:

- **Lid region (a < r < b):** vertical modes `cos((n - 1/2) pi (z + h)/h)` with radial functions `I_m` and `K_m`; `dphi/dr = 0` at r = a.
- **Outer region (r > b):** the propagating mode `cosh k(z+h)` with `J_m` + `H_m^(1)`, and evanescent modes `cos k_n (z + h)` with `K_m`.
- **Matching:** phi and dphi/dr at r = b, projected on the outer modes, N terms.
- **Output:** surge force and elevation at (0, 8, 0) and (8, 0, 0), with time factor and normalisation matching HAMS output, read from `docs/output_files.md`. It prints N = 20, 40, 80 to show convergence, which is slow because of the edge singularity. Extrapolate in 1/N and report the spread.
- **Self-check:** with b = a the script must reproduce MacCamy-Fuchs to 1e-8.

- [ ] **Step 1:** Write the script and its self-check, and run it. Expected: the self-check passes and the extrapolated values are stable to 3 digits.
- [ ] **Step 2: Failing tests:**
  - `mm_nolid_matches_maccamy_fuchs`: surge force of `mm_cyl_nolid` within 2% of MacCamy-Fuchs. This checks the hull mesh on the sea bed before any lid comparison.
  - `mm_lid_matches_reference`: surge force and both elevations of `mm_cyl_lid_b` within 3% of the reference, and closer than `mm_cyl_lid_a`.
  - The reference values are stored as literals in the test, with the script and N in a comment.
- [ ] **Step 3:** Run. If `mm_lid_matches_reference` fails while the no-lid case passes, check the signs in Task 4 first. Record the numbers in the README.
- [ ] **Step 4:** Commit `test: compare the damping lid with eigenfunction matching`.

### Task 8: Multiple bodies

**Files:**
- Create: `src/DampingLidMulti.f90`
- Modify: `src/HAMS_Prog.f90` (multi-body loop), `src/AssbMatxMulti.f90` (`FILL_HULL_LEFT_MULTI`, as in Task 4), `src/PressureElevationMulti.f90`
- Tests and fixtures

**Interfaces:**
- Consumes: the combined arrays `XYZ_GLOBAL_MULTI_COMB_P`, `NCN_MULTI_COMB`, `NCON_MULTI_COMB`, `GEONRM_MULTI_COMB`, `DS_MULTI_COMB`, `PNSZ_MULTI_COMB`, `CGRN_MULTI_COMB`, `RKBN_MULTI_COMB`; `SGLINTBD_MULTI_*2` with `IRR = 5`; `LID_SINGLE_LAYER` (geometry of the lid is global already).
- Produces: `LID_ASSEMBLE_MULTI`, `LID_RADIATION_MULTI` and `LID_DIFFRACTION_MULTI` fill `MXPOT_MULTI_COMB` and `MXQ` with columns `6*NBODY+1`, as `RADIATION_SOLVER_MULTI` and `DIFFRACTION_SOLVER_MULTI` do. They use the per-body radiation right-hand sides as the `ASSB_RBC_MULTI` calls in `HAMS_Prog.f90:630-643`.

- [ ] **Step 1:** Read the multi-body `IRSP` path in `HAMS_Prog.f90:560-720` and write in the README whether it is used here. The spec says reuse it only after checking it. Use the plain multi-body path for hull rows; add interior rows only if Task 6's approach carries over without change.
- [ ] **Step 2: Fixtures:**
  - `lid_two_far`: `image_cos_two` hulls with the `lid_ann_160` annulus around body 1.
  - `lid_gap`: two boxes 20 x 8 x 4 m, 4 m apart, side by side, with the gap lid 20 x 4 m at 20 x 4 panels.
- [ ] **Step 3: Failing tests:**
  - `lid_multi_matches_single`: body 1 of `lid_two_far` within 1e-3 of `lid_ann_160`, the limit of `multi_matches_single`.
  - `lid_gap_symmetric_ab`: 12 x 12 A and B symmetric within 1e-2 of max.
  - `lid_gap_damping_nonnegative`: eigenvalues of B >= -1e-6 max(abs(B)).
  - `lid_gap_elevation_zero`: gap elevation at lid centroids below 1e-6.
- [ ] **Step 4:** Implement, run, record, tighten, and commit `feat: damping lid for multiple bodies`.

### Task 9: Developer documentation and final check

**Files:** `docs/developer_documentation.md`, `readme.md`, `test/data/damping_lid/README.md`

- [ ] **Step 1:** Add a "Damping lid" section to `docs/developer_documentation.md`:
  - the formulation and why the `IRSP` rows alone fail, with the measured table from the spec
  - the self-term scheme and its measured accuracy
  - the convergence rate and the eigenfunction-matching comparison
  - the limits: no GMRES, the finite-depth small-panel caveat, a reflecting rather than absorbing lid

  Keep measured and derived separate. Add one line to `readme.md` on the new input.
- [ ] **Step 2:** Run the full suites one last time from a clean build (`rm -rf build install` first). Expected: all pass. Record pass counts against Task 0.
- [ ] **Step 3:** Commit `docs: damping lid formulation and validation`. Report to the user; no push.
