# Dissipative Damping Lid Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the phi = 0 lid by the WAMIT-style dissipative lid p/rho = epsilon dphi/dz, with epsilon per damper range, for one and several bodies, without changing runs without a lid.

**Architecture:** The lid unknown becomes phi (was q). The lid columns of the existing assembly are scaled by -alpha_l, alpha = i omega epsilon nu / (g - i omega epsilon), the lid rows gain 4 pi on the diagonal, the diffraction right-hand side gains sum_l alpha_l S phi_I instead of 4 pi phi_I, and the field and elevation use the new lid term. Everything stays in `src/DampingLid.f90`, `src/DampingLidMulti.f90` and the lid branches of `PressureElevation*.f90` and `HAMS_Prog.f90`.

**Tech Stack:** Fortran 90 (ifort 2022.1), MKL LAPACK, fpm with test-drive, Python 3.7 + numpy 1.21 + scipy 1.7.

**Spec:** `docs/superpowers/specs/2026-11-06-epsilon-damping-lid-design.md`

## Global Constraints

- Branch `feature/damping-lid` (worktree `C:\Users\vraghavan1\Downloads\MREL_HAMS-MREL\HAMS-MREL-damping-lid`). Never push or open a PR. Never touch `C:\Users\vraghavan1\Downloads\MREL_HAMS-MREL\HAMS-MREL`. Never stash, reset or switch branches.
- First-order formulation unchanged: do not edit CALGREEN, the hull D and S, panel integrals, `SingularIntgr*`, `AssbMatx*`, Green function files, solver files, `FILL_HULL_LEFT*`. Lid code only adds terms.
- Runs without `Damping_lid 1` give bit-identical results. Never regenerate stored expected outputs of other suites.
- Build and test only with ifort through `cmd //c "C:\Users\vraghavan1\Downloads\MREL_HAMS-MREL\validation-runs\tools\fpmenv_baseline.bat" C:\Users\vraghavan1\Downloads\MREL_HAMS-MREL\HAMS-MREL-damping-lid <fpm args>`.
- Build: `<env> build`, then `<env> install --prefix install`. Unit: `<env> test --target unit-tester`. Integration: `<env> test --target integration-tester -- install/bin/hamsmrel`. Fixtures: `python -I pre_postprocessing_python/DampingLid/make_lid_fixtures.py`.
- No full example cases on the laptop. Snellius only by key login (`ssh -o BatchMode=yes -o PasswordAuthentication=no vraghavan@snellius.surf.nl`), never a password file, and only after the controller has shown job size and cost to the user.
- Do not read Nemoh source.
- Commit trailer: `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Commit only files the task changed. Do not commit `Damping_lid_report.docx` or `.pdf`.
- New lid arrays are `REAL*8` / `COMPLEX*16`.
- Every limit in a test is set from a measurement recorded in `test/data/damping_lid/README.md`, with the measured value next to the limit in the test source. Docs separate measured from derived.
- Error messages start with ` Error: ` and name the damper line or panel number.

## Review Focus

1. A lid split into two damper ranges with different epsilon must use each panel's own epsilon in the matrix, the field term and the elevation. Test: `lid_two_dampers_split` in Task 2 (two ranges with equal epsilon equal one range, to 1e-12).
2. Symmetric meshes (ISY = 1) must give the images the epsilon of their source panel. Test: `lid_half_matches_full` in Task 2 with epsilon = 0.2.
3. A field point on z = 0 on a lid node or edge must give a finite elevation. Test: `lid_field_off_centroid` in Task 2.
4. epsilon = 0 must reproduce the no-lid hull results exactly, for one and two bodies. Tests: `lid_eps0_matches_no_lid` (Task 2), `lid_gap_eps0_matches_no_lid` (Task 3).
5. A zero frequency in a lid run must stop, not give NaN. Test: `stops_lid_zero_frequency` in Task 1.

---

### Task 1: Damper input, stops, rename

**Files:**
- Modify: `src/WavDynMods.f90` (`Lid_mod`), `src/DampingLid.f90` (`READ_LID_MESH`, `ALLOCATE_LID_UNKNOWNS`), `src/DampingLidMulti.f90`, `src/PressureElevation.f90`, `src/PressureElevationMulti.f90`, `src/HAMS_Prog.f90`
- Modify: `pre_postprocessing_python/DampingLid/make_lid_fixtures.py`, `test/test_integration_lid.f90`, `test/data/damping_lid/*`

**Interfaces:**
- Produces: `Lid_mod` gains `REAL*8,ALLOCATABLE:: LEPS(:)` (epsilon of each lid panel, m/s) and `INTEGER:: NDAMP`, `INTEGER,ALLOCATABLE:: IDAMP(:,:)` (NDAMP x 2, IP1 IP2), `REAL*8,ALLOCATABLE:: EDAMP(:)`. `MXQ` renamed `MXLPOT` everywhere (same shape). The solver keeps the phi = 0 equations in this task.
- Fixture script: `lid_file(nodes, panels, isx=0, isy=0, dampers=None)`; `dampers` is a list of `(ip1, ip2, eps)`, default `[(1, len(panels), 0.2)]`; `None` passed as the string `"none"` omits the block (for the missing-block stop). The block text:

```
    # Start Definition of Damping   ! number of dampers, then: first panel  last panel  epsilon (m/s)
         2
         1        40    0.2000000000
        41        80    0.1000000000
    # End Definition of Damping
```

- [ ] **Step 1: Fixtures.** Add `dampers` to `lid_file` and to `annulus`/`rectangle` callers through `write_case`; add stop cases `lid_stop_damp_missing` (no block), `lid_stop_damp_count` (NDAMPERS 0), `lid_stop_damp_range` (IP2 = LNELEM + 1), `lid_stop_damp_order` (IP1 > IP2), `lid_stop_damp_overlap` (1-25 and 20-40 of a 40-panel lid), `lid_stop_damp_gap` (1-20 and 25-40), `lid_stop_damp_negative` (eps -0.1), `lid_stop_damp_value` (eps `abc`), `lid_stop_zero_frequency` (lid on, frequency list with 0). Change `lid_limits`/`lid_limits_off` so they contain the infinite-frequency limit but not zero frequency (read how `0_inf_frequency_limits` selects them in `InputFiles.f90`; if both come together, give the infinite limit through the frequency list instead and document it). Regenerate with the script.
- [ ] **Step 2: Failing tests.** In `test_integration_lid.f90` add `stops_lid_damp_missing`, `_count`, `_range`, `_order`, `_overlap`, `_gap`, `_negative`, `_value`, `stops_lid_zero_frequency`, each with `run_case_expect_stop(name, message, error)` and these messages: `no damping block`, `number of dampers`, `damper panel range`, `damper panel range`, `lid panel in two dampers`, `lid panel in no damper`, `negative damping`, `cannot read the damping`, `zero frequency with a damping lid`. Run: they fail.
- [ ] **Step 3: Implement** the block read at the end of `READ_LID_MESH` (skip lines until one containing `Start Definition of Damping`; EOF before it gives the missing stop), fill `LEPS` from the ranges with the checks of the spec, print one line per damper (` Damper k: panels IP1 to IP2 (n panels), epsilon = E m/s`). Add the zero-frequency stop in `HAMS_Prog.f90` in both the single- and multi-body frequency loops, before the solve, when `LIDON.EQ.1` and `ABS(V).LT.1.D-8`. Rename `MXQ` to `MXLPOT` and update the `Lid_mod` comment to "lid potential phi" (it still holds q until Task 2; say so in the comment).
- [ ] **Step 4: Run** unit and integration suites: all pass, including every pre-existing lid test (they now read fixtures with the block).
- [ ] **Step 5: Commit** `feat: damping block in DampingLid.pnl`.

### Task 2: The epsilon lid for one body

**Files:**
- Modify: `src/DampingLid.f90` (`LID_ASSEMBLE`, `LID_RADIATION`, `LID_DIFFRACTION`, `LID_UNPACK`, comments), `src/PressureElevation.f90` (`CalPotential` lid term, `CalElevation`), `src/WavDynMods.f90` comment
- Modify: `pre_postprocessing_python/DampingLid/make_lid_fixtures.py`, `test/test_integration_lid.f90`, `test/data/damping_lid/README.md`

**Interfaces:**
- Consumes: `LEPS`, `MXLPOT` from Task 1.
- Produces in `DampingLid` (public):
  - `COMPLEX*16 FUNCTION LID_ALPHA(EPS)`: `CI*W1*EPS*V/(G-CI*W1*EPS)`, 0 when `V` is 0 or -1. `W1` is in `WaveDyn_mod`; `G` and `CI` as `PressureElevation.f90` uses them (`Const_mod`).
  - `REAL*8 FUNCTION LID_EPS_AT(X, ON)`: `ON` (logical, out) true when X(3) is within 1e-6 of 0 and (x, y) of X, or of one of its symmetry images (`REFLECTED`), lies in or on the edge of a lid panel (tolerances as the on-body test in `CalPotential`); returns that panel's `LEPS`, else 0 and `ON` false.
  - Private storage `LCOL(NELEM+NROW, LNELEM, NSYS)`: the lid column block sum_IS RXY(IS,IP) S before scaling.
- MXLPOT now holds the lid phi: scattered for the diffraction column, as MXPOT.

- [ ] **Step 1: Fixtures.** All solve fixtures keep epsilon 0.2 (Task 1 default). Add `lid_ann_160_eps0` (epsilon 0), `lid_ann_160_rigid` (epsilon 1e6), `lid_ann_160_two_dampers` (ranges 1-80 and 81-160, both 0.2), `lid_ann_160_split` (1-80 at 0.2, 81-160 at 1.0), and `lid_nolid_160` (the same control file with `Damping_lid 0`, if `lid_flag_off` does not already have the same hull, frequency and field points; check and reuse when it does).
- [ ] **Step 2: Failing tests.** Replace the phi = 0 checks in `test_integration_lid.f90`:
  - `lid_eps0_matches_no_lid`: A, B, X of `lid_ann_160_eps0` equal the no-lid run to 1e-12 relative to the largest entry; the field values at points off the lid too.
  - `lid_rigid_elevation_small` replaces `lid_elevation_zero_at_centroids`: elevation at the lid centroids of `lid_ann_160_rigid` over the incident amplitude below a limit set from the measurement.
  - `lid_two_dampers_split`: `lid_ann_160_two_dampers` equals `lid_ann_160` to 1e-12; `lid_ann_160_split` differs from both by more than 1e-3 in B11.
  - `lid_absorbs`: keep `lid_damping_nonnegative` (all B_ii >= 0) and assert that the lid power at the centroids, 1/2 rho omega^2 sum epsilon_l |eta_l|^2 A_l from the diffraction elevations, is >= 0 and > 0 for epsilon 0.2.
  - Remove `lid_haskind_heave` (Haskind does not hold with epsilon > 0, spec); add `lid_haskind_eps0`: the Haskind residual of heave for `lid_ann_160_eps0` below the measured no-lid value times 2.
  - Keep and re-measure: `lid_converges`, `lid_symmetric_ab`, `lid_order_free`, `lid_tri_close`, `lid_shrink_to_none`, `lid_limits_finite` (infinite frequency equals no lid, 1e-10), `lid_field_off_centroid` (finite; relative difference from the nearest centroid value below a measured limit), `lid_field_matches_collocation`, `lid_half_matches_full`, `lid_irsp_matches`.
  Run: the new and changed tests fail.
- [ ] **Step 3: Implement** per spec "Formulation" and "Symmetry": in `LID_ASSEMBLE` fill `LCOL` where the lid columns are filled now (with +S, not -S), then `ALID(:,NELEM+L,IP) = -LID_ALPHA(LEPS(L))*LCOL(:,L,IP)` and `ALID(NELEM+L,NELEM+L,IP) += 4 pi` for the lid rows only (not interior rows), before ZGETRF. In `LID_DIFFRACTION` remove the 4 pi VINP term; for every row R of 1..NELEM+NROW and class IP add `sum_L LID_ALPHA(LEPS(L))*LCOL(R,L,IP)*UI(L,IP)` with `UI(L,IP) = sum_J RXY(IP,J)*VINP(REFLECTED(LXYZ_P(L,:),J), ...)`, skipped at V = 0. `LID_RADIATION` needs no new term (radiation phi has no incident part). In `CalPotential` replace `SLID*MXQ` by `LID_ALPHA(LEPS(LEL))*SLID*(MXLPOT(LEL,MD,IS) + PHII)` with `PHII = VINP` at `LXYZ_P(LEL,:)` reflected for image IS when MD = 7 and V is neither 0 nor -1, else 0. In `CalElevation` use `CI*W1*XPOT/(G-CI*W1*EPS)` with `EPS = LID_EPS_AT(XP, ON)` when `LIDON.EQ.1` and ON. Rewrite the routine header comments for the new unknown.
- [ ] **Step 4: Measure, set limits, run.** Run the suites; record every measured value in the README section "Single-body solve with a lid" (replace the phi = 0 numbers, say they are replaced) and set each limit from it, value next to the limit. All suites pass; no-lid suites bit-identical.
- [ ] **Step 5: Commit** `feat: dissipative damping lid for one body`.

### Task 3: The epsilon lid for several bodies

**Files:**
- Modify: `src/DampingLidMulti.f90`, `src/PressureElevationMulti.f90`, fixture script, `test/test_integration_lid.f90`, README section "Several bodies"

**Interfaces:**
- Consumes: `LID_ALPHA`, `LID_EPS_AT` from Task 2 (add them to the `USE DampingLid, ONLY:` lists).

- [ ] **Step 1: Fixtures.** `lid_gap` keeps epsilon 0.2; add `lid_gap_eps0` (epsilon 0), `lid_gap_nolid` (no lid, same field points), `lid_gap_eps1` (epsilon 1.0); `lid_two_far` keeps 0.2.
- [ ] **Step 2: Failing tests.** `lid_gap_eps0_matches_no_lid` (12 x 12 A, B and X to 1e-12 relative); `lid_gap_elevation_zero` replaced by `lid_gap_damps`: the largest gap-centroid diffraction elevation falls from `lid_gap_nolid` to `lid_gap` to `lid_gap_eps1`; keep and re-measure `lid_multi_matches_single`, `lid_gap_symmetric_ab`, `lid_gap_damping_nonnegative`.
- [ ] **Step 3: Implement** the Task 2 changes in `LID_ASSEMBLE_MULTI` (store `LCOLM(NC,LNELEM,NSYS)`), `LID_DIFFRACTION_MULTI`, the field term and elevation in `PressureElevationMulti.f90` (diffraction column `6*NBODY+1` gets the incident lid term).
- [ ] **Step 4: Measure, set limits, run** as in Task 2.
- [ ] **Step 5: Commit** `feat: dissipative damping lid for several bodies`.

### Task 4: Independent reference with epsilon

**Files:**
- Modify: `pre_postprocessing_python/DampingLid/annular_lid_reference.py`, fixture script, `test/test_integration_lid.f90`, README section "Independent reference"

**Interfaces:**
- Reference gains `--epsilon E` (m/s, default 0). The lid region uses k tanh(kh) = nu g/(g - i omega E): E = 0 is open water, a large E a rigid lid. Lid-region roots: Newton iteration on f(k) = k tanh(kh) - K with K = nu g/(g - i omega E), started from the open-water roots k_0 and i k_n, continued in small steps of E from 0 so that each root follows its open-water branch; report the residual of every root.

- [ ] **Step 1:** Implement the complex roots and use them in the existing matching. Check in `main()`: E = 0 gives MacCamy-Fuchs to the current accuracy; E = 1e6 gives the rigid-lid case; the energy balance 1/2 B_11 |U|^2 = far-field power + lid power holds to 1e-6 relative for E = 0.2 and 1 (lid power from the lid-region solution, 1/2 rho E integral |phi_z|^2 dA). Print reference X1 and the two elevations, extrapolated in 1/N as now, for E = 0.2, 1, 1e6.
- [ ] **Step 2: Fixtures.** `mm_cyl_lid_{a,b,c}` at epsilon 0.2, and `mm_cyl_lid_c_eps1`, `mm_cyl_lid_c_rigid` (1e6); lid-centroid field points added to `mm_cyl_lid_c` for the lid power.
- [ ] **Step 3: Tests.** `mm_lid_matches_reference` for the three epsilon values at `_c`; `mm_lid_converges` (error at a, b, c falling); `mm_lid_power_matches` (HAMS-MREL lid power from the centroid elevations against the reference lid power, limit from the measurement). Reference values as `complex(8), parameter` arrays with their source comment, as now.
- [ ] **Step 4: Run, record in README, set limits.**
- [ ] **Step 5: Commit** `test: eigenfunction reference for the dissipative lid`.

### Task 5: Documentation and convergence study

**Files:**
- Modify: `docs/input_files.md` (DampingLid.pnl section, Damping_lid line), `docs/developer_documentation.md` (Damping lid section), `test/data/damping_lid/README.md` (top), `readme.md` only if it mentions the lid.

- [ ] **Step 1:** Rewrite the docs for the epsilon lid: block format, stops, starting values table (WAMIT tests 0.1 to 1.0, Gourlay 2019 0.2, must be calibrated), conversion epsilon ~ g epsilon_Chen / omega (information only), outputs (hull only; elevation on the lid uses g - i omega epsilon), limits (Haskind invalid for epsilon > 0, zero frequency stops, GMRES stops, IRSP single body only). Mark the phi = 0 material as replaced, pointing to the old spec.
- [ ] **Step 2: Convergence study.** Using `mm_cyl_lid` with 16 and 32 rings and epsilon 0.2: if one run takes more than 10 minutes on the laptop, stop and report the job size to the controller (Snellius needs user approval); otherwise run locally, outside the test suite, and record the errors against the reference in the README.
- [ ] **Step 3: Commit** `docs: dissipative damping lid`.
