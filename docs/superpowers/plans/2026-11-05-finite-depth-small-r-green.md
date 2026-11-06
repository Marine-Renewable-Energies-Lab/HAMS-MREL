# Finite-Depth Green Function at Small R: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stop the `Limes` scatter in `PIDCOCK` from reaching results, by using `LINTON` for R/h < 0.02 when nu*h <= 50, where nu = omega^2/g.

**Architecture:** One changed branch test in `FINGREEN3D_R`. `PIDCOCK` keeps R/h from 0.02 to 0.05, and all of R/h < 0.05 above the nu*h limit. `EIGENE`, `EIGEN`, `Limes` and `LINTON` do not change.

**Tech Stack:** Fortran, fpm with test-drive. ifort 2022 locally, CI on Linux and Windows. fpm is at `C:\Users\vraghavan1\Downloads\MREL_HAMS-MREL\tools\fpm\fpm.exe`.

**Spec:** The measurement in `C:\Users\vraghavan1\AppData\Roaming\Claude\side-session-notes\local_f93f9438-3cb3-4638-9dde-5b195bdd546b\limes-measurement-result.md`, as revised by "Why these limits" below.

**Branch:** `fix/finite-depth-small-r-green`, cut from `main` at 6f4cc3b.

## Why these limits

All numbers are from scratch runs of the rule below, with the direct sum as reference. They replace the option A in the measurement, which switched at R/h = 0.05 with no limit.

**R/h = 0.02.** Scaled error in dG/dR (times h^2), over 6 depth pairs, `LINTON` / HAMS today:

| R/h | nu*h = 1, h = 10 | nu*h = 20, h = 100 | nu*h = 50, h = 100 |
| :--- | ---: | ---: | ---: |
| 0.0006 | 6e-6 / 8e2 | 4e-4 / 4e3 | 2e-4 / 8e3 |
| 0.005 | 2e-7 / 2 | 0.1 / 50 | 0.1 / 100 |
| 0.01 | 2e-7 / 6e-2 | 0.1 / 3 | 0.2 / 8 |
| 0.02 | 1e-7 / 2e-2 | 0.2 / 0.8 | 0.3 / 0.1 |
| 0.04 | 2e-7 / 1e-4 | 0.3 / 3e-3 | 0.3 / 2e-2 |

`PIDCOCK` is the better routine from R/h of about 0.02 upward, so a switch at 0.05 gave up accuracy for nothing. The two are about equal at 0.02.

**nu*h = 50.** Line 755 of `src/FinGreen3D.f90` sets `A*k*h = 0.7*sqrt(nu*h)`, clamped at 5, and the clamp is reached at nu*h = 51. The largest `LINTON` error in dG/dz below R/h = 0.05 rises steeply above that, at h = 10: 0.2 at nu*h = 50, 0.7 at 60, 5 at 80, 28 at 100, 530 at 200. The clamp is the likely cause, not proven. The measurement named the cubic rule as first suspect, and that was wrong: line 755 overrides it.

**On real cases**, worst of A11, A33, A55, B11, B33, B55, |X1|, |X3|, |X5|, percent against the direct sum:

| Case | nu*h | HAMS today | This plan | R/h < 0.02, no limit |
| :--- | ---: | ---: | ---: | ---: |
| a = 1, h = 10 | 20 | 0.018 | 0.035 | 0.035 |
| a = 1, h = 30 | 30 | 0.030 | 0.045 | 0.045 |
| a = 1, h = 30 | 60 | 0.19 | 0.19 | **1.18** |
| a = 1, h = 100 | 50 | 0.53 | 0.32 | 0.32 |
| a = 1, h = 100 | 100 | 0.81 | 0.81 | 0.62 |
| a = 1, h = 100 | 200 | 4.2 | 4.2 | 4.9 |
| a = 1, h = 300 | 150 | 4.2 | 4.2 | 0.32 |
| a = 1, h = 300 | 300 | 9.2 | 9.2 | 0.89 |
| example 1 body 1, h = 0.9 | 20 | 0.25 | 0.011 | 0.011 |

Term I, 32 x 16, a = 5, depth 10, against MacCamy-Fuchs:

| Case | HAMS today | This plan | Direct sum |
| :--- | ---: | ---: | ---: |
| Offset 0.01, 1.0 rad/s | -11.40% | -1.94% | -1.96% |
| Offset 0.01, 1.6 rad/s | -3.56% | -3.46% | -3.49% |
| Radius scaled 1 + i x 1e-10, i = 1, 2, 8, 1.0 rad/s | up to -13.77% | -1.94% each | not run |
| 48 x 16, offset 0.01, 1.6 rad/s | -4.38% | -3.53% | -3.56% |

**What this does not fix.**
- Above nu*h = 50, results stay as they are today, with the jumps. On a = 1, h = 300 that means up to 9%, where the deep-water route is correct.
- `EIGENE` has a bias of up to 13 in scaled dG/dR at R/h = 0.05 to 0.08, for nu*h of 0.1 or less. Its effect on results has not been separated.
- `LINTON` errors at h = 100 are larger than at h = 10 for the same nu*h and R/h, so something in it does not scale with h. The cause is not known. The errors are still at or below the `PIDCOCK` ones below R/h = 0.02.

## Global Constraints

- Limits, exact: `RHP < 0.02D0` and `V*H <= 50.D0`. The existing `RHP < 0.0005D0` stays unconditional.
- Do not change `LINTON`, `PIDCOCK`, `EIGENE`, `EIGEN` or `Limes`. A separate session is working on NaN results at very large kh in those routines' depth factors, so keep this change to the branch test and the header comments.
- Do not touch `second-order-drift` from this branch. Its README statement and its `waterline_method = field` warning are updated there after this merges.
- Do not regenerate a stored test reference without the user's agreement.
- Long runs do not go on the laptop. Snellius by SSH key only, with the job size shown to the user first.
- Nothing is pushed and no PR is opened until the user asks.
- PR bodies and release notes carry no "Generated with Claude Code" line.

## Review Focus

1. A frequency sweep crossing nu*h = 50 steps by up to about 0.8% on panels this small, as the routine changes. Expected, and it goes in the README text (Task 2). Pinned by `linton_limit`.
2. Very low frequency: example 5 reaches nu*h = 0.018. Swept down to 0.001 with no problem. Pinned by `small_r_matches_direct_sum`.
3. Irregular-frequency lid points at z = 0, which `LINTON` has only met at R/h < 0.0005 until now. Not yet measured. Pinned by the z = 0 pairs in `small_r_matches_direct_sum`, and by example 5 in Task 2.
4. A pair exactly at R/h = 0.02. Pinned by `linton_pidcock_continuity`, moved to the new boundary.
5. Run time on meshes with many close pairs: `LINTON` costs about 4 times `PIDCOCK` per call. Measured on example 1 body 1 and example 5 in Task 2.

---

### Task 1: Tests and the branch test

**Files:**
- Create: `test/direct_sum_green.f90`
- Modify: `test/test_fingreen3d.f90`
- Modify: `src/FinGreen3D.f90:179-191`, and the region lines in the headers of `EIGENE` (about line 367), `PIDCOCK` (about line 512) and `LINTON` (about line 666)

**Interfaces:**
- Produces: `subroutine direct_sum_green(R, ZF, ZP, V, WK, H, GRN)` in module `direct_sum_green_mod`, with `real(8), intent(in)` arguments and `complex(8), intent(out) :: GRN(3)`. Same return convention as `FINGREEN3D_R` with TAG = 1.
- Consumes: `IK01A`, `JY01B` from `FinGrnExtSubs`. `FINGREEN3D_R`, `LINTON`, `PIDCOCK` from `FinGreen3D_Open`, all public because the module has no `PRIVATE`.

- [ ] **Step 1: Write `direct_sum_green`.** The tested version is `C:\Users\VRAGHA~1\AppData\Local\Temp\claude\7ad205b7dbee283c\95449479-1bac-4a72-a29e-80d701cfbf32\scratchpad\limes\h1\ref_green.f90`, which may not outlive the session, so the formulas are given here. It sums the eigenfunction series with no acceleration, with Neumaier compensation, and stops where um*R > 40. Roots come from Newton iteration on `um*tan(um*h) = -V`. Then subtract 1/SR and 1/SR1, as `EIGEN` does:

```fortran
GRN(1) = P0*(CI*J0 - Y0)  + sum(PM*K0*CF*CP)        - (1/SR + 1/SR1)
GRN(2) = -WK*P0*(CI*J1 - Y1) - sum(PM*UM*K1*CF*CP)  + (R/SR**3 + R/SR1**3)
GRN(3) = WK*DP0*(CI*J0 - Y0) - sum(PM*UM*K0*SF*CP)  + ((ZF-ZP)/SR**3 + (ZF+ZP)/SR1**3)
! with NM0 = H/2*(1 + sinh(2*WK*H)/(2*WK*H)), P0 = PI/NM0*cosh(WK*(ZF+H))*cosh(WK*(ZP+H)),
! DP0 = PI/NM0*sinh(WK*(ZF+H))*cosh(WK*(ZP+H)), PM = 4/(1 + sin(2*UM*H)/(2*UM*H))/H,
! CF = cos(UM*(ZF+H)), CP = cos(UM*(ZP+H)), SF = sin(UM*(ZF+H))
```

- [ ] **Step 2: Write the failing tests** in `test/test_fingreen3d.f90`, registered in `collect_fingreen3d`. Grid: h in [10, 100]; nu*h in [0.018, 1, 20, 50]; R/h in [0.0006, 0.002, 0.005, 0.01, 0.015]. Depth pairs (zf/h, zp/h): (-0.003, -0.03), (-0.01, -0.003), (-0.03, -0.03), (-0.1, -0.2), (-0.3, -0.5), (-0.06, -0.06). Set V = nu*h/H and W = sqrt(V*GRAV), then call `DISPERSION`.

```fortran
! small_r_no_scatter: 20 samples of R*(1 + k*2d-11), on the grid
call check(error, maxval(dGdR_samples*H**2) - minval(dGdR_samples*H**2) < 1d0, ...)  ! same for dG/dz, real and imaginary

! small_r_matches_direct_sum: on the grid, plus depth pairs (0, -0.03) and (0, 0)
call check(error, abs(g(1) - ref(1))*H    < 0.05d0, ...)
call check(error, abs(g(2) - ref(2))*H**2 < 0.5d0,  ...)
call check(error, abs(g(3) - ref(3))*H**2 < 0.5d0,  ...)

! linton_limit: H = 100, zf/h = -0.01, zp/h = -0.03; each exact, all three components
! nu*h = 49.9, R/h = 0.01    -> equals LINTON
! nu*h = 60,   R/h = 0.01    -> equals PIDCOCK
! nu*h = 60,   R/h = 0.0004  -> equals LINTON
! nu*h = 20,   R/h = 0.0199  -> equals LINTON
! nu*h = 20,   R/h = 0.0201  -> equals PIDCOCK
```

Move `linton_pidcock_continuity` from `5d-4*H` to `2d-2*H`, and update the header comment of the test file.

- [ ] **Step 3: Run the tests to verify they fail.**

Run: `..\..\..\..\tools\fpm\fpm.exe test --target unit-tester` from the worktree, with the Intel oneAPI environment set (`setvars.bat intel64 vs2019`, with `NoDefaultCurrentDirectoryInExePath` cleared first).
Expected: `small_r_no_scatter`, `small_r_matches_direct_sum` and `linton_limit` FAIL. The other tests pass.

- [ ] **Step 4: Change the branch test in `FINGREEN3D_R`.** Declare `RLIN` and `VHLIN` as named constants, with a comment giving the reason for each and pointing at the clamp in `LINTON`.

```fortran
IF (RHP.LT.0.0005D0.OR.(RHP.LT.RLIN.AND.V*H.LE.VHLIN)) THEN
  CALL LINTON(RR,ZF,ZP,V,WK,WVN,NK,H,FG)
ELSEIF (RHP.LT.0.05D0) THEN
  CALL PIDCOCK(RR,ZF,ZP,V,WK,WVN,NK,H,FG)
```

Correct the region lines in the three headers to what the code does.

- [ ] **Step 5: Run the tests to verify they pass.**

Run: the same command.
Expected: all unit tests PASS. If only the z = 0 pairs fail, stop and report the size. That case has not been measured, and the tolerance must not be widened to fit it.

- [ ] **Step 6: Commit.**

```bash
git add test/direct_sum_green.f90 test/test_fingreen3d.f90 src/FinGreen3D.f90
git commit -m "Use LINTON below R/h 0.02 up to nu*h 50 in the finite-depth Green function"
```

### Task 2: Regression, real cases and the note

**Files:**
- Modify: `readme.md`, next to the other v2.3.3 corrections
- Test data: `test/data/config1`, `test/data/config1_gmres`, `test/data/image_term`, only if the user agrees

**Interfaces:**
- Consumes: the `FINGREEN3D_R` rule from Task 1.

- [ ] **Step 1: Run the integration tests.**

Run: `fpm.exe install --prefix install`, then `fpm.exe test --target integration-tester -- install/bin/hamsmrel`.
Expected: PASS at the stored tolerances (A 5e-5, B 2e-3). The census puts config1 entries at about 2e-5, so a pass is likely but not certain. If one fails, stop and report the relative change for the user to decide.

- [ ] **Step 2: Run A, B, X on the branch build.** Use the scratch run folders for a = 1 at h = 10, 30, 100 and example 1 body 1, then examples 2, 4 and 5 from the repo. Compare with HAMS today and record run times.
Expected: the "This plan" column above within 0.001%. Examples 2, 4, 5 move by at most the census size, about 0.06%. Example 5 gives no NaN. Report the run-time ratio.

- [ ] **Step 3: Check term I against the real source.** In a throwaway `git worktree` at `origin/second-order-drift` under the scratchpad, apply this `src/FinGreen3D.f90` without committing, build, and run the 32 x 16 and 48 x 16 cases at offset 0.01, plus the three scaled radii.
Expected: -1.94%, -3.46%, -3.53%, and -1.94% on all three radii. Then remove the worktree.

- [ ] **Step 4: Ask before full example 1.** It has 5 bodies and 3635 panels, and has never been run in this comparison. It goes to Snellius, with the job size shown to the user first. If declined, say in the PR that body-to-body effects at high frequency are untested.

- [ ] **Step 5: Write the README note**, in the style of the other v2.3.3 corrections: who is affected (finite depth, close panels), the numbers measured in Step 2, that results above nu*h = 50 do not change, and that those runs can use the deep-water solver. Give the step of up to about 0.8% at nu*h = 50. Say which results were checked. Follow the manuscript standard: cite or justify each claim, in simple language.

- [ ] **Step 6: Commit.**

```bash
git add readme.md
git commit -m "Note the finite-depth Green function change in the readme"
```

## After this branch

- On `second-order-drift`, once this reaches `main`: merge, re-run term I, and update the README statement of at most 0.02 points and the `field` warning.
- Above nu*h = 50: option C (direct sum above R/h = 0.01) or a fix to the clamp in `LINTON`. Neither measured.
- `EIGENE` bias at R/h 0.05 to 0.08 and nu*h of 0.1 or less.
