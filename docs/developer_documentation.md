# Development Documentation

This page contains information for the developers and mantainers of HAMS-MREL:

- [Tests](https://github.com/vaibhavraghavan/HAMS-MREL/edit/develop/docs/developer_documentation.md#tests)
- [Continous Integration](https://github.com/vaibhavraghavan/HAMS-MREL/edit/develop/docs/developer_documentation.md#continous-integration-in-github-actions)
- [Damping lid](https://github.com/vaibhavraghavan/HAMS-MREL/edit/develop/docs/developer_documentation.md#damping-lid)
- [Profiling Results](https://github.com/vaibhavraghavan/HAMS-MREL/edit/develop/docs/developer_documentation.md#profiling-results)

## Tests

HAMS-MREL includes both unit tests and integration tests written using the [Test-Drive](https://github.com/fortran-lang/test-drive) framework. The tests are designed to ensure that both individual components (subroutines, functions) and full simulation workflows are correct.

The tests are defined in [fpm.toml](../fpm.toml) as two separate executables: unit-tester and integration-tester. The unit-tester runs small-scale tests on individual subroutines and functions. The integration-tester runs full simulations for pre-defined test configurations and compares outputs with reference results. All of the data needed for the tests in included in the tests/data folder.

### Running the tests

HAMS-MREL can be built using multiple methods, but tests must always be run via **fpm** because they rely on the fpm build/test environment configured in `fpm.toml`. fpm will automatically install and link the Test-Drive framework to the test executables.

```bash
# FPM: set the compiler to IFX
export FPM_FC=IFX

# FPM: set the compiler flags
# `-heap-arrays 0` is required to keep large multi-body cases from overflowing the per-thread
# stack — see the installation docs for the full explanation.
export FPM_FFLAGS="-qmkl -qopenmp -heap-arrays 0"
export FPM_LDFLAGS="-qmkl -qopenmp"

# Compilation step. Build is handled by FPM (HAMS-MREL and test executables will be placed in build directory)
fpm build

# Place the executable in 'install' directory (an exact locaion is needed by integration tests)
fpm install --prefix install

# Run unit tests
fpm test --target unit-tester

# Run integration tests (they expect the path to the executable as an argument)
fpm test --target integration-tester -- install/bin/hamsmrel
```

The integration tests compare against the outputs stored in `test/data/config1/expected_output` with a relative tolerance, because results change slightly between compilers and flags. The measured differences behind each limit are listed in [test_integration.f90](../test/test_integration.f90). Regenerate the stored outputs only when a change is meant to alter the results.

### Free-surface image term

The Green function holds the image term 1/r1 of the free surface. Until this change it was taken at the centroid of the source panel and multiplied by the area, also for panels next to the field point, while the 1/r term was integrated over the panel. On a panel close to z = 0 compared with its size, the centroid is too close to its image: for the top row of a wall graded to the waterline, with height 0.02 m and width 0.65 m, the one-point value is 4.7 times the integral. Now, for every pair within 50 panel sizes (`FLAG = 1`), `IMAGE_TERM_EXACT` and `IMAGE_TERM_EXACT_MULTI` move the term and its gradient out of the Green function and add their integral over the panel. They are called in `CALGREEN`, `CALGREEN_IRR`, `CALGREEN_MULTI`, `CALGREEN_MULTI_IRR`, `CalPotential` and `CalPotentialMulti`, so collocation and field points are treated the same. This changes first-order results on all meshes, a little: see [the fixtures](../test/data/image_term/README.md). The stored outputs of `config1` were regenerated for it.

One-point and integrated image terms converge to the same limit at second order. On half-spheroids in deep water, from 144 to 9216 panels, with even and thin top rows, A and B from both agree with the exact limits to within 0.13% at the finest mesh, and neither is the more accurate in every mode. The integrated term makes the results depend less on the lid: at 576 panels the median change from adding one is 0.03 percentage points, against 0.08 with the one-point term.

The rest of the Green function stays at the centroid: the logarithmic part of the wave term, also singular at the free surface; in deep water the term `2V/R1` of `GRN(4)`, the z gradient of the wave part (`InfGreen_Appr.f90`), which is stronger near z = 0 than the logarithm; and with finite depth the image in the sea bed.

A hull panel with its centroid on z = 0 has no usable image term: the distance from the centroid to its own image is zero. `CalNormals` and `CalNormalsMulti` stop with a message on such a mesh. Before this check, the same mesh ran to the end and wrote NaN.

### Triangle panels

The triangle routines (`SGLINTBD_TRI`, `SGLINTBD_TRI2` and their multi-body copies) work from `ALFA`, the projection of side 1-3 on side 1-2, in units of side 1-2. Until this change they integrated over the wrong shape when `ALFA > 1`, an obtuse angle at node 2, so the result depended on the node order. This affected the direct term 1/r as well as the image term, in every earlier version. Results now change on any mesh with such a triangle as a hull panel. Lid panels are never sources, so their node order does not matter. Quadrilaterals were not affected. See [test_panel_integrals.f90](../test/test_panel_integrals.f90).

After executing these commands you should see the following output (or similar):

<img width="557" height="280" alt="hams-tests" src="https://github.com/user-attachments/assets/1bb39a74-593d-498b-8a6d-40f20fcc36f9" />

## Continous Integration in GitHub Actions

The project has [Continuous Integration (CI)](https://www.software.ac.uk/guide/how-continuous-integration-can-help-you-regularly-test-and-release-your-software) workflows in place to verify the program compiles and runs correctly. The workflows are configured with [GitHub Actions](https://docs.github.com/en/actions). They automatically compile HAMS-MREL and run the tests.

There are two workflows: one for Windows and another for Linux. Both run independently of each other and results can be seen in the [Actions tab](https://github.com/vaibhavraghavan/HAMS-MREL-DCC/actions) of the GitHub repository. The workflows are automatically triggered whenever a new commit is pushed to the repository (main or develop branch) or a pull request is opened. Running the workflows when a pull request is opened helps us verify that the changes that will be integrated compile correctly and don't cause tests to fail.

Both workflows define a single job each: *Build and Test*. The job:
- Configures the GitHub runner (server) where the workflow is executed. This means installing the Intel compilers and Intel OneAPI Math Kernel Library, Fortran Package Manager and Test-Drive framework.
- Compiles HAMS-MREL
- Runs the tests

To simplify the installation of the dependencies mentioned above, we use Actions from the [GitHub Actions marketplace](https://github.com/marketplace?type=actions). Actions can be thought of as mini-apps, extensions or plugins that do one thing specifically. In this case, the Actions install different software on the GitHub Actions runners. The following Actions are used:
- [setup-fortran](https://github.com/marketplace/actions/setup-fortran): install the Fortran Intel compiler and Intel OneAPI Math Kernel Library (Linux only)
  Although the setup-fortran action can install the Intel Fortran compiler on both Windows and Linux, installation of the oneMKL library is only available for Linux systems. Intel dependencies are manually installed on the Windows runner.
- [setup-fpm](https://github.com/marketplace/actions/setup-fpm): install the Fortran Package Manager

**The Windows workflow configuration is defined in [build_and_test_windows.yml](../.github/build_and_test_windows.yml).**

**The configuration for Linux is in [build_and_test_linux.yml](../.github/build_and_test_linux.yml).**

> [!NOTE]
> **Windows runner pinned to `windows-2022`.** Around mid-June 2026, GitHub migrated
> `windows-latest` to `windows-2025` / `windows-2025-vs2026`. The new image places
> Git for Windows' coreutils `link.exe` on PATH ahead of MSVC's linker, which breaks
> fpm + ifx's link step (`link: unknown option -- s`). The workflow is therefore
> pinned to `windows-2022` until upstream fpm or setup-fortran handles the new image
> cleanly. When you eventually migrate, the build step will need to strip
> `C:\Program Files\Git\usr\bin` (and possibly `C:\Program Files\Git\bin`) from
> PATH before `fpm build`, OR you'll need to use a different invocation that
> ensures MSVC tools are resolved first.

### About GitHub Actions

A GitHub Actions workflow is defined in a `<workflow-name>.yml` file. There can be multiple workflows, each one defined in its own file. All workflow files must be located in a `.github/workflows` folder, at the top level of the repository. 

Workflow files are written in YAML. [See the GitHub specification for workflow files](https://docs.github.com/en/actions/writing-workflows/workflow-syntax-for-github-actions) (**important link**). 

Workflows run on virtual machines provided by GitHub. These machines are called `runners`. GitHub allows us to choose among different pre-configured runners, each with a different operating system and hardware. All runners come with pre-installed compilers, tools and software. Unfortunately, none of the HAMS-MREL dependencies come by default with the runners so we must install them ourselves. Runners are created when a workflow is 

Workflows define `triggers` that tell GitHub when the workflow needs to be executed. A workflow defines one or more `jobs`, and each job is made up of one or more `steps`. Steps within a job are executed sequentially. By default, all jobs run in parallel unless dependencies are defined between them.

## Iterative Solvers

HAMS-MREL supports two linear-solver paths for the BEM system, dispatched at runtime by the `ISOLV` parameter in `ControlFile.in` (1 = direct LU, default; 2 = GMRES). The full user-facing description is in [iterative_solvers.md](iterative_solvers.md); this section covers the developer-side layout.

**Source files:**

| File | Contents |
|---|---|
| [`src/IterativeSolvers.f90`](../src/IterativeSolvers.f90) | Restarted GMRES(m), block-diagonal LU preconditioner, Block GMRES with deflated initial guesses. Public entry points: `ZGMRES_SOLVE`, `ZBLOCK_GMRES_SOLVE`, `ZITER_SOLVE_MULTI`, `SETUP_BLOCK_PRECONDITIONER`. Uses inline BLAS (`ZDOTC`/`DZNRM2`/`ZGEMV` unrolled) to sidestep known ifort/MKL calling-convention issues with the compiled BLAS bindings. |
| [`src/HMatrix.f90`](../src/HMatrix.f90) | Hierarchical-matrix compression: cluster tree, block cluster tree with admissibility test, Adaptive Cross Approximation (ACA), fast H-matvec. Public entry points: `HMATRIX_BUILD`, `HMATVEC`, `HMATRIX_DESTROY`. |

**Dispatch pattern:**

Each of the four `AssbMatx*.f90` files gates the LAPACK path on `ISOLV`. `ASSB_LEFT*` skips the `ZGETRF` factorization when `ISOLV > 1`. Each `RADIATION_SOLVER*` / `DIFFRACTION_SOLVER*` branches: `ISOLV = 1` uses the existing in-place `ZGETRS` path (Phase 2.4 — passes `AMAT`/`CMAT` directly, no full-size workspace copy); `ISOLV = 2` allocates a scratch matrix, builds the H-matrix + block preconditioner, and calls `ZITER_SOLVE_MULTI`.

**Interaction with `GRN_KIND=4`:**

The `ACA_TOL` and `ITER_TOL` constants are both `1e-4`, matched to the ~`1e-7` per-entry precision floor of an `AMAT` assembled from single-precision Green's-function values. Iterative results agree with the direct-LU reference at 4–5 significant figures on physically meaningful entries. For bit-identity against WAMIT references, users should stay on the default (`GRN_KIND = 8`, `ISOLV = 1`).

**Adding a new iterative test:**

Small multi-body test cases with `ISOLV = 2` can be added to `test/integration_tester.f90`. The recommended check is that `ISOLV = 2` output agrees with an `ISOLV = 1` baseline to `1e-3` relative on the WAMIT `.1` (`AmssDamp.1`) file peak-magnitude entries.

## Damping lid

The damping lid holds the total potential at zero on a user-given set of panels on the free surface z = 0 (`DampingLid.pnl`, see [input_files.md](input_files.md#dampinglidpnl)). Design: [the spec](superpowers/specs/2026-11-05-damping-lid-design.md). Measurements: [test/data/damping_lid/README.md](../test/data/damping_lid/README.md). Code: `src/DampingLid.f90` (one body), `src/DampingLidMulti.f90` (several bodies). Each claim below is marked measured (source in the README) or derived.

### What the lid does

Derived: the mean power through the free surface per unit area is proportional to Im(phi* dphi/dz), which is zero where phi = 0. The lid absorbs no energy. It reflects the wave, as a pressure-release surface does, despite the name of the file. Energy balance and the Haskind relation still hold, and the tests use them.

Measured (README, "Single-body solve with a lid"): the diagonal B_ii > 0 in every single-body case (the 12 x 12 B of the two-body `lid_gap` has a smallest eigenvalue of -5.5e-5 of max(abs(B)), a hull-mesh effect, see Several bodies), and the Haskind heave residual is 0.32% at 640 lid panels. For the ring r 6 to 10 around the draft 5 cylinder, the lid cuts B11 by a factor of 100 and abs(X1) by 10.

### Formulation

Derived (spec, "Formulation"). The wave Green function obeys dG/dz = nu G on z = 0, so the free-surface integral leaves a single layer on the lid with the unknown q = dphi/dz - nu phi. With D the double layer, S the single layer, j over hull panels and l over lid panels:

    hull row i:  2 pi phi_i + sum_j D_ij phi_j - sum_l S_il q_l = sum_j S_ij V_j
    lid row l:               sum_j D_lj phi_j - sum_l' S_ll' q_l' = sum_j S_lj V_j + T_l

Radiation has T_l = 0. Diffraction solves for the scattered potential, which cancels the incident one on the lid, so T_l = 4 pi phi_I(p_l). The unknowns are hull phi and lid q. Without `IRSP` the system is square and is solved by direct LU (`ZGETRF`). With `IRSP` the interior rows make it overdetermined, and it is solved by least squares (`ZGELS`).

The sign of q is not observable in A, B or X (README, "Field output with a lid"). It rests on the derivation and on the test `lid_elevation_zero_at_centroids`: with the opposite sign the elevation at a lid centroid would be about twice the incident one, not about zero. Measured: 7.0e-8 for diffraction. The independent check of the signs and of the size of the lid's effect is the eigenfunction-matching comparison below.

### Why the IRSP rows alone fail

Derived: `IRSP` adds rows and no unknowns. Inside a body the rows state the null-field identity, which the solution satisfies. On open water the same rows state 4 pi phi(p) = 0, which the hull solution does not satisfy. The lid needs the extra column q and the incident term T_l.

Measured (spec, "Why the IRSP rows cannot be reused unchanged"; binary from 0a77224, not a build of main): `IRSP = 1` with an open-water annulus r 6 to 10 as the waterplane mesh; cylinder r = 5, draft 5, depth 50, 1.35 rad/s.

| Lid panels | A11 | B11 | A33 | B33 |
| ---: | ---: | ---: | ---: | ---: |
| none | 253.4 | 206.3 | 210.2 | 21.9 |
| 40 | 434.3 | -186.3 | 88.1 | 14.3 |
| 160 | -255.9 | -374.7 | 65.1 | 1327 |
| 640 | -337.2 | -102.6 | 602.9 | 239.1 |

The results do not converge with refinement, B11 is negative, and the runs gave no warning.

### Self term and its accuracy

Derived (spec, "Panel integrals"): the Rankine and image parts are exact panel integrals. For panels apart the wave term is one point at the centroid. For a field point inside the panel, the deep-water wave term is integrated by a Duffy rule after splitting the panel at the field point, and the finite-depth remainder G_finite - G_deep is added as a smooth term sampled at R/h of at least 0.03.

Measured (README, "Single-layer integral of G over a lid panel"):

- Duffy rule, deep water, nu = 0.19, 6 x 6 against 24 x 24 (the rule used, u = s**2 with the split): 6.4e-7 to 8.6e-7 for the square and 6.3e-7 to 7.9e-7 for the triangle over the field points tried. The plain rule gave 6.8e-5 at the square's centroid.
- One-point wave term for panels apart: the error is about (nu side)**2/24, 1.5e-3 at side 1, and does not fall with distance. Derived: the hull's one-point wave term has the same error.
- Finite depth: `FINGREEN3D` matches the series to 3e-10 or better from rho/H = 0.05. Between 0.0005 H and 0.05 H it is off by a near constant 7.9e-6, and below 0.003 H it scatters, to 1.1e-2 at 0.001 H and nu = 1. It returns NaN from kh = 360.
- Electrified disc (nu = 0, relative error of the density against 1/pi): -7.80e-2, -2.58e-2, -9.31e-3 for nt = nr = 8, 16, 32. The error falls by 3.0 and 2.8 per doubling. The density has a 1/sqrt singularity at the rim, so the ring count matters as much as the panel count.

### Convergence

Measured (README, "Single-body solve with a lid"): from 40 to 160 to 640 lid panels the change falls by 2.1 (A11), 4.6 (A33), 7.6 (B11) and 2.7 (abs(X1)). The ring count drives the results. From 4 rings on, the changes of B33 fall as 0.15, 0.084, 0.043, about first order in the ring width. Derived: first order is expected of constant panels with the 1/sqrt edge singularity of q. The approach to no lid is slow for a narrow ring (a ring 0.25 m wide still changes B11 by 83%). A likely reason, not checked, is that a phi = 0 strip acts like a conductor, whose capacity falls only as 1/log(1/width).

Conditioning, measured (README, "Condition number"): the 1-norm condition number of the assembled lid system is 231, 118 and 98 for 40, 160 and 640 lid panels (`lid_ann_*`, `NSYS = 1`, 1.35 rad/s), against 11.1 for the hull-only matrix at the same frequency. It does not grow with lid refinement. The ratio of extreme singular values is 48, 30 and 51 (hull only 4.6).

Symmetry, measured: the half hull with the half lid gives A, B and X of the full model to 5.4e-7, 6.6e-6 and 1.8e-5 (relative to the largest entry).

### Independent check: eigenfunction matching

Measured (README, "Independent reference"): a bottom-mounted cylinder (radius 2 m, depth 10 m, kh = 1.5) with an annular lid r 2 to 4, against `pre_postprocessing_python/DampingLid/annular_lid_reference.py`, which shares no code with the solver. The reference converges as 1/N, and its projections agree to 8e-5. The relative error of HAMS-MREL in X1 is 3.14%, 1.07% and 0.28% for 2, 4 and 8 lid rings (in eta(8, 0, 0): 6.54%, 2.95%, 1.48%). The reference gives abs(X1) 0.5445 of the bare cylinder and HAMS-MREL 0.5499 at 8 rings. The phase shift of X1 caused by the lid agrees to 0.3 deg. This checks the signs and the size of the lid's effect. It does not check the factor of 100 on B11 of the `lid_ann_*` cases, which have another geometry.

### Several bodies

One lid file serves all bodies, since the lid belongs to no body. Measured (README, "Several bodies"): with the second body 5 km away, body 1 differs from the single-body lid result by 5.0e-4 at most in A11, A33, B11, B33 and abs(X1). The elevation on the lid between two boxes (`lid_gap`) is 9.4e-15 at most. The diagonal of B of `lid_gap` is positive, but the full 12 x 12 B has a smallest eigenvalue of -5.5e-5 relative to max(abs(B)); the same boxes without the lid have -2.3e-5. Derived: these come from the hull discretisation, since they fall by 4 to 9 per halving of the panel size.

### Limits

- No GMRES: `Solver_type 2` with a lid stops (`damping lid needs Solver_type 1`). Derived: the H-matrix assumes hull geometry, and the lid rows would need their own preconditioner.
- `IRSP` with a single-body lid runs but is experimental: the least-squares solve holds phi = 0 on the lid only to about 1e-5 (Open questions). With several bodies it stops (`damping lid with several bodies needs If_remove_irr_freq 0`), because the multi-body interior rows are not combined with the lid.
- The lid reflects and absorbs no energy (above).
- Finite depth: every integral over a lid panel samples `FINGREEN3D` at R/h of at least 0.03, so lid-to-lid pairs are protected. Measured on z = 0 against the series (README, "Remainder"): `FINGREEN3D` is off by a near constant 7.9e-6 between 0.0005 h and 0.05 h, and below 0.003 h it scatters, to 6.7e-5 at 0.003 h and 1.1e-2 at 0.001 h (nu = 1). Derived: the exposed paths are the hull panels in the lid rows (`HULL_LAYERS` calls `FINGREEN3D` at the true distance) and in the field output, so lid centroids or field points within about 0.01 h of a hull panel centroid get the scatter. Hull centroids lie below z = 0, where it was not measured.
- The inside-hull check (`CHECK_LID`) tests lid centroids only, so a lid panel that straddles the waterline with its centroid outside passes. The in-panel test of `LID_DEEP_WAVE` (the signs of edge cross products) assumes convex panels.
- Accuracy is set by the number of lid rings, first order in the ring width (measured, above). The solver does not flag an under-resolved lid.
- Field points at z = 0 get radiation elevation rows in `PressureElevation.6p` only when the lid is on.

### Open questions

- IRSP with a single-body lid (experimental, see [input_files.md](input_files.md#dampinglidpnl)): against the same run without IRSP, A11 changes by 0.09% and B33 by 0.06%, but B11 and abs(X1) change by 1.1% for `lid_ann_160`. The Haskind heave residual with IRSP grows with the lid mesh, 4.4e-3 at 160 and 8.1e-3 at 640 panels, against 4.6e-4 and 3.1e-3 without IRSP (measured, README "Symmetry and IRSP"). The least-squares solve does not hold phi = 0 on the lid as the square solve does. The diffraction elevation at the lid centroids is 1.6e-5 to 1.8e-5 with IRSP at 160 lid panels and 1.2e-5 at most at 640, against 9.3e-9 and 7.0e-8 without IRSP. For radiation it is 2.0e-5 and 1.3e-5 at most, against 1.1e-8 and 1.5e-7 (measured, README "Symmetry and IRSP"). The cause is open. Hypothesis, not checked: the lid rows have a different scale from the hull and interior rows, so the least-squares fit puts part of its residual on them. Derived and unproven as well: a discretisation difference of the least-squares system. Checks that would settle it: hulls finer than 390 panels with IRSP, and a least-squares solve with the rows scaled to equal norms.
- A multi-body lid with symmetry (`ISX` or `ISY`) is untested.
- Why the elevation at lid centroids is about 1e-8 in the annular cases and about 1e-14 in `lid_gap` was not looked for.

## Profiling Results

This section documents the primary sources of memory consumption in HAMS-MREL and how they scale with problem size. It is intended to guide future development and optimization efforts aimed at reducing RAM usage.

Profiling of HAMS-MREL shows that memory usage is dominated by a small number of large, multi-dimensional arrays whose sizes scale with the number of bodies and the number of panels per body. These arrays are dynamically allocated before the main solver iteration loop begins and remain allocated for the duration of the simulation.

The table below lists the highest-memory-consumption variables and their sizes for systems involving 2, 3, and 6 bodies. For all cases shown, each body is discretized using 234 panels. For all cases shown, each body is discretized using 234 panels.

| System | Variable Name | Type | Size | Number of Elements |
| :--- | :--- | :--- | :--- | :--- |
| 2-bodies (234 panels each) | `CGRN_MULTI_COMB` & `RKBN_MULTI_COMB` | `COMPLEX(8)` | (1:468,1:468,1,4) | 876096 |
| 2-bodies (234 panels each) | `AMAT_MULTI` | `COMPLEX(8)` | (1:468,1:468,1:1) | 219024 |
| 3-bodies (234 panels each) | `CGRN_MULTI_COMB` & `RKBN_MULTI_COMB` | `COMPLEX(8)` | (1:702,1:702,1:1,1:4) | 1971216 |
| 3-bodies (234 panels each) | `AMAT_MULTI` | `COMPLEX(8)` | (1:702,1:702,1:1) | 492804 |
| 6-bodies (234 panels each) | `CGRN_MULTI_COMB` & `RKBN_MULTI_COMB` | `COMPLEX(8)` | (1:1404,1:1404,1:1,1:4) | 7884864 |
| 6-bodies (234 panels each) | `AMAT_MULTI` | `COMPLEX(8)` | (1:1404,1:1404,1:1) | 1971216 |

### Observations

- **Dominant memory consumers**
  - The first two dimensions of these arrays scale linearly with the total number of panels across all bodies, resulting in a quadratic growth in memory usage as the number of bodies increases.
  - The variables are allocated once per program execution, before the solver iteration begins.
- **Overall memory behavior**
  - For a given simulation, memory usage is largely constant throughout execution.
  - Peak memory consumption is approximately 2 GB for the cases studied.
  - Memory does not accumulate across solver iterations.



