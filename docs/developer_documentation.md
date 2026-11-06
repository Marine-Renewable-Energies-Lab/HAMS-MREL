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

The damping lid is the dissipative lid of WAMIT (WAMIT User Manual v7.5, Section 12.8 and eqs. 15.91 to 15.93, `IDAMPER = 1`) on a user-given set of panels on the free surface z = 0 (`DampingLid.pnl`, see [input_files.md](input_files.md#dampinglidpnl)). Design: [the spec](superpowers/specs/2026-11-06-epsilon-damping-lid-design.md). Measurements: [test/data/damping_lid/README.md](../test/data/damping_lid/README.md). Code: `src/DampingLid.f90` (one body), `src/DampingLidMulti.f90` (several bodies), the field terms in `src/PressureElevation.f90` and `src/PressureElevationMulti.f90`. Each claim below is marked measured (source in the README) or derived.

This section replaces the one on the lid that held phi = 0 on its panels ([its spec](superpowers/specs/2026-11-05-damping-lid-design.md), superseded). Derived (spec, "Why the phi = 0 lid is replaced"): phi = 0 makes the mean energy flux through the lid, proportional to Im(phi* dphi/dz), zero, so that lid reflected the wave and absorbed nothing. The phi = 0 measurements are in the git history of this file and of the README (commit 27cb37e and earlier); the panel-integral and hull measurements below do not depend on the lid condition and are kept.

### What the lid does

Derived (spec, "Formulation"). On the lid the free surface carries the pressure p/rho = epsilon dphi/dz, epsilon in m/s. With the time factor e^{-i omega t} and nu = omega**2/g, the dynamic and kinematic conditions give

    dphi/dz = beta phi,    beta = nu g/(g - i omega epsilon)

The power absorbed by the lid is P = 1/2 rho sum_l epsilon_l abs(dphi/dz)_l**2 A_l >= 0. epsilon = 0 gives beta = nu, the open free surface; epsilon to infinity gives dphi/dz = 0, a rigid lid.

Measured (README, "Single-body solve with a lid", `lid_ann_*`: cylinder r = 5, draft 5, lid r 6 to 10, 50 m of water, 1.35 rad/s):

- epsilon = 0 gives the results without a lid to the last printed digit, for one body and for two (`lid_gap_eps0`, 5.0e-16 in A, B, X).
- The lid power from the centroid elevations is positive: 4.84e4 W per m**2 of incident amplitude in diffraction, and 4.8% (surge) and 4.5% (heave) of the power the body radiates.
- B11 falls from 206.25 without the lid to 199.71 at 0.2 m/s and 48.89 for the rigid lid (1e6 m/s). Derived: B is the far-field power plus the lid power, and the lid can cut the far-field radiation by more than it absorbs, so B can fall.
- At 0.2 m/s and this frequency the lid is weak: omega epsilon/g = 0.028, and A, B, X move by about 3% from no lid.

### Formulation

Derived (spec, "Formulation"). The wave Green function obeys dG/dz = nu G on z = 0, so the free-surface integral of Green's identity leaves only the lid, with q = dphi/dz - nu phi = alpha phi and

    alpha = beta - nu = i omega epsilon nu/(g - i omega epsilon)

The unknowns are phi on the hull panels and phi on the lid panels. D and S are the double- and single-layer panel integrals, j runs over hull panels, l over lid panels, and alpha_l uses the epsilon of panel l:

    hull row i:  2 pi phi_i + sum_j D_ij phi_j - sum_l alpha_l S_il phi_l = sum_j S_ij V_j
    lid row l:   4 pi phi_l + sum_j D_lj phi_j - sum_l' alpha_l' S_ll' phi_l' = sum_j S_lj V_j

The lid row has the 4 pi free term of WAMIT eq. 15.93. Diffraction solves for the scattered potential. The incident potential obeys dphi_I/dz = nu phi_I on z = 0, so the scattered potential has q_S = alpha (phi_S + phi_I) on the lid, and every row gains + sum_l alpha_l S_xl phi_I,l on the right-hand side. The system is square and is solved by direct LU (`ZGETRF`); with `IRSP` the interior rows make it overdetermined, and it is solved by least squares (`ZGELS`).

Field points: 4 pi phi(x) = sum_j (S_xj V_j - D_xj phi_j) + sum_l alpha_l S_xl (phi_l + phi_I,l), with phi_I,l in diffraction only. The elevation on a lid panel l is eta = i omega phi_T/(g - i omega epsilon_l), phi_T the total potential, and i omega phi_T/g elsewhere.

Symmetry: the lid columns of symmetry class IP are sum_IS RXY(IS,IP) S(x reflected for IS, l) times alpha_l, and the incident term uses the class components of phi_I at the reflected lid centroids. The images of a panel take its epsilon.

Limits that follow (derived): at epsilon = 0, alpha = 0, so the lid rows only evaluate the representation and do not feed the hull rows. At V = 0 alpha = 0 as well, so at the zero frequency of `0_inf_frequency_limits` the lid is skipped, as at the infinite frequency. A and B stay symmetric by reciprocity, and B is positive semi-definite. The Haskind relation over the hull does not hold for epsilon > 0, since the lid adds sum_l alpha_l phi_j phi_I A_l.

Sign of alpha, measured (README, "Independent reference"): HAMS-MREL at 1 m/s is 8.4e-3 from the eigenfunction reference and 5.4e-2 from the reference computed with the conjugate lid condition (Im(alpha) flipped). A binary built with the sign flipped fails five integration tests. The reference's own energy balance (net power into r = 1000 m equals the lid power, to 2.6e-10) shows that the reference has the absorbing sign.

### Why the IRSP rows alone fail

Derived: `IRSP` adds rows and no unknowns. Inside a body the rows state the null-field identity, which the solution satisfies. On open water the same rows state 4 pi phi(p) = 0, which the hull solution does not satisfy. The lid needs its own unknowns and lid term.

Measured (spec of the phi = 0 lid, "Why the IRSP rows cannot be reused unchanged"; binary from 0a77224, not a build of main): `IRSP = 1` with an open-water annulus r 6 to 10 as the waterplane mesh; cylinder r = 5, draft 5, depth 50, 1.35 rad/s.

| Lid panels | A11 | B11 | A33 | B33 |
| ---: | ---: | ---: | ---: | ---: |
| none | 253.4 | 206.3 | 210.2 | 21.9 |
| 40 | 434.3 | -186.3 | 88.1 | 14.3 |
| 160 | -255.9 | -374.7 | 65.1 | 1327 |
| 640 | -337.2 | -102.6 | 602.9 | 239.1 |

The results do not converge with refinement, B11 is negative, and the runs gave no warning.

### Self term and its accuracy

Derived (spec of the phi = 0 lid, "Panel integrals"; unchanged by the epsilon lid): the Rankine and image parts are exact panel integrals. For panels apart the wave term is one point at the centroid. For a field point inside the panel, the deep-water wave term is integrated by a Duffy rule after splitting the panel at the field point, and the finite-depth remainder G_finite - G_deep is added as a smooth term sampled at R/h of at least 0.03.

Measured (README, "Single-layer integral of G over a lid panel"):

- Duffy rule, deep water, nu = 0.19, 6 x 6 against 24 x 24 (the rule used, u = s**2 with the split): 6.4e-7 to 8.6e-7 for the square and 6.3e-7 to 7.9e-7 for the triangle over the field points tried. The plain rule gave 6.8e-5 at the square's centroid.
- One-point wave term for panels apart: the error is about (nu side)**2/24, 1.5e-3 at side 1, and does not fall with distance. Derived: the hull's one-point wave term has the same error.
- Finite depth: `FINGREEN3D` matches the series to 3e-10 or better from rho/H = 0.05. Between 0.0005 H and 0.05 H it is off by a near constant 7.9e-6, and below 0.003 H it scatters, to 1.1e-2 at 0.001 H and nu = 1. It returns NaN from kh = 360.
- Electrified disc (nu = 0, relative error of the density against 1/pi): -7.80e-2, -2.58e-2, -9.31e-3 for nt = nr = 8, 16, 32. The error falls by 3.0 and 2.8 per doubling. This is a test of S alone, with the 1/sqrt singularity of a phi = 0 strip at the rim.

### Convergence

Measured (README, "Single-body solve with a lid"): from 40 to 160 to 640 lid panels of the `lid_ann_*` cases at 0.2 m/s, the change falls by 3.1 (A11), 4.0 (A33), 8.5 (B11), 4.3 (B33) and 8.6 (abs(X1)). A narrowing ring approaches no lid about in proportion to its width (distance 5.5e-3 at 0.5 m, 2.8e-3 at 0.25 m).

Lid convergence against the eigenfunction reference, measured (README, "Lid convergence with the 48 x 20 hull"; bottom-mounted cylinder of the next section, 48 x 20 hull, lid r 2 to 4 with 48 panels around). The worst lid-effect error is 8.15e-3, 7.97e-3 and 7.90e-3 for 8, 16 and 32 rings at 0.2 m/s, and 8.10e-3, 7.88e-3 at 8 and 32 rings at 1 m/s, while the raw error stays at the bare-hull error (2.1e-3 in X1). Derived: the rings are resolved at 8. What remains is set by the panels around and by the hull rows near the waterline: 96 lid panels around give 5.81e-3, and the 72 x 30 hull with a 72-around lid 3.71e-3. The floor falls with finer circumferential and hull resolution; the rate, and the split between the two, are not determined by these runs (doubling the panels around alone gives a factor of 1.4 to 1.6). The lid power error is of the order of the area deficit of the polygonal lid (-3.2e-3 measured against 2.9e-3 at 48 panels around, -1.1e-3 against 7.1e-4 at 96).

Symmetry, measured (README, "Symmetry and IRSP"): the half hull with the half lid (`ISY = 1`) gives A, B and X of the full model to 1.4e-6 relative to the largest entry, and the field values to 6.7e-7.

Conditioning: measured with the phi = 0 lid only (README, "Condition number"), 1-norm condition number 98 to 231 against 11.1 for the hull alone. Not repeated for the epsilon lid. Derived, not measured: its lid block 4 pi I - alpha S is close to the identity for moderate epsilon.

### Independent check: eigenfunction matching

Measured (README, "Independent reference"): a bottom-mounted cylinder (radius 2 m, depth 10 m, kh = 1.5) with an annular lid r 2 to 4, against `pre_postprocessing_python/DampingLid/annular_lid_reference.py`, which shares no code with the solver. In the lid region the vertical modes satisfy k tanh(kh) = nu g/(g - i omega epsilon), with complex roots by Newton iteration. The reference passes seven self-checks, among them MacCamy-Fuchs through the lid modes at epsilon = 0 (3.0e-16), the rigid lid by a separate code path (1.2e-6) and the energy balance (2.6e-10 or less). Its values are extrapolated in 1/N with a spread of 1.5e-7 or less at 0.2 and 1 m/s.

HAMS-MREL with the 24 x 10 hull and the lid of 48 x 8 panels, complex relative error in X1, eta(0, 8, 0), eta(8, 0, 0):

| epsilon (m/s) | X1 | eta(0, 8, 0) | eta(8, 0, 0) | Lid effect, worst | Lid power, relative error |
| ---: | ---: | ---: | ---: | ---: | ---: |
| none (against MacCamy-Fuchs) | 8.50e-3 | 6.53e-4 | 8.90e-4 | | |
| 0.2 | 8.49e-3 | 6.46e-4 | 9.10e-4 | 2.21e-2 | -4.14e-3 |
| 1 | 8.42e-3 | 6.19e-4 | 9.79e-4 | 2.20e-2 | -4.12e-3 |
| 1e6 | 7.07e-3 | 6.75e-4 | 7.83e-4 | 1.97e-2 | |

The lid effect is the change the lid makes in HAMS-MREL against the change it makes in the reference, relative to the reference change; the hull error drops out of it to first order. Derived from the table: the error against the reference is the hull error of the bare cylinder (8.5e-3 with and without the lid), so it cannot show lid convergence; the lid effect can. Of the -0.41% in lid power, about 0.29% is the polygonal lid area (README).

### Several bodies

One lid file serves all bodies, since the lid belongs to no body. Measured (README, "Several bodies"):

- With the second body 5 km away, body 1 differs from the single-body lid result by 4.5e-4 at most in A11, A33, B11, B33 and abs(X1).
- Two boxes 20 x 8 x 4 m with a 4 m gap and a lid over the whole gap (`lid_gap`, heading 90 degrees, 1.35 rad/s): A and B are symmetric to 8.0e-3, and the smallest eigenvalue of (B + B^T)/2 is -4.6e-6 of max(abs(B)) (-2.3e-5 without the lid). Derived: these come from the hull discretisation, as measured with the phi = 0 lid, where they fell by 4 to 9 per halving of the panel size.
- The largest gap elevation in diffraction falls with epsilon: 2.967 without the lid, 2.430 at 0.2 m/s, 1.306 at 1 m/s.

Measured (README, "Gap resonance"): a sweep of `lid_gap_nolid` from 0.8 to 2.0 rad/s puts the epsilon = 0 gap resonance at 1.29 to 1.30 rad/s (largest gap elevation 3.89), so 1.35 rad/s is near it, at 76% of the peak. At 1.29 rad/s the largest gap elevation is 3.893 without the lid, 3.048 at 0.2 m/s and 1.560 at 1 m/s.

### Limits

- epsilon must be chosen by the user and calibrated against model tests or CFD; the code has no default. Starting values are in [input_files.md](input_files.md#dampinglidpnl).
- A, B, X and the RAOs are integrated over the hull only, as in WAMIT (WAMIT User Manual v7.5, Section 12.8). The lid power is not written.
- The Haskind relation over the hull does not hold for epsilon > 0 (derived above; measured heave residual 4.6e-2 at 0.2 m/s against 2.1e-3 without the lid).
- A zero frequency in the frequency list stops a run with the lid on (`zero frequency with a damping lid`). The zero frequency of `0_inf_frequency_limits` is exempt, and the lid is skipped there.
- No GMRES: `Solver_type 2` with a lid stops (`damping lid needs Solver_type 1`). Derived: the H-matrix assumes hull geometry, and the lid rows would need their own preconditioner.
- `IRSP` with a single-body lid runs but is experimental: it changes `lid_ann_160` by 0.32% at most (B11). With several bodies it stops (`damping lid with several bodies needs If_remove_irr_freq 0`), because the multi-body interior rows are not combined with the lid.
- Finite depth: every integral over a lid panel samples `FINGREEN3D` at R/h of at least 0.03, so lid-to-lid pairs are protected. Derived: the hull panels in the lid rows (`HULL_LAYERS` calls `FINGREEN3D` at the true distance) and in the field output are not, so lid centroids or field points within about 0.01 h of a hull panel centroid get the scatter measured above.
- The inside-hull check (`CHECK_LID`) tests lid centroids only, so a lid panel that straddles the waterline with its centroid outside passes. The in-panel tests of `LID_DEEP_WAVE` and `LID_EPS_AT` (the signs of edge cross products) assume convex panels.
- The solver does not flag an under-resolved lid.
- Field points at z = 0 get radiation elevation rows in `PressureElevation.6p` only when the lid is on.

### Open questions

- No WAMIT or Gourlay geometry is available, so there is no cross-code comparison of the epsilon lid.
- IRSP with a single-body lid: with the phi = 0 lid IRSP changed B11 and abs(X1) of `lid_ann_160` by 1.1%, and the cause was not found (README, "Symmetry and IRSP"). With the epsilon lid at 0.2 m/s the change is 0.32%, but that lid is close to no lid, where IRSP changes `image_cos` by 0.1 to 0.3%, so this does not settle the question.
- A multi-body lid with symmetry (`ISX` or `ISY`) is untested.
- The condition number of the epsilon-lid system was not measured.

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



