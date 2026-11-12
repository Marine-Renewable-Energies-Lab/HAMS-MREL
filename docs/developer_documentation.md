# Development Documentation

This page contains information for the developers and mantainers of HAMS-MREL:

- [Tests](https://github.com/vaibhavraghavan/HAMS-MREL/edit/develop/docs/developer_documentation.md#tests)
- [Continous Integration](https://github.com/vaibhavraghavan/HAMS-MREL/edit/develop/docs/developer_documentation.md#continous-integration-in-github-actions)
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

### Mean drift tests

The mean drift code (`SecondOrder*.f90`, `FarFieldDrift.f90`, `NearFieldDrift.f90`, `ControlSurfaceDrift.f90`, `ControlSurfaceMesh.f90`, `FieldVelocity.f90`, `SurfaceVelocity.f90`, `SurfaceGeometry.f90`, `WaterlineField.f90`, `FirstOrderState.f90`, `QTFDriver.f90`, `QTFPotential.f90`, `QTFField.f90`, `SecondOrderIncident.f90`, `SharpEdge.f90`) is covered by unit tests in `test/test_*.f90` and by integration cases in `test/test_integration_drift.f90`, whose fixtures and measured values are in [test/data/drift_cyl/README.md](../test/data/drift_cyl/README.md). Usage and the method choice are in [SecondOrder.in](input_files.md#secondorderin-optional), and the studies behind them in [Second-order validation](second_order_validation.md).

With the QTF potential part (Phase 5) and the partial second-order field (Phase 5b, `qtf_field`), the unit tester has 250 tests and 1 expected failure, and the integration tester 183, all passing with ifort 2022.1 on Windows (clean build of the merge, 2026-11-11). Before the merge, with the potential part alone, they were 226 and 169 (see [QTF potential part](#qtf-potential-part)). With the quadratic QTFs (Phase 4), the unit tester had 185 tests and the integration tester 147, all passing with ifort 2022.1 on Windows. With the control-surface method, before them, it was 177 and 141. At commit dff7a32, before the control-surface method, they had 133 and 118, all passing with ifort 2022.1 on Windows and with ifx 2024.2 on Linux; the times below are from then. Measured run time: the integration tests take 148 s with ifort on a laptop with 8 logical cores and 177 s with ifx at 32 threads on Snellius; the unit tests take 19 s with ifx. All cases stay in CI. The time on the GitHub runners is not measured here; GitHub's default job limit is 6 hours.

`QTFDriver.f90` keeps, for the QTFs, the converted state of every frequency and heading: φ_I, φ_P, ∇φ_I, ∇φ_P (8 complex values per panel, 128 bytes) and the motions, and the waterline potential. The whole store is allocated by `qtf_init` at set-up, so a run that does not fit stops before the first-order solves. `qtf_keep` must be called for every heading of a frequency, in order (heading 1 starts a new frequency), and only for the frequencies of the QTF grid, which leaves out the limits of `SYBO = 1`; it stops if it is called out of order. Its state is module state, so the unit tester runs the `qtf` suite serially.

### QTF potential part

The body forcing of the second-order potential in the difference QTF (`qtf_terms = body`, Phase 5). Usage is in [QTFs, potential part](input_files.md#qtfs-potential-part), the files in [Output Files](output_files.md#wamit-format), and the validation in [Second-order validation](second_order_validation.md#qtf-potential-part).

**Modules.**

| File | Contents |
| :--- | :--- |
| [`src/QTFPotential.f90`](../src/QTFPotential.f90) | The body-forcing terms of one pair of states (`potential_pair`): Froude-Krylov and diffraction of the second-order incident wave, B1, B2, B3. ψ at the waterline (`psi_waterline`), its tangential gradient on the hull (`psi_surface_gradient`), u at the waterline (`waterline_velocity`), and the control-surface form of B2 (`b2_body_term`, `b2_control_term`, `b2_control_check`). |
| [`src/QTFDriver.f90`](../src/QTFDriver.f90) | The store of the first-order states (Phase 4), and for the potential part: the list of difference frequencies (`qtf_radiation_list`, `qtf_radiation_omegas`), ψ of every Ω (`qtf_keep_psi`), the fields on C and F (`qtf_set_control_surface`, `qtf_keep_cs_velocity`, `qtf_keep_cs_psi`), the edge weights (`qtf_set_edge_weights`), the potential part Ω by Ω (`potential_part`) and the files (`qtf_finish`). Counters of NaN or Inf in ψ and T_pot set the exit code. |
| [`src/SecondOrderIncident.f90`](../src/SecondOrderIncident.f90) | The second-order incident potential φ_I2 of a pair of waves at Ω = ω_i − ω_j, any depth and headings (Molin 1979; Sharma and Dean 1981). Zero on the exact diagonal. |
| [`src/SharpEdge.f90`](../src/SharpEdge.f90) | Sharp edges of the hull and the product weights of the panels next to them, for `qtf_edge_weights = on` (off by default). |
| [`src/SecondOrderDriver.f90`](../src/SecondOrderDriver.f90) | Set-up (`body_terms_setup`, `b2_control_setup`), ψ after each radiation solve (`so_after_radiation`) and after each radiation-only solve (`so_after_radiation_only`), ∇ψ on C and F (`keep_cs_psi`), and `QTF_PsiCheck.txt`. |
| [`src/HAMS_Prog.f90`](../src/HAMS_Prog.f90) | `RadiationOnlySolves`: the radiation problems at the difference frequencies that are not on the main grid, after the frequency loop. |

**Formulation.** Time factor e^{−iωt}, normals out of the fluid, per unit wave amplitudes. At Ω = ω_1 − ω_2, with the free-surface integral left out,

T_pot,k = iΩρ [ ∬_S φ_I2 n_k dS + ∬_S ψ_k (B − ∂φ_I2/∂n) dS ],
B = (Ẋ − u)·(α × n) [B1] − n·[(X·∇)u] [B2] + n·d/dt(H x) [B3], H x = ½ α × (α × x),

with X = ξ + α × (x − x_r) the first-order motion, u the total first-order velocity and ψ_k the radiation potential of mode k at Ω (∂ψ_k/∂n = n_k). This is the indirect (Haskind) method of Molin (1979) with the second-order body condition on the mean surface of Ogilvie (1983); see also WAMIT v7.5 (ch. 15) and Chen (1988). Moments come from ψ_4 to ψ_6, so no moment arm enters. Every product of two first-order quantities uses the difference kernel of `NearFieldDrift` (`pair_product`, `pair_cross`), as the quadratic part does. For Ω < 0, ψ(−|Ω|) = ψ(|Ω|)* is used, so T_pot(2, 1) = T_pot(1, 2)*.

**B2.** The direct integrand n·[(X·∇)u] is not integrable at a sharp edge, where u ~ ρ^(−1/3). Stokes' theorem on S, bounded by the waterline Γ, gives

∬_S ψ n·[(X·∇)u] dS = R + S_B2 + L, R = ∬_S ψ n·(α × u) dS, S_B2 = ∬_S (n × ∇ψ)·(X × u) dS, L = −∮_Γ ψ (X × u)·t dl,

so B2 = −(R + S_B2 + L), and only first derivatives of φ and the tangential gradient of ψ enter ([QTFPotential.f90](../src/QTFPotential.f90)). S_B2 has a ρ^(−2/3) integrand at the edge and converges slowly with one point per panel. With `qtf_b2 = control` it is moved off the hull. For a rigid field A and harmonic φ, ψ (v = ∇φ), the tensor T_ab = v_a ∂_bψ + v_b ∂_aψ − δ_ab (v·∇ψ) is divergence-free and div(T A) = 0; this extends the Gauss identity of Chen (2007, App. C, eq. 63 to 65) to two potentials, as Cong et al. (2020, eq. 35 to 39) use it for the second-order body condition. On the fluid between the hull S, a closed control surface C (wall and bottom) and the free surface F between the waterline and C:

S_B2 = H + CF_T, H = ∬_S ψ_n (A·v) dS, CF_T = ∬_{C∪F} [ (A·v) ψ_no + (A·∇ψ) v_no − (A·n_o)(v·∇ψ) ] dS,

with n_o out of the control volume on C and F (the opposite of `cs_points%nrm`). No line integral enters. H needs no ∇ψ and is singular like u only; ∇ψ on C and F is smooth. R and L stay on the hull. One rigid A per control volume, so each C must enclose its own body only, and must have a bottom: a C that reaches the sea bed would leave the sea bed as a boundary of the volume, whose term is not implemented. Moving L to C ∪ F as well, and R, was tried on the exact reference and not kept ([validation](second_order_validation.md#tried-and-not-kept)).

**Fields on C and F.** `fields_green` needs the solver globals at the state's frequency. So the total velocity on C and F is evaluated at the end of each kept frequency (`so_end_frequency`, once for the control-surface drift and B2), and ∇ψ_k when ψ is taken at each Ω (`keep_cs_psi`, one pseudo-state per mode with no incident wave and ∂ψ_k/∂n = n_k on its body). Both go to the QTF store.

**Threads.** `fields_green` ([FieldVelocity.f90](../src/FieldVelocity.f90)) runs its loop over the points in parallel: each point writes only its own values, with no reduction, so the result does not depend on the thread count (unit test `green_field_thread_invariant` in `test/test_field_velocity.f90`). The panel loop of `field_influence_row` inside it is also an OpenMP loop and runs serially there while nested parallelism is off, the default. With nested parallelism on, the threads multiply (up to `NTHREAD`²); nothing guards against that. Measured on the 448-panel fixture with 686 points on C and F, one pass: about 150 s on 1 thread, 45 s on 4 threads (55 s before the point loop was threaded) and 28 s on 8; the outputs did not change.

**Notes for development.**

- HAMS's panel centroids are vertex averages. On strongly trapezoidal panels they leave a large O(h²) constant in a nearly cancelling closed-surface integral: on a sphere of 2048 panels, the Stokes form of the B2 integral differed from the direct form by 1.9e-2 (relative) with vertex-average centroids, against 5.5e-4 with area centroids (`test/test_qtf_potential.f90`). The unit test uses area centroids. Whether the vertex-average centroids of the cylinder meshes' bottom rings add a constant of this kind to their results was not measured.
- u at the waterline: the line term L with u extrapolated along a differenced panel is only O(h). `waterline_velocity` builds u at each waterline segment from the incident wave (exact), the free-surface condition (u_z), the body condition (normal part) and differencing of φ along the waterline (tangential part); B2 is then O(h²) (5.42e-3, 1.48e-3, 3.87e-4, ratios 3.66 and 3.83, `test_qtf_potential.f90`). This was measured on a regular cylinder only; an irregular spacing of the waterline segments was not measured.
- The tangential differencing of φ along the waterline runs across sharp waterline corners, which makes it O(h) next to the corners of a box barge [derived, not measured]. Flared (sloped) waterline panels are not covered by a test.
- The diagonal test is exact equality of the two frequencies, so the callers pass the same grid value for both states.
- The unit tests pin the signs of the terms as specified; the physical signs are checked only against the exact reference (`truncated_cylinder_reference.py qtf-potential`), whose signs are derived from the boundary-value problem.
- After a change to the layout of a type in `SecondOrderTypes` (such as `surface_geom`), delete the stale objects of fpm's build before testing.

**Tests.** With the potential part, the unit tester has 226 tests and 1 expected failure, and the integration tester 169 [measured on this branch at 4436f34, ifort 2022.1, Windows]. The unit tests of the potential part are in `test/test_second_order_incident.f90`, `test/test_qtf_potential.f90` and `test/test_sharp_edge.f90`; the integration cases (`qtf_body_*`, `qtf_*_body`, `qtf_round_body_irsp`, `qtf_b2c_*`, `qtf_sharp_b2_*`, `qtf_round_low*`) in `test/test_integration_drift.f90`.

Sources: Molin, B. (1979). Second-order diffraction loads upon three-dimensional bodies. Appl. Ocean Res. 1(4), 197-202. Ogilvie, T. F. (1983). Second-order hydrodynamic effects on ocean platforms. Int. Workshop on Ship and Platform Motions, Berkeley. Chen, X.-B. (1988). Etudes des réponses du second ordre d'une structure soumise à une houle aléatoire. PhD thesis, Université de Nantes. Sharma, J. N., Dean, R. G. (1981). Second-order directional seas and associated wave forces. SPE J. 21(1), 129-140. WAMIT v7.5 user manual, ch. 15. Chen, X.-B. (2007). Middle-field formulation for the computation of wave-drift loads. J. Eng. Math. 59, 61-82. Cong, P., Teng, B., Chen, L., Gou, Y. (2020). A novel solution to the second-order wave radiation force on an oscillating truncated cylinder based on the application of control surfaces. Ocean Eng. 204, 107278.

Two limits have thin margins: `source_two_body_total` measures 4.69% against 5% on both compilers, and `source_near_far_low_frequency` 4.10% (ifx) and 4.11% (ifort) against 5%. If one of them fails after a compiler or flag change, look at the printed value before changing the limit.

Do not pass an array constructor of a derived type with allocatable components, such as `[s1, s2]` of `first_order_state`, as an argument. Under ifx 2024.2 with `-heap-arrays 0` it crashes at run time; copy the values into a named array first, as `far_field_pair` does. With that change both testers pass under ifx 2024.2 with `-heap-arrays 0` on Linux (Snellius job 28933784, 133 unit and 118 integration tests).

After executing these commands you should see the following output (or similar):

<img width="557" height="280" alt="hams-tests" src="https://github.com/user-attachments/assets/1bb39a74-593d-498b-8a6d-40f20fcc36f9" />

### Direct method (QTFDirect.f90)

The direct solution of the difference-frequency second-order potential with body forcing (`qtf_direct = 1`, Phase 7a). Usage is in [QTFs, potential part](input_files.md#qtfs-potential-part), the file in [Output Files](output_files.md#wamit-format), and the validation in [Second-order validation](second_order_validation.md#direct-method-against-indirect-method).

**Formulation.** At Ω = ω_1 − ω_2, φ_D solves Laplace's equation with the bottom and radiation conditions, −Ω²φ_D + g ∂φ_D/∂z = 0 on z = 0 (no free-surface forcing in 7a), and ∂φ_D/∂n = b on S with b = B − ∂φ_I2/∂n. B is the body forcing of the [potential part](#qtf-potential-part). HAMS's integral equation for a given normal derivative is the one of `ASSB_LEFT`, so φ_D = AMAT⁻¹ r with r_i = ∬_S G(x_i, ξ) b(ξ) dS_ξ, the same right side as `ASSB_RBC` forms for ∂φ/∂n = n_k with b in place of n_k. The load is

T_dir,k = iΩρ ∬_S (φ_I2 + φ_D) n_k dS.

For the same B, Green's theorem on φ_D and ψ_k (both harmonic in the fluid, with the same free-surface and radiation conditions at Ω) gives ∬_S φ_D n_k dS = ∬_S ψ_k ∂φ_D/∂n dS, so T_dir,k = T_pot,k of the potential part [derived by Green's second identity; the indirect form is Molin 1979]. Nemoh's `.12d` (direct) and `.11d` (indirect) are the same pair, and the design spec (`notes/specs/2026-11-04-second-order-waves-design.md`, Nemoh outputs) requires that they converge to each other. That is what the `.12d` against `.11d` check tests.

**Modules.**

| File | Contents |
| :--- | :--- |
| [`src/QTFDirect.f90`](../src/QTFDirect.f90) | The right side of the B2 part (`direct_rhs_b2`, `direct_rhs_b2_block`, `b2_stokes_densities`, `b2_stokes_at_point`), the panel integrals of G and ∇_ξ G (`green_panel_integrals`, `stored_panel_integrals`, `rankine_moments`, `green_line_integral`), the sharp-edge mask (`sharp_edge_panels`), the solve (`direct_setup`, `direct_factor`, `direct_solve`) and the load of all pairs of one Ω (`direct_at_omega`). |
| `ASSB_RHS_FIELD` and `ASSB_RHS_FIELD_IRR` ([AssbMatx.f90](../src/AssbMatx.f90) and the multi-body and lid versions) | ∬ G b dS for a panel field b, per symmetry image (`V(NELEM, NRHS, NSYS)`), as `ASSB_RBC` and `ASSB_RBC_IRR` do for n_k. |
| [`src/QTFDriver.f90`](../src/QTFDriver.f90) | The store of T_dir (`tdir`, filled by `qtf_direct_at_omega`), `tdir_kept` per Ω, `QTF_Direct.12d` and `QTF_DirectTerms.txt` (`write_direct_terms`). The row writer `write_sum` is the one of `.11d`. |
| [`src/HAMS_Prog.f90`](../src/HAMS_Prog.f90) | `RadiationOnlySolves` runs over every difference frequency with `qtf_direct = 1`, including those on the main grid, since the direct solve needs the matrix at Ω after the frequency loop. |

**B2 by the Stokes form with G.** B2 = −n·[(X·∇)u] has the same sharp-edge problem as in the potential part, so the Stokes form of [B2](#qtf-potential-part) is used with G(x_i, ·) in place of ψ. G(x_i, ·) is singular at x_i on S: ∇_ξ G goes like 1/ρ², with a leading part odd in the direction, so the surface term is a principal value and the identity keeps its form [derived, QTFDirect.f90 header]. One point per panel for the density on the 1/ρ² kernel converged only as h^0.65 at collocation points on S (0.186, 0.116, 0.075 at 128, 512 and 2048 panels on a sphere; the same rule gives O(h²) off the surface). So on the panels that `CALGREEN` takes as near (within 50 panel sizes) the Rankine part 1/r + 1/r₁ is integrated against a density linear on the panel, with closed-form first moments over flat panels (`rankine_moments`); the density gradient comes from `tangential_fit`. The wave part and the far panels keep the one-point rule. Panels that share a node with a panel across a sharp edge (normals more than `qtf_edge_angle` apart) keep the one-point rule, because u is singular there and the fitted gradient means nothing. The mask is built once per run (`direct_setup`).

The panel integrals of G and ∇_ξ G come from the stored arrays of the left side at Ω (`CGRN`, `RKBN`, `DS`, and the multi-body versions), so `CALGREEN` must have run at Ω before the right side is formed. They are single precision (`GRN_KIND = 4`); `green_panel_integrals` computes them afresh and is used by the tests.

**Solve.** At each Ω the matrix is the radiation matrix of `ASSB_LEFT` (LU factors with `ISOLV = 1`). With `ISOLV = 2` it is left unfactored for GMRES, and `direct_factor` factors it once per Ω (the direct solve is its last use at that Ω; a second factorisation would give a wrong φ_D with no error, so `direct_factor` records the Ω and stops on a repeat). The right sides of all pairs of one Ω are one block, `ZGETRS` with `NRHS` = the number of pairs, per symmetry block, in blocks of at most `DIRECT_BLOCK_BYTES` (256 MB) of work arrays. A pair with Ω < 0 is treated as `potential_pair` treats it: φ_D solves the problem at −|Ω|, whose matrix and Green function are the conjugates of those at |Ω|. The equal-frequency pairs (Ω = 0) need no solve: `.12d` equals `.10d` there, as `.11d` does.

**Irregular frequency removal (`If_remove_irr_freq = 1`).** The lid rows add r_l = ∬_S G(x_l, ξ) b dS at the lid points, with the stored lid integrals `DGRN` and `PKBN`, and the B2 part by the same Stokes form with G(x_l, ·), which is smooth since x_l is off the hull. The system is solved in the least-squares form of HAMS, (AᵀA) x = Aᵀ r with `CMAT`, and Aᵀ r formed as `ASSB_RBC_IRR` forms `DRMAT`. A mesh with no lid points stops the run (`check_lid`).

**Memory.** Besides the arrays of the potential part, `qtf_direct = 1` holds T_dir, `tdir(6·nbody, nfreq, nfreq, heading pairs)`, complex, 16 · 6 · nbody · nfreq² · (heading pairs) bytes, the size of T_pot. With `write_terms` 1 or 2 it also holds the arrays `tdir_i2` and `tdir_d`, of the same size each, for `QTF_DirectTerms.txt`. So the direct method adds 1 such array, or 3 with `write_terms`. The log line `QTF: with psi and the potential part the store takes ... MB` counts T_dir (it doubles the heading pairs, see `SecondOrderDriver.f90`) but not `tdir_i2` and `tdir_d`. For 101 frequencies, one body and one heading each array is 1.0 MB, so this matters only for many frequencies and heading pairs (for 1000 frequencies and 4 headings with `mean_drift = 2`, 16 pairs: 1.5 GB each [derived]). The work arrays of `direct_at_omega` are about 25 complex values per panel and pair (with the lid rows), in blocks of at most 256 MB.

**Notes for development.**

- The B2 test of the production path (`rhs_b2_stored_matches_direct`) compares the stored integrals with fresh ones: they agree to the single-precision storage (1.6e-8 relative).
- Finite depth h = 50 m on the 2048-panel sphere gives B2 Stokes against direct within 1.2e-3 with observed order 1.7 over 128, 512 and 2048 panels (`rhs_b2_stokes_vs_direct_smooth`). In deep water the Green function and its gradient are inconsistent at 0.46%, which is in the code before this phase; the test prints it and does not gate on it.
- The first-order code builds the lid point differently in two places: `CALGREEN_IRR` uses the mirrored lid point for the near and far choice, `ASSB_LEFT_IRR` and `ASSB_RBC_IRR` the unmirrored one (`CalGreenFunc.f90`:236 against `AssbMatx_irr.f90`:123 and :259). The direct solve follows the assembly. It was not changed in 7a; it is to be revisited if a second-order result with `If_remove_irr_freq = 1` looks wrong.
- After a change to the layout of a type in `SecondOrderTypes`, delete the stale objects of fpm's build before testing, as for the potential part.

**Tests.** Unit tests in `test/test_qtf_direct.f90` (collected as `qtf_direct`; right side, B2 forms, solve, lid), `test/test_qtf.f90` (`qtf_direct_omegas_all`, `qtf_direct_missing`) and `test/test_second_order_input.f90` (`qtf_direct_key`, `qtf_direct_needs_body`). Integration tests in `test/test_integration_drift.f90`: `run_qtf_direct_cases` and `run_qtf_direct_sharp_m2` make the runs; `qtf_direct_file_layout`, `direct_vs_indirect_round`, `direct_vs_indirect_sharp`, `direct_half_matches_full`, `direct_two_bodies`, `direct_diagonal`, `direct_gmres` (`ISOLV = 2`), `direct_irr_runs`, `direct_vs_indirect_irr_file` and `direct_no_lid_stops` check the file, the agreement with `.11d`, the symmetry images, two bodies, the diagonal, the LU-per-Ω path, the lid and the stops. The limits are measured values with a margin; the printed values are in the test output. With the direct method the unit tester has 271 tests and the integration tester 196, all passing [measured on this branch at 493674a plus these docs, ifort 2022.1, Windows, 2026-11-12].

Sources: Molin (1979), Ogilvie (1983) and Chen (1988) as under the potential part. Kurnia, R., Ducrozet, G. (2023). NEMOH: open-source boundary element solver for first- and second-order hydrodynamic loads. Comput. Phys. Commun. 292, 108885 (the direct and the indirect load). The equality of the two loads is derived above by Green's second identity.

### Partial second-order field (QTFField.f90)

`QTFField.f90` writes `QTF_FieldPartial.txt` (`qtf_field = 1`): the second-order pressure and elevation at the field points of `ControlFile.in`, from the quadratic part and the second-order incident potential of `SecondOrderIncident.f90` only. The second-order diffracted potential and, for a free body, the potential radiated by its second-order motion are omitted. Near the bodies, where the free-surface forcing drives the second-order scattered wave, the omitted part is not small ([Partial second-order field](input_files.md#partial-second-order-field-at-the-field-points)). The complete field is a later phase.

`field_pair` uses Bernoulli's equation expanded about z = 0 (Dean and Dalrymple 1991), time factor $e^{-i\omega t}$: $p^{(2)} = -\rho[\phi^{(2)}_t + \frac{1}{2}|\nabla\phi^{(1)}|^2]$ at a fixed point below z = 0, and $\eta^{(2)} = -\frac{1}{g}[\phi^{(2)}_t + \frac{1}{2}|\nabla\phi^{(1)}|^2 + \eta^{(1)}\phi^{(1)}_{zt}]$ on z = 0, with $\eta^{(1)} = (i\omega/g)\phi^{(1)}$ and $\phi^{(1)}_{zt} = -i\omega\,\phi^{(1)}_z$. For the ordered pair (i, j) of unit waves and the kernel P of `NearFieldDrift` (`pair_minus` for the difference, `pair_plus` for the sum), with $u_m$ the components of $\nabla\phi^{(1)}$:
- quad_p = −(ρ/2) Σ_m P(u_im, u_im, u_jm, u_jm);
- quad_η = −(1/g) [½ Σ_m P(u_im, u_im, u_jm, u_jm) + P(η_i, −iω_i u_i3, η_j, −iω_j u_j3)];
- pot_p = iρΩ φ₂(x), pot_η = (iΩ/g) φ₂(x), with Ω = ω_i − ω_j (difference) or ω_i + ω_j (sum) and φ₂ the second-order incident potential of the pair (`incident2_phi`), made with the same kernel.

These are T⁻_ij and T⁺_ij in the sense of `QTF_Quadratic.10d` and `.10s`: the physical $e^{-i\Omega t}$ amplitude is, for i ≠ j, 2 a_i a_j* (quad + pot) for the difference and 2 a_i a_j (quad + pot) for the sum; for i = j it is a_i² (quad + pot) for the sum and |a_i|² (quad + pot) for the mean. For the sum kernel, `incident2_pair` (optional argument `kernel`, default the difference) takes Ω = ω_1 + ω_2, K = k_1 e_1 + k_2 e_2 and the free-surface forcing with `pair_plus`, with no diagonal short-cut. The incident potential takes the gravity of `wave_gravity`, as the QTFs do; Bernoulli takes `G` of `Const_mod`, the one of the first-order elevation. The pressure and the elevation are conjugated and scaled only in `write_field_row` (`SecondOrderOutput.f90`).

`build_field_points` turns the points into a `cs_points` set of kind `CS_FREE`, so that `fields_green` (`FieldVelocity.f90`) computes φ at every point. Its rules: |z| ≤ 1e-6 is an elevation point (z set to 0), z < −1e-6 a pressure point, z > 1e-6 rejected; a point within `psize(i)` of a hull centroid i is on the hull; a point is inside body b when the winding number of body b (`body_winding`) is more than 0.5. It sums the solid angles (Van Oosterom and Strackee 1983, IEEE Trans. Biomed. Eng. 30(2), 125–126) of the triangles of the hull panels, symmetry images included, and of their mirror in z = 0, each turned to face out of the body; the hull and its mirror form a closed surface for a surface-piercing body, and the mirror of a submerged body lies above z = 0. An earlier rule, the side of the nearest centroid's panel plane, stopped valid points below the hull near an edge where fine wall panels meet a coarse bottom (unit test `points_mixed_panel_sizes`). The points are checked in order, and the first that fails gives the message. `psize` and `body` of a point are those of the nearest hull panel; the finite-difference step is `FIELD_FD_STEP` = 0.02 of that size.

At the end of each frequency of the QTF grid, `field_frequency` (`SecondOrderDriver.f90`) calls `fields_green` for all headings and `qfield_keep` stores φ and ∇φ: 4 complex values per point, frequency and heading. `fields_green` builds one influence row per stencil point and applies it to every heading: 7 rows for a point below z = −δ (the centre and ±δ along x, y, z), 5 for a point from z = −δ up (the centre and ±δ along x, y, with ∂φ_P/∂z = ν φ_P). After the frequency loop, `qfield_finish` forms the pairs in the loops of `qtf_finish` (kernel, β1, β2, i, j, then the point), with `qtf_max_dw` on the difference pairs, so that the pairs are those of `.10d` and `.10s`. Its state is module state, as for `QTFDriver.f90`.

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



