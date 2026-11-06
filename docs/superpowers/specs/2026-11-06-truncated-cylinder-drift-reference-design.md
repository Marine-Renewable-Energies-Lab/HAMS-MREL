# Semi-analytical mean drift reference for a truncated cylinder

Design, 2026-11-06. Branch `validation/truncated-cylinder-reference`, cut from `second-order-drift` at ea8cff1.

## 1. Purpose

HAMS-MREL's near-field mean drift converges slowly on hulls with a sharp edge. On the sharp free truncated cylinder at 7168
panels, `differencing` is -136% and `source` +10.3% from far-field at 0.70 rad/s ([mean_drift.md](../../mean_drift.md)). Nemoh is a
second panel code with the same error sources, and it failed at 4320 panels. No WAMIT drift values are available.

This tool gives a reference for that hull that is exact up to series truncation, with a stated error bar:

- first order: A, B, X and RAOs in surge, heave and pitch;
- mean drift: F_x by far-field and by near-field, and F_z and M_y by near-field, each near-field term separately, fixed and free.

It checks HAMS near-field (`differencing` and `source`), far-field, and the control-surface method of
`feature/control-surface-drift`, against the existing Snellius results. It needs no new HAMS runs.

Out of scope: the rounded bilge, several bodies, headings other than 0 (the body is axisymmetric), and CI.

## 2. Problem

Truncated vertical cylinder, radius a, draft d, in water of depth h, sharp bilge. Defaults are the study's sharp hull: a = 10 m,
d = 10 m, h = 100 m (`DRIFT_CYL_MESH` and `HULLS["sharp"]` in `drift_convergence.py`). ρ and g are HAMS's values. Linear
potential flow, time factor e^(-iωt), plane incident wave of unit amplitude along +x. Mass, centre of gravity, inertia and
hydrostatic stiffness of the free body follow `make_mesh.hydrostatics`, the rule that wrote HAMS's `Hydrostatic.in` (a uniform
solid at its draft, z_G = z_B = -d/2), applied to the exact circle. HAMS's hull is a polygon of 32 × multiple sides, so its
mass differs by the polygon's area error, which the comparison prints per mesh; that is part of HAMS's discretization error. The study's result
folders keep only the output, so the input files are not read. Rotation centre and `XR` are at the origin, at the still waterline, as in the study.

The conversion to HAMS's output normalisation (WAMIT format, ρgA, ρgA²L^k, phases) is taken from the HAMS source and stated in
the script's docstring, as `annular_lid_reference.py` on `feature/damping-lid` does.

## 3. First-order solution

Each azimuthal order m is solved on its own: m = 0 ... M for diffraction, m = 1 for surge and pitch, m = 0 for heave.

- **Exterior**, r > a, -h < z < 0: the propagating mode Z_0 with J_m(kr) (incident) and H_m^(1)(kr) (scattered), and evanescent
  modes Z_n with K_m(k_n r), n = 1 ... N.
- **Interior**, r < a, -h < z < -d: modes cos(λ_n (z + h)), λ_n = nπ/(h - d), with I_m(λ_n r), and r^m for n = 0. For heave and
  pitch, the standard particular solution meets the condition on the cylinder bottom.
- **Matching at r = a.** On -d < z < 0 the radial velocity is the body's. On the gap -h < z < -d the radial velocity u(z) is the
  unknown. The fluid corner at (a, -d) is a 270° wedge, so u behaves like (−d − z)^(-1/3). After reflection about the sea bed
  (∂φ/∂z = 0 there), u is expanded in even Gegenbauer polynomials C_2p^(1/6)(t) with weight (1 - t²)^(-1/3),
  t = (z + h)/(h - d), p = 0 ... P. Continuity of φ on the gap, projected on the same functions, gives a P + 1 system per m and
  per problem. Inner products are closed form or one-dimensional Gauss-Jacobi integrals.

Outputs: A and B (11, 13, 15, 33, 35, 55 and the symmetric partners), X (1, 3, 5), and for the free body the RAOs from the
equation of motion with the mass and stiffness above. The potentials are series that can be evaluated anywhere in the fluid,
including on the body and at the waterline.

## 4. Mean drift

**Far-field F_x.** Maruo's formula with the outgoing amplitude a(θ) = Σ_m a_m e^(imθ) from the exterior coefficients. For an
axisymmetric body it reduces to sums of a_m a*_(m+1). Same normalisation as `MeanDriftMomentum.8`.

**Near-field F_x, F_z and M_y.** Terms I to V of spec section 5.3 of
[2026-11-04-second-order-waves-design.md](2026-11-04-second-order-waves-design.md) (WAMIT eqs. 15.54-15.57): waterline,
velocity squared, displacement times gradient, rotation, hydrostatic. Written from those equations, not from the Fortran. M_y
about `XR`. Each term is reported separately, in the order and meaning of `MeanDriftTerms.txt`, and the totals in the
normalisation of `MeanDriftPressure.9`. F_y, M_x and M_z are zero by symmetry and are not computed.

**Corner.** |∇φ|² grows like ρ^(-2/3) at the corner. On the wall (exterior series) and the bottom (interior series), each face
is integrated with Gauss-Jacobi quadrature in the distance from the corner, weight ρ^(-2/3), and Gauss-Legendre away from it.
The series are evaluated with up to a few hundred modes. Convergence is shown with a sequence of N and P, with Richardson
extrapolation where the sequence is regular.

**Free body.** Terms III, IV and the X₃ parts of term I use the reference's own RAOs. Term V uses the same formula as HAMS
(plan Task 7 of `2026-11-04-mean-wave-drift.md`).

**Frequency pairs.** Each quadratic term is a function of two first-order solutions (ωᵢ, ωⱼ), following the kernel T⁻ᵢⱼ, T⁺ᵢⱼ of
spec section 5.4. Only i = j is computed and checked in this work. The far-field route stays diagonal.

**Frequencies.** 0.40, 0.60, 0.70, 0.80, 0.95, 1.00, 1.10, 1.20, 1.40, 1.60, 1.80 rad/s (`WIDE_FREQUENCIES`). HAMS's
irregular frequencies do not exist in this solution.

## 5. Checks

The script exits with 1 if one fails. Each check states in the docstring what it covers and what it does not.

1. **Bottom-mounted limit.** Gap 10⁻³h: X₁ and far-field F_x agree with MacCamy-Fuchs (`maccamy_fuchs_fx` in
   `drift_convergence.py`) to 10⁻³.
2. **Reciprocity.** X from diffraction equals the Haskind value from radiation, to 10⁻⁶. A₁₅ = A₅₁ and B₁₅ = B₅₁ to 10⁻⁸.
3. **Energy.** B₁₁, B₃₃ and B₅₅ from body pressure equal the radiated energy flux at infinity, to 10⁻⁶. Fixed body: the
   scattered energy equals the optical-theorem term, to 10⁻⁶.
4. **Matching.** The pointwise mismatch of φ and ∂φ/∂r at r = a on the gap falls with P (P = 10, 20, 40), away from the corner.
5. **F_x identity.** Near-field F_x equals far-field F_x to 5 × 10⁻⁴ relative, fixed and free, at every frequency. This is the
   gate for F_z and M_y. If it is not met, the fallback is approach B: F_z and M_y by a middle-field formula on a control
   surface away from the corner, as a separate design.
6. **Reference point.** Moving `XR` by d changes M_y by (d × F)_y, to round-off.
7. **Convergence.** Every output is printed for three N and three P, and for M and 2M, with the spread. The spread is the error
   bar of the reference values, not a check of the formulation.

## 6. Comparison with HAMS

`python -I truncated_cylinder_reference.py compare RESULTS` reads the HAMS case folders of the study (for example
`validation-runs/snellius-drift/results`, outside the repository) for the sharp hull, fixed and free, all meshes (448 to 7168
panels), both velocity methods. For each case and frequency it prints the percentage difference in A, B, X, RAOs, F_x near
and far, F_z, M_y, and each near-field term. It reads `MeanDriftControlSurface.7` the same way where it exists. Percentages use
the same floor as `drift_convergence.py` where a value is near zero.

`python -I truncated_cylinder_reference.py reference` writes the reference values as CSV and runs the checks.

## 7. Files

- `pre_postprocessing_python/SecondOrder/truncated_cylinder_reference.py`: the solver, the checks, `reference` and `compare`.
  numpy and scipy only. Expected run time: seconds to a few minutes on the laptop.
- `docs/mean_drift.md`: a section with the reference values and the comparison table for the sharp hull. These are our own
  numbers.
- No Fortran changes. No CI test: a HAMS mesh small enough for CI cannot meet a tight tolerance against this reference.

## 8. Sources

Verified on 2026-11-06 (Task 8). Authors, years, titles, volumes and pages were correct; issue numbers and DOIs added,
Newman (1967) added, WAMIT section number added.

- Garrett, C.J.R. (1971). Wave forces on a circular dock. J. Fluid Mech. 46(1), 129-139. doi:10.1017/S0022112071000430
  (Cambridge Core article page).
- Yeung, R.W. (1981). Added mass and damping of a vertical cylinder in finite-depth waters. Appl. Ocean Res. 3(3), 119-133.
  doi:10.1016/0141-1187(81)90101-2 (Elsevier record via Scilit; ASCE reference list).
- Mavrakos, S.A. (1988). The vertical drift force and pitch moment on axisymmetric bodies in regular waves. Appl. Ocean Res.
  10(4), 207-218. doi:10.1016/S0141-1187(88)80005-1 (TRID, Scilit, NTUA repository).
- MacCamy, R.C. and Fuchs, R.A. (1954). Wave forces on piles: a diffraction theory. Technical Memorandum 69, Beach Erosion
  Board, December 1954 (ERDC library record, DTIC scan).
- Maruo, H. (1960). The drift of a body floating on waves. J. Ship Res. 4, 1-10. Not in Crossref; volume and pages from
  reference lists (Duan and Zhao 2011, an Osaka University dissertation). The issue is uncertain: one list gives No. 3, but
  Crossref has another paper at 4(3) pp. 1-16, and 4(4) has a gap at pp. 1-10; so the issue is left out.
- Newman, J.N. (1967). The drift force and moment on ships in waves. J. Ship Res. 11(1), 51-60. doi:10.5957/jsr.1967.11.1.51
  (Crossref).
- WAMIT User Manual, Version 7.5, section 15.8 "Mean drift forces based on pressure integration", eqs. (15.54) to
  (15.57) (checked in https://www.wamit.com/manual7.x/v75_manual.pdf: 15.54 force, 15.55 moment, 15.56 and 15.57 the
  hydrostatic parts).
