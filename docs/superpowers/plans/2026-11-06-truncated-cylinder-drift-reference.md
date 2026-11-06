# Truncated Cylinder Drift Reference Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A semi-analytical (eigenfunction matching) reference for the first-order coefficients and mean drift (F_x by two
routes, F_z, M_y, terms I to V) of a fixed and a free truncated cylinder with a sharp bilge, and a comparison against the
existing HAMS-MREL Snellius results.

**Architecture:** One Python script. Each azimuthal order is solved by matching an exterior and an interior expansion at r = a,
with the gap velocity in a corner-singular Gegenbauer basis. Drift comes from Maruo's far-field formula and from the near-field
terms integrated with corner-weighted quadrature. Self-checks are registered functions run by a `check` subcommand, which exits
with 1 on a failure.

**Tech Stack:** Python 3, numpy, scipy (special, integrate, linalg). No other packages.

**Spec:** `docs/superpowers/specs/2026-11-06-truncated-cylinder-drift-reference-design.md`

## Global Constraints

- numpy and scipy only. Run as `python -I truncated_cylinder_reference.py ...` from `pre_postprocessing_python/SecondOrder/`;
  the script puts its own folder on `sys.path` itself to import `make_mesh` and `drift_convergence`.
- Linear potential flow, time factor e^(-iωt), incident wave of unit amplitude along +x, heading 0 only.
- Defaults: a = 10 m, d = 10 m, h = 100 m; ρ and g are `make_mesh.RHO` and `make_mesh.G`.
- Free body: `make_mesh.hydrostatics` rule applied to the exact circle (uniform solid, z_G = z_B = -d/2). Rotation centre and
  `XR` at the origin.
- Output normalisation is HAMS's WAMIT format, derived from the HAMS source (`src/PrintOutput.f90`, `src/PotentWavForce.f90`,
  `src/SecondOrderOutput.f90`) and stated in the module docstring with the file and routine each rule comes from.
- Near-field terms are written from spec section 5.3 of `2026-11-04-second-order-waves-design.md`, not from the Fortran.
- Frequencies: `drift_convergence.WIDE_FREQUENCIES` (0.40 ... 1.80 rad/s).
- Tolerances: bottom-mounted limit 10⁻³ (gap 10⁻³h); Haskind 10⁻⁶; A₁₅ = A₅₁ and B₁₅ = B₅₁ 10⁻⁸; energy 10⁻⁶; near-field F_x
  against far-field F_x 5 × 10⁻⁴ relative, fixed and free, every frequency; reference point to round-off (10⁻¹⁰ relative).
- No Nemoh or WAMIT source, files or numbers in the repository. No Fortran changes. No CI test. No HAMS runs.
- Commit trailer: `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Review Focus

1. **Sign and conjugation against HAMS.** A flipped sign or a conjugate in the conversion looks plausible. Expected: at the
   finest HAMS mesh with `source`, F_x, F_z and M_y have the reference's sign wherever |reference| is above the floor. Test in
   Task 7.
2. **Overflow at small gaps and high orders.** The bottom-mounted limit (gap 10⁻³h) and orders up to 2M need scaled Bessel
   functions (`ive`, `kve`, ratios). Expected: no inf or nan anywhere. Test in Task 2 (bottom limit) and Task 6 (M and 2M).
3. **HAMS rows near its irregular frequencies** (1.548 and 1.60 rad/s for this hull). Expected: the comparison marks rows within
   5% of `drift_convergence.IRREGULAR_FREQUENCIES` instead of reporting them as HAMS errors. Test in Task 7.
4. **Values near zero** (free surge F_x at 0.40 rad/s, M_y near its sign change). Expected: percentages use the floor of
   `drift_convergence.percent` and print an empty field. Test in Task 7.
5. **Result folders with missing files** (no `MeanDriftControlSurface.7` or `MeanDriftTerms.txt`). Expected: the comparison
   prints "missing" for those columns and continues. Test in Task 7.

---

### Task 1: Modes, Gegenbauer basis, inner products

**Files:**
- Create: `pre_postprocessing_python/SecondOrder/truncated_cylinder_reference.py`

**Interfaces:**
- Produces:
  - `Cylinder(radius=10.0, draft=10.0, depth=100.0)`: dataclass, with `gap = depth - draft`.
  - `wave_numbers(omega, depth, count) -> (k, kn)`: k real, kn the first `count` roots of kn tan(kn h) = -ν, increasing.
  - `gegenbauer_basis(p_max) -> Basis`: evaluates u_p(t) = (1 - t²)^(-1/3) C_2p^(1/6)(t) on t ∈ [-1, 1], with Gauss-Jacobi
    nodes and weights for α = β = -1/3.
  - `inner_exterior(cyl, k, kn, basis) -> ndarray`: ∫ over the gap of u_p times each exterior vertical mode.
  - `inner_interior(cyl, n_max, basis) -> ndarray`: ∫ over the gap of u_p times each interior mode cos(λ_n (z + h)).
  - `CHECKS: dict[str, Callable[[], tuple[bool, str]]]` and `main(argv)` with subcommand `check [NAME ...]`.

- [ ] **Step 1: Write the checks `roots` and `inner_products`.** `roots`: for ω in `WIDE_FREQUENCIES`, |k tanh(kh) - ν| and
  |kn tan(kn h) + ν| below 10⁻¹² ν, and kn in ((n - ½)π/h, nπ/h). `inner_products`: the matrices of both functions equal
  `scipy.integrate.quad` of the same integrands (weight handled with `weight="alg"`) to 10⁻¹⁰ for p ≤ 10, n ≤ 20.
- [ ] **Step 2: Run** `python -I truncated_cylinder_reference.py check roots inner_products`. Expected: FAIL (functions missing).
- [ ] **Step 3: Implement the module docstring (conventions, sources), `Cylinder`, `wave_numbers`, `gegenbauer_basis`,
  `inner_exterior`, `inner_interior`, `CHECKS`, `main`.** Inner products by Gauss-Jacobi on the reflected interval; closed form
  where it exists.
- [ ] **Step 4: Run the same command.** Expected: both PASS, exit 0.
- [ ] **Step 5: Commit** `feat: add modes and corner basis for the truncated cylinder reference`.

### Task 2: Diffraction, excitation and far-field amplitude

**Interfaces:**
- Consumes: Task 1.
- Produces:
  - `Solution`: per order m, exterior coefficients (propagating and evanescent), interior coefficients, gap coefficients;
    `potential(r, theta, z)` and `gradient(r, theta, z) -> (3,)` complex, anywhere in the fluid.
  - `solve_diffraction(cyl, omega, orders, n_modes, p_max) -> Solution`.
  - `excitation(cyl, sol) -> dict[int, complex]`: X for modes 1, 3, 5, in HAMS's normalisation.
  - `far_field_amplitude(sol) -> ndarray`: a_m per order, so that η_scattered ≈ Σ_m a_m e^(imθ) e^(ikr)/√r.

- [ ] **Step 1: Write the checks `bottom_limit_x1`, `matching`, `scattering_energy`.**
  - `bottom_limit_x1`: gap 10⁻³h, ω in `WIDE_FREQUENCIES`: |X₁| equals the MacCamy-Fuchs value 4 tanh(kh)/(k² |H₁'(ka)|) in
    units of ρgA, to 10⁻³.
  - `matching`: at ω = 0.95, orders 0, 1, 2, 5: max pointwise |Δφ| and |Δ∂φ/∂r| at r = a on the gap, at depths 10%, 50%, 90%
    of the gap from the corner, falls monotonically for p_max = 10, 20, 40.
  - `scattering_energy`: fixed body, ω in `WIDE_FREQUENCIES`: scattered energy flux equals the optical-theorem term to 10⁻⁶.
- [ ] **Step 2: Run** `check bottom_limit_x1 matching scattering_energy`. Expected: FAIL.
- [ ] **Step 3: Implement `Solution`, `solve_diffraction`, `excitation`, `far_field_amplitude`.** One complex system of size
  p_max + 1 per order. Scaled Bessel functions (`ive`, `kve`, `h1vp`/`hankel1` ratios) everywhere.
- [ ] **Step 4: Run the same command.** Expected: PASS.
- [ ] **Step 5: Commit** `feat: add diffraction solution and excitation`.

### Task 3: Radiation, added mass and damping

**Interfaces:**
- Consumes: Task 2 `Solution`.
- Produces:
  - `solve_radiation(cyl, omega, mode, n_modes, p_max) -> Solution` for mode 1, 3, 5 (unit velocity), with the heave and pitch
    particular solutions in the interior.
  - `added_mass_damping(cyl, omega, rads) -> (A, B)`: 6 × 6 real arrays, modes 1, 3, 5 and their couplings filled, in HAMS's
    normalisation.

- [ ] **Step 1: Write the checks `reciprocity` and `radiation_energy`.** `reciprocity`: A₁₅ = A₅₁ and B₁₅ = B₅₁ to 10⁻⁸; X from
  diffraction equals the Haskind value from the radiation solutions to 10⁻⁶, modes 1, 3, 5. `radiation_energy`: B₁₁, B₃₃, B₅₅
  equal the radiated energy flux at infinity to 10⁻⁶. All at every ω in `WIDE_FREQUENCIES`.
- [ ] **Step 2: Run** `check reciprocity radiation_energy`. Expected: FAIL.
- [ ] **Step 3: Implement `solve_radiation` and `added_mass_damping`.**
- [ ] **Step 4: Run the same command.** Expected: PASS.
- [ ] **Step 5: Commit** `feat: add radiation solutions, added mass and damping`.

### Task 4: Free body, motions and far-field drift

**Interfaces:**
- Consumes: Tasks 2 and 3.
- Produces:
  - `body_matrices(cyl) -> (M, C)`: 6 × 6, the `make_mesh.hydrostatics` rule for the exact circle.
  - `motions(cyl, omega, diff, rads) -> dict[int, complex]`: ξ for modes 1, 3, 5, in HAMS's normalisation; zero for `fixed`.
  - `total_solution(diff, rads, xi) -> Solution`: φ_total = φ_I + φ_D + Σ ξ_j φ_j (with the -iω factor of the convention).
  - `far_field_fx(cyl, omega, total) -> float`: Maruo, over ρgA².

- [ ] **Step 1: Write the checks `bottom_limit_drift` and `body_matrices`.** `bottom_limit_drift`: gap 10⁻³h, fixed:
  `far_field_fx` equals `drift_convergence.maccamy_fuchs_fx` to 10⁻³ at every ω. `body_matrices`: `body_matrices(cyl)` equals
  `make_mesh.hydrostatics` on a 4096-sided polygon of the same radius to 10⁻⁵ in every entry.
- [ ] **Step 2: Run** `check bottom_limit_drift body_matrices`. Expected: FAIL.
- [ ] **Step 3: Implement `body_matrices`, `motions`, `total_solution`, `far_field_fx`.**
- [ ] **Step 4: Run the same command.** Expected: PASS.
- [ ] **Step 5: Commit** `feat: add free body motions and far-field drift`.

### Task 5: Near-field drift with corner quadrature

**Interfaces:**
- Consumes: Task 4.
- Produces:
  - `body_quadrature(cyl, n_corner, n_far) -> Quadrature`: points, weights and normals on the wall and the bottom, Gauss-Jacobi
    with weight ρ^(-2/3) in the distance from the corner within one gap height of it (or a/2 on the bottom, whichever is
    smaller), Gauss-Legendre beyond.
  - `near_field_terms(cyl, omega_i, sol_i, xi_i, omega_j, sol_j, xi_j, quad, xr=(0, 0, 0)) -> dict[str, ndarray]`: keys
    "I", "Ia", "II", "III", "IV", "V", each a (6,) complex array of T⁻ᵢⱼ (spec 5.4), components 1, 3, 5 filled. For i = j the
    real part is the mean drift over ρgA².
  - `near_field(cyl, omega, total, xi, quad, xr) -> (terms, total)`: the diagonal, real.

- [ ] **Step 1: Write the checks `fx_identity`, `reference_point`, `pair_symmetry`.**
  - `fx_identity`: fixed and free, every ω in `WIDE_FREQUENCIES`: |near F_x - far F_x| ≤ 5 × 10⁻⁴ |far F_x|.
  - `reference_point`: xr = (0, 0, -3) against (0, 0, 0): ΔM_y = (Δxr × F)_y to 10⁻¹⁰ relative, at 0.70 and 1.40 rad/s.
  - `pair_symmetry`: T⁻ᵢⱼ = conj(T⁻ⱼᵢ) for ωᵢ = 0.70, ωⱼ = 0.95, to 10⁻¹².
- [ ] **Step 2: Run** `check fx_identity reference_point pair_symmetry`. Expected: FAIL.
- [ ] **Step 3: Implement `body_quadrature`, `near_field_terms`, `near_field`.** Term V and the waterline unit normal as in
  spec 5.3. If `fx_identity` cannot reach 5 × 10⁻⁴ with n_modes ≤ 400 and p_max ≤ 60, stop and report BLOCKED with the
  measured errors: the spec's fallback (approach B) is a separate design.
- [ ] **Step 4: Run the same command.** Expected: PASS.
- [ ] **Step 5: Commit** `feat: add near-field mean drift with corner quadrature`.

### Task 6: Reference values and convergence

**Interfaces:**
- Consumes: Tasks 1-5.
- Produces:
  - `reference(cyl, frequencies, motion, n_modes, p_max, orders) -> list[dict]`: one row per ω with A, B (11, 13, 15, 33, 35,
    55), |X| (1, 3, 5), |RAO| (1, 3, 5, free only), F_x far, F_x near, F_z, M_y, terms I-V of F_x, F_z, M_y.
  - Subcommand `reference [--out CSV]`: writes fixed and free rows at the default resolution, then the convergence table, then
    runs all checks. Exit 1 if a check fails.
  - `DEFAULT_RESOLUTION = dict(n_modes=..., p_max=..., orders=...)`: chosen in Step 3.

- [ ] **Step 1: Write the check `convergence`.** For F_x far, F_z, M_y and B₁₁ at 0.70, 0.95 and 1.40 rad/s, fixed and free:
  values at three (n_modes, p_max) and at orders M and 2M are finite, and the spread between the two finest is below 10⁻⁴
  relative. Prints the table.
- [ ] **Step 2: Run** `check convergence`. Expected: FAIL.
- [ ] **Step 3: Implement `reference` and the subcommand; choose `DEFAULT_RESOLUTION` as the coarsest that meets the
  `convergence` spread; record the run time in the docstring.**
- [ ] **Step 4: Run** `python -I truncated_cylinder_reference.py reference --out tc_reference.csv` (in a scratch folder).
  Expected: all checks PASS, CSV written, run time printed.
- [ ] **Step 5: Commit** `feat: add reference values and convergence table`.

### Task 7: Comparison with HAMS results

**Interfaces:**
- Consumes: Task 6 `reference`; `drift_convergence.read_drift`, `read_terms`, `read_hams_first_order`, `percent`,
  `IRREGULAR_FREQUENCIES`.
- Produces: subcommand `compare RESULTS [--out CSV]`. For every folder under RESULTS whose name starts with `a_sharp_`, reads
  `output/Wamit_format/` and prints per case and ω: % difference of A, B, |X|, |RAO|, F_x near, F_x far, F_z, M_y, each term,
  and control-surface F_x where `MeanDriftControlSurface.7` exists. The case's panel count and the polygon area error
  (32 × multiple sides) head each block.

- [ ] **Step 1: Write the checks `compare_fixture` and `compare_signs`.**
  - `compare_fixture`: builds a temporary RESULTS with one case folder whose files are written from the reference itself
    (in HAMS's format), plus one folder missing `MeanDriftTerms.txt` and `MeanDriftControlSurface.7`, plus one row at
    1.60 rad/s. Expected: 0% differences (to 10⁻⁹) for the first, "missing" columns for the second, the 1.60 row marked
    "irregular", and empty fields where |reference| is below the floor.
  - `compare_signs`: if `../../../validation-runs/snellius-drift/results/a_sharp_free_src_m4` exists, F_x, F_z and M_y have the
    reference's sign wherever |reference| is above the floor; otherwise PASS with a "skipped" message.
- [ ] **Step 2: Run** `check compare_fixture compare_signs`. Expected: FAIL.
- [ ] **Step 3: Implement `compare`.**
- [ ] **Step 4: Run the checks, then**
  `python -I truncated_cylinder_reference.py compare ../../../validation-runs/snellius-drift/results --out tc_compare.csv`
  (in a scratch folder). Expected: checks PASS; the table covers all 14 `a_sharp_*` cases.
- [ ] **Step 5: Commit** `feat: compare HAMS results with the truncated cylinder reference`.

### Task 8: Documentation

**Files:**
- Modify: `docs/mean_drift.md` (new section "Semi-analytical reference: sharp truncated cylinder", after "Comparison with
  Nemoh"), `docs/python_scripts.md` (one entry).

- [ ] **Step 1: Verify the sources.** Confirm each citation in spec section 8 (authors, year, title, journal, volume, pages)
  against the publisher's page; correct the spec where they differ.
- [ ] **Step 2: Write the section:** method in four sentences, the checks with their measured values, reference values at
  0.70, 0.95, 1.40 rad/s fixed and free, and a table of HAMS against the reference at the finest mesh of each velocity method
  (F_x near, F_x far, F_z, M_y), from Task 7's CSV. Measured values only, marked as such.
- [ ] **Step 3: Run** `python -I truncated_cylinder_reference.py reference` once more. Expected: all checks PASS.
- [ ] **Step 4: Commit** `docs: add the truncated cylinder drift reference and comparison`.
