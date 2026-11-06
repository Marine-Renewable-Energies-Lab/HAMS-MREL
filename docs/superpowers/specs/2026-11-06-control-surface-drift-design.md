# Control-surface (middle-field) mean drift: design

Date: 2026-11-06. Branch base: `second-order-drift`. Status: design, no solver code yet.

This adds a third mean drift method next to near-field pressure integration (`NearFieldDrift.f90`) and far-field momentum (`FarFieldDrift.f90`). It is Phase 3 of `2026-11-04-second-order-waves-design.md`, written out in full.

Evidence labels: **[paper]** read in a published paper, with equation or page. **[doc]** read in the WAMIT v7.5 manual. **[measured]** from runs of this code. **[derived]** worked out here, not yet run.

## 1. Why

The Snellius convergence study of 2026-11-05 (`../validation-runs/snellius-drift/tables.md`, outside the repo) shows that near-field drift converges on rounded hulls but not on sharp edges [measured]:

| Case | Method | Panels | Error against far-field |
| :--- | :--- | ---: | ---: |
| Rounded cylinder, free, 0.95 rad/s, F_x | differencing | 7680 | +0.03% |
| Sharp cylinder, free, 0.95 rad/s, F_x | differencing | 7168 | +17.4% |
| Sharp cylinder, free, 0.95 rad/s, F_x | source | 7168 | +0.53% |
| Sharp cylinder, free, 0.70 rad/s, F_x | both | 7168 | -136% and +10.3% |
| Sharp box, fixed, heading 30, 1.0 rad/s, F_y | source | 7040 | +12.5% |
| Sharp box, free, heading 30, 1.0 rad/s, F_y | source | 1760 | -85% |

The velocity is singular at a sharp edge, and near-field integrates |∇φ|² on the hull. Chen (2007) puts the same point this way: near-field "suffers from poor convergence, as singularities are present in the velocity field" at corners [paper, p. 66].

Far-field is accurate but gives only F_x, F_y and M_z, and only the total over all bodies. The source velocity route cannot use the irregular-frequency lid (`IRSP = 1` is rejected). At the irregular frequencies 1.548 and 1.939 rad/s, on the free cylinders at 1792 and 1920 panels, it is off by 32 to 198% [measured, `tables.md`, irregular-frequency table].

A control-surface method evaluates the flow on a surface C at a distance from the body, where the field is smooth. Chen reports that it converges like far-field and, like near-field, gives loads per body [paper, abstract and §4.2].

## 2. Scope

In scope:
- Mean drift, unidirectional (`mean_drift = 1`) and all heading pairs (`mean_drift = 2`).
- Six components per body, for a single body or several bodies, in deep or finite depth.
- Two phases:
  - **Phase A:** F_x, F_y and M_z.
  - **Phase B:** F_z, M_x and M_y.

Out of scope:
- QTFs. The general formula (Chen eqs. 35 and 36) covers them [paper, §3.2], but they come later.
- WAMIT's second alternative (`IALTCSF = 2`), for C very close to the body [doc, ch. 11].
- Generalized modes. They are rejected, as for near-field.
- Exploiting symmetry on C: C is always a full surface.

## 3. Formulation

### 3.1 Source

The formulas are Chen, X.-B. (2007), *Middle-field formulation for the computation of wave-drift loads*, J. Eng. Math. 59:61-82, eqs. (37) to (41). WAMIT v7.5 eqs. 15.58 and 15.59 (Alternative 1, `IALTCSF = 1`) have the same terms [doc]. WAMIT is the reference for behaviour, not for code.

Notation, as in Chen:
- D is the fluid volume bounded by the hull H, the control surface C, and the part F of the mean free surface between the hull waterline Γ and the line Γc where C meets z = 0.
- n points into the fluid domain on every boundary [paper, App. C, eq. 61]. On C this is towards the body, as in WAMIT [doc, ch. 11].
- n̄ = (n̄₁, n̄₂, 0) is the unit horizontal normal on Γ and Γc, in the free-surface plane [paper, eq. 60].
- r = x − x₀, with x₀ space-fixed, and n₆ = (x−x₀)n₂ − (y−y₀)n₁.
- ⟨·⟩ is the mean over one period, Φ is the total first-order potential and η = −Φ_t/g [paper, eq. 10].

### 3.2 Phase A: horizontal force and yaw moment

Chen eqs. (39) to (41) [paper, p. 69]:

    F_x = < (ρg/2) ∮_Γc η² n̄₁ dl + ρ ∬_C [ Φ_n Φ_x − (∇Φ·∇Φ/2) n₁ ] ds >
    F_y = < (ρg/2) ∮_Γc η² n̄₂ dl + ρ ∬_C [ Φ_n Φ_y − (∇Φ·∇Φ/2) n₂ ] ds >
    M_z = < (ρg/2) ∮_Γc η² n̄₆ dl + ρ ∬_C { Φ_n [ (x−x₀)Φ_y − (y−y₀)Φ_x ] − (∇Φ·∇Φ/2) n₆ } ds >

What this means for the code:
- No hull integral appears. Chen shows that the hull terms have zero mean, because they group as ∂[∇Φ(X·n)]/∂t [paper, p. 68-69]. The body motions enter only through Φ, which includes the radiation potentials.
- For multiple bodies, C surrounds one body and gives the load on that body [paper, §3.2 and §4.3]. Φ is the total field of the whole system.
- x₀ is the body reference point `mesh%xr(:, b)`, the same point near-field uses. M_z is then about a space-fixed point, which Chen requires for the equivalence with far-field [paper, §5].

### 3.3 Phase B: vertical force and horizontal moments

Chen eqs. (37) and (38) [paper, p. 69], all components:

    F = < −(ρg/2) ∮_Γ Ξ k dl + ρ ∬_F (η Φ_zt + ∇Φ·∇Φ/2) k ds
          + (ρg/2) ∮_Γc η² n̄ dl + (ρ/2) ∬_C [ 2Φ_n ∇Φ − (∇Φ·∇Φ) n ] ds >
    M = < −(ρg/2) ∮_Γ Ξ (r∧k) dl + ρ ∬_F (η Φ_zt + ∇Φ·∇Φ/2)(r∧k) ds
          + (ρg/2) ∮_Γc η² (r∧n̄) dl + (ρ/2) ∬_C [ 2Φ_n (r∧∇Φ) − (∇Φ·∇Φ)(r∧n) ] ds >

with Ξ = (η² − 2ηζ₃) n₃/|cos γ| − 2η (X·n̄) on the hull waterline Γ [paper, p. 67]. Here X = T + R∧r is the body displacement, ζ₃ its vertical part, and γ the angle between the hull normal and the horizontal.

To these add the second-order hydrostatic loads (F_S, M_S) [paper, eq. 14]. They are non-zero only for F_z, M_x and M_y, and depend only on hull geometry and motions [paper, p. 65]. The Γ and F terms also contribute only to F_z, M_x and M_y [paper, p. 69], so Phase B does not change the Phase A rows.

**Open item before Phase B is built.** Chen takes (F_S, M_S) from Molin and Hairault (1983) [paper, ref. 6]. `NearFieldDrift` uses WAMIT eq. 15.56 and 15.57 for its hydrostatic term. Check whether these are the same quantity, with the same split between the waterline and hydrostatic terms, before reusing `T_HYDROSTATIC`. The gate in section 8.3 tests this.

### 3.4 Complex amplitudes and the mean

The time factor is e^{−iωt}, as in `first_order_state`. For two products of first-order quantities, the existing `pair_minus(f1, g1, f2, g2)` gives the mean for one state (f1 = f2, g1 = g2) and the pair term for two headings. Every product in sections 3.2 and 3.3 goes through it.

η = iωφ/g at z = 0, from η = −Φ_t/g. This sign is to be checked against `wave_gravity` and `incident_phi` in `FirstOrderState`, not assumed.

On z = 0, Φ_z = νΦ with ν = ω²/g, which is the linear free-surface condition [paper, eq. 3]. So the vertical velocity there is exact once φ is known, and η Φ_zt follows from φ alone [derived].

### 3.5 Quadrature

As in WAMIT's low-order method:
- On C and F: centroid value times panel area [doc, ch. 11].
- On Γc and Γ: segment midpoint value times segment length.

## 4. Field quantities on C, Γc and F

Each method needs the total φ and ∇φ at the points of section 3.5 for every heading. The incident part is analytic (`incident_phi` and the incident gradient in `FirstOrderState`). Only the perturbation part φ_P is computed. Two routes, chosen by `cs_velocity`.

### 4.1 `cs_velocity = green` (default)

**Potential.** φ_P(x) = (1/4π) Σ [G ∂φ_P/∂n − (∂G/∂n) φ_P] over the hull panels, as `field_influence_row` and `field_potential` in `WaterlineField.f90` compute it. With `IRSP = 1` the lid has already corrected the hull potentials, so the formula with hull panels only holds in the fluid. `WaterlineField` already relies on this, and it is tested by `drift_field_irsp` [measured, fixture].

**Velocity.** Central differences of φ_P with step δ:
- Below z = 0: 6 points, which is second order in δ.
- On z = 0 (Γc and F): the 4 horizontal points, with φ_z = νφ from section 3.4. No stencil point is above z = 0.

The analytic gradient of the Green formula would need ∇ of ∂G/∂n, which means second derivatives of G. `FINGREEN3D` and `INFGREEN3D` give only first derivatives, so differencing is the route that needs no new Green function code (the earlier spec, D2, said the same).

**Step size.** `cs_fd_step` is δ as a fraction of the local C panel size.
- Truncation error grows as δ², and round-off in the Green rows grows as ε_G/δ.
- The default is set by the sensitivity test in section 8.1, over δ from 0.005 to 0.1. The starting value is 0.02.

**One near-or-far decision per stencil.** `field_influence_row` takes the 1/r panel integral exactly when the point is within 50 panel sizes of a panel, and at the centroid otherwise. If two points of one stencil fell on either side of that switch, the difference would contain the jump between the two integrals, divided by δ [derived]. So all points of a stencil take the decision of its centre. This needs an optional argument on `field_influence_row` that gives the point used for the decision.

**Advantages.** It works with `IRSP = 1` and `ISOLV = 2`, and it uses the potential of the direct solve, the most accurate quantity in HAMS.

### 4.2 `cs_velocity = source`

This is Chen's route [paper, eq. 8]: φ_P = (1/4π) ∬ σ G, ∇φ_P = (1/4π) ∬ σ ∇G.
- σ comes from `velocity_source`, through a new optional output argument (σ is a local variable there today).
- The field gradient comes from `source_panel_integrals(mesh, xf, jel, flag, ...)`, which returns ∫G and ∫∇_p G over a panel at any point.
- φ on Γc and F comes from the same σ, so the route is consistent.

It keeps the existing rules: `ISOLV = 1` only, and `IRSP = 1` rejected. It costs one Green call per point per panel, against 4 to 6 for the green route.

### 4.3 When the fields are computed

The green route's rows for all points would need (stencil points × N_C) × N_hull × nsys × 2 complex numbers: about 5 GB at 7168 hull panels with 4000 C panels [derived]. So rows are not stored.
- The states of all headings of a frequency are kept, as `mean_drift = 2` already does. One state is a few complex arrays of length n.
- In `so_end_frequency`, a single pass over the points computes each point's rows and applies them to every heading at once.
- Memory per thread is O(nbase × nsys).
- The fields stored per heading, φ and ∇φ at the points, are small.

### 4.4 Dependency: smooth G in finite depth

Differencing turns any jump of G in the field point into a velocity error of order jump/δ. In finite depth, the `Limes` series acceleration scatters dG/dR by up to about 1e5 (scaled by h²) below R/h = 0.002 [measured, side-session note `limes-measurement-result.md`; the plan `2026-11-05-finite-depth-small-r-green.md` on that branch has the numbers].
- The fix (`LINTON` below R/h = 0.02 up to νh = 50) is on `fix/finite-depth-small-r-green`, merged into `second-order-drift-fixes`, not yet into `second-order-drift`.
- Bring it in before any finite-depth validation.
- Points of the bottom cap of C right under the hull have small R.

## 5. Surfaces

### 5.1 Control surface C

`cs_shape` sets how C is made: `cylinder`, `box`, or `file`.

**Generated (`cylinder`, `box`).**
- **Shape.** C is a closed vertical wall plus a flat bottom at depth d_C, one for each body, made in that body's frame (`LCS` for multi-body).
  - The cylinder is centred on the centre of the body's waterplane bounding box.
  - The box is aligned with the body axes.
  - Chen used a box and an ellipsoid [paper, §4.1-4.2]. WAMIT generates a cylinder or a box [doc, ch. 11].
- **Clearance.** `cs_margin` is the horizontal gap between the hull's largest horizontal extent and C. `cs_bottom_margin` is the gap below the deepest hull point. Both are in metres.
  - The default for both is 4 × the median hull panel size of that body. Chen gives a distance of 3 to 5 panel sizes [paper, p. 76].
- **Panel size.** `cs_panel_size` defaults to the median hull panel size, since Chen gives "a finesse equivalent to that of the body hull" [paper, p. 76].
  - The wall rows are graded towards z = 0, with the top row about half the nominal size. Chen notes that higher frequencies need a finer mesh near the free surface [paper, p. 76].
- **Reaching the seabed.** If d_C would reach the seabed (d_C ≥ h − cs_bottom_margin), the wall goes down to z = −h and there is no bottom panel. This is needed for bodies standing on the seabed.
  - On the seabed patch Φ_n = 0 and n = ±k, so it adds nothing to F_x, F_y or M_z [derived].
  - It does add ½|∇Φ|² n_z terms to the vertical components. Phase B integrates the seabed patch between the body and C for that case.
- **Normals.** They point towards the body (into D), as in section 3.1.

**User file (`file`).**
- C comes from `Input/CSMesh_{b}.pnl` in the `HullMesh` format, read with the existing panel reader.
- The normal must point towards the body. The run checks this with the sign of Σ n·(x_c − x_b) A, where x_b is the body centre, and stops with a message if it is wrong.

**Γc.** It is found by running `build_surface_geom` on C as a `panel_mesh`. Its waterline segments and `wl_nh` give the midpoints, lengths and n̄ on Γc. Whether `wl_nh` points into F when the normals point towards the body is to be confirmed by the incident-wave test (section 8.1), not assumed.

**Checks.**
- Generated C: every node of every other body must be outside C, by an exact point-in-cylinder or point-in-box test, and the body itself must be inside.
- File: no C centroid may be closer than one local panel size to any hull panel. Containment is the user's job, and this is documented. WAMIT has the same rule: a body's control surface "must not include or intersect any other body" [doc, ch. 11].
- In both cases, a C point inside a body stops the run.

**Output for inspection.** The generated C and F are written to `Output/` as `.pnl` files, as WAMIT writes `gdf_csf.dat` [doc, ch. 11].

### 5.2 Free-surface patch F (Phase B)

**Generated.** F is made of radial strips from each hull waterline node out to Γc, with `cs_fs_rows` rows graded towards the body.
- This works when the waterline is star-shaped about the body centre, meaning each ray from the centre crosses it once.
- F and C need not share nodes, since their integrals are separate.

**File.** F comes from `Input/FSMesh_{b}.pnl`, needed for waterlines that are not star-shaped (moonpools, several waterlines on one body). WAMIT's automatic surfaces also fail on concave waterlines and moonpools [doc, ch. 11].

**Seabed.** For a body standing on the seabed, the seabed annulus between the body and C is generated in the same way as F.

## 6. Code structure

New modules, on explicit arrays with an `NBODY`-independent interface (earlier spec, D7):

| Module | Does | Uses |
| :--- | :--- | :--- |
| `ControlSurfaceMesh.f90` | Generate or read C, F and the seabed patch per body; checks; Γc via `build_surface_geom`; write for inspection | `SecondOrderTypes`, `SurfaceGeometry`, `ReadPanelMesh` reader |
| `FieldVelocity.f90` | φ and ∇φ at points for all headings, by the green or source route | `WaterlineField` (`field_influence_row`, `field_potential`), `SurfaceVelocity` (`source_panel_integrals`, σ from `velocity_source`), `FirstOrderState` (incident φ and ∇φ) |
| `ControlSurfaceDrift.f90` | `control_surface_pair(cs, fields1, fields2, ...)`: the six components per body for one state pair, and each term apart | `NearFieldDrift` (`pair_minus`; Phase B: `waterline_point`, waterline φ, hydrostatic term if section 3.3 confirms it) |

Changes to existing code:
- `SecondOrderTypes`: options, a `control_surface` type and a field-values type.
- `SecondOrderInput`: keys and validation (section 7).
- `WaterlineField`: an optional point for the near-or-far decision in `field_influence_row` (section 4.1).
- `SurfaceVelocity`: an optional σ output from `velocity_source`.
- `SecondOrderDriver`: keep all heading states when `control` is on; call the control-surface step in `so_end_frequency`.
- `SecondOrderOutput`: `MeanDriftControlSurface.7`.
- Build lists: `CompileHAMS.bat` and `HAMS_MREL.vfproj` by hand (`fpm` picks up `src/`).

Phase B needs the hull waterline φ at Γ for Ξ. It reuses `waterline_extrapolate` or `waterline_field`, as chosen by `waterline_method`.

## 7. Input and output

`SecondOrder.in`. New keys are optional except where stated:

```text
drift_methods    = far control   ! far, near, control, any combination
cs_velocity      = green         ! green, source
cs_shape         = cylinder      ! cylinder, box, file
cs_margin        = 0             ! m; 0 means 4 x median hull panel size
cs_bottom_margin = 0             ! m; 0 means 4 x median hull panel size
cs_panel_size    = 0             ! m; 0 means median hull panel size
cs_fd_step       = 0.02          ! fraction of local C panel size, green only
cs_fs_rows       = 0             ! Phase B; 0 means set from cs_panel_size
```

A run is rejected when:
- `control` is chosen and a generalized mode file is present;
- `cs_velocity = source` with `ISOLV = 2` or `IRSP = 1` (the existing rules);
- `cs_shape = file` and `CSMesh_{b}.pnl` is missing for any body.

Output:
- `Output/Wamit_format/MeanDriftControlSurface.7`, in the format of `MeanDriftPressure.9`: `PER BETA1 BETA2 I Mod Pha Re Im`, index 6(b−1)+I, scaled by ρgA²L^k. Phase A writes I = 1, 2, 6 and Phase B adds 3, 4, 5.
- `MeanDriftTerms.txt` gains the control-surface terms per body: the Γc line integral, the C flux and C pressure terms, and in Phase B the Γ, F, seabed and hydrostatic terms.
- The log gives the number of C, Γc and F points, and the time spent on field rows per frequency.

## 8. Validation

Meshes and frequencies are those of the Snellius study, so results compare one-to-one:
- sharp and rounded cylinders: a = 10, draft 10, h = 100, m1 to m4 meshes;
- boxes: 20 × 10 × 5, h = 50, heading 30;
- bottom-mounted cylinder.

The reference is the far-field value on the finest mesh.

### 8.1 Unit tests

| Test | Expectation | Why |
| :--- | :--- | :--- |
| Incident wave alone, closed C, no body | All six components → 0 as C is refined | Momentum flux in equals flux out. This checks the signs of n, n̄ and η, and the quadrature |
| Stencil on the incident potential | Differenced ∇φ_I matches the analytic value to O(δ²) | Checks the stencil and the z = 0 rule |
| Differenced field rows against `source_panel_integrals` | d(row_g)/dx matches ∫∇_p G for the same panel and point, including near panels | Checks the near-or-far rule (section 4.1) and the smoothness of G (section 4.4) |
| Green against source route at C points, rounded hull, away from irregular frequencies | Agree to a tolerance set from the measured difference | Two independent routes to the same field |
| `cs_fd_step` from 0.005 to 0.1 | A plateau of the drift values. The default goes in the middle of it | Sets the default |
| Generated C | Areas and normals exact on a cylinder and a box; another body inside C is caught | Mesh generator |

### 8.2 Phase A integration tests

| Case | Expectation |
| :--- | :--- |
| Independence of C: cylinder against box, margin 3, 5 and 8 panel sizes, rounded and sharp cylinder | Same F_x, F_y and M_z within the mesh tolerance. Gauss's theorem makes the exact result independent of C [paper, §3.2] |
| C made large | Approaches the far-field values. Chen shows (39-41) become the far-field form (19-21) at infinity [paper, p. 69] |
| Bottom-mounted cylinder, a = 5, h = 10 | F_x against MacCamy-Fuchs (`maccamy_fuchs_fx` in `drift_convergence.py`); F_y and M_z zero |
| Sharp cylinder, fixed and free, 0.40 to 1.80 rad/s | Against far-field, m1 to m4 |
| Sharp box, fixed and free, heading 30, 0.80, 1.00 and 1.40 rad/s | F_x, F_y and M_z against far-field, m1 to m4 |
| Rounded cylinder and rounded box | Against far-field and the converged near-field |
| Two bodies (`drift_two_round`) | Sum of F_x and F_y over the bodies against the far-field total. Each body against near-field |
| Bodies far apart (`drift_far_apart`) | Each against the isolated body |
| `IRSP = 1` against `IRSP = 0` (`drift_irsp`) | Agree away from irregular frequencies. At 1.548 and 1.939 rad/s, green + `IRSP = 1` against far-field + `IRSP = 1`, on the free rounded and sharp cylinders |
| Half mesh against full mesh | Same values |
| `ISOLV = 2` against `ISOLV = 1` (green route) | Same to solver tolerance |
| `mean_drift = 2`, diagonal pairs | Equal the `mean_drift = 1` values |

**Phase A gate.** F_x, F_y and M_z by `cs_velocity = green` must be within 1% of far-field on the m2 meshes (1760 to 1920 panels) for the sharp cylinder and the sharp box, at every study frequency. Where the far-field value is small, the error is measured against 0.001 ρgA²a, with a the cylinder radius or the box half-length (0.01 ρgA² for a = 10 m). This covers the free cylinder at 0.40 to 0.70 rad/s, where |far| is 3e-6 to 0.11 ρgA² against 3 to 6 ρgA² at higher frequencies.

### 8.3 Phase B integration tests

| Case | Expectation |
| :--- | :--- |
| Rounded cylinder and rounded box, fixed and free | F_z, M_x and M_y against near-field on the finest mesh. Near-field converges to within 0.24% on rounded hulls [measured] |
| Independence of C, as in Phase A | Same F_z, M_x and M_y |
| Sharp cylinder and sharp box | Converge with the mesh. There is no independent reference, so convergence is the result |

**Phase B gate.** Within 1% of near-field on the rounded hulls at m3. A failure points first to the open item of section 3.3.

### 8.4 Where the runs go

- The unit tests, and integration fixtures of a few hundred panels, run on the laptop as part of `fpm test`.
- The m1 to m4 sweeps go to Snellius with the existing `drift_convergence.py` driver, extended with a `control` column.
- Every Snellius job is put to the user first, with its node hours and memory estimated from the study's measured run times (`tables.md`, run-time table).

## 9. Cost

| Route | Green calls per frequency | Extra memory |
| :--- | :--- | :--- |
| green | about (6 N_C + 5 N_Γc + 5 N_F) × N_hull × nsys | O(N_hull nsys) per thread |
| source | about (N_C + N_Γc + N_F) × N_hull × nsys, plus the source factorization | the source operator, 64 nb² ns bytes (existing) |

Example [derived]: 7168 hull panels with no symmetry and N_C = 4000 give about 1.7e8 Green calls per frequency for the green route. That is about 3 times the 5e7 of the first-order influence matrix. Chen reports a small increase in run time for the source route [paper, p. 76]. Measure both on the m2 meshes before choosing any default other than `green`.

## 10. Risks

| Risk | Effect | Response |
| :--- | :--- | :--- |
| G not smooth in the field point (section 4.4) | Noisy velocity on C, finite depth | Merge the small-R fix first; unit test 3 |
| Cost of the green route | Slow runs on large meshes | Measure; the source route as fallback away from irregular frequencies; symmetry on C later |
| F next to the hull (Phase B) | Field points within a panel of the hull are less accurate. WAMIT warns about this [doc, ch. 11] | Grade F, measure, report |
| Hydrostatic term mismatch (section 3.3) | Wrong F_z, M_x, M_y | Resolve before Phase B |
| Sign of η or n̄ | Wrong signs of the Γc term | Incident-wave unit test |

## 11. References

- Chen, X.-B. (2007). Middle-field formulation for the computation of wave-drift loads. *J. Eng. Math.* 59:61-82. doi:10.1007/s10665-006-9074-x.
- WAMIT v7.5 manual: [chapter 11, control surfaces](https://www.wamit.com/manual7.x/html/wamit_v75manualch11.html), [chapter 15, eqs. 15.56 to 15.61](https://www.wamit.com/manual7.x/html/wamit_v75manualch15.html). The equations are images on the page; recheck them against the page before coding.
- Molin, B., Hairault, J.P. (1983). On second-order motion and vertical drift forces for three-dimensional bodies in regular waves. Proc. Intl Workshop on Ship & Platform Motions, Berkeley, pp. 344-362. Cited by Chen for (F_S, M_S); not read here.
- MacCamy, R.C., Fuchs, R.A. (1954). Wave forces on piles: a diffraction theory. US Army Beach Erosion Board, Tech. Memo 69.
- `docs/superpowers/specs/2026-11-04-second-order-waves-design.md`, sections D2, D6, D7 and Phase 3.
