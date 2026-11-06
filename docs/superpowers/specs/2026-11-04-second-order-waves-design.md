# Second-order wave loads in HAMS-MREL: solver review and implementation plan

Date: 2026-11-04. Branch reviewed: `Bugfixes` at `95950e6` (v2.3.1), then rechecked against `origin/main` at `58526e6`. The later commits up to `c40ad4d` (v2.3.2) changed only `readme.md` and `test/unit_tester.f90`. File and line references still hold there, apart from `src/HAMS_Prog.f90:735`, where issue 9 below has been fixed.
Status: proposal for review. No source code has been changed.

## 1. Summary

**Recommendation.** Build second-order loads in stages, from cheap and self-checking to expensive:

1. Mean drift by the far-field method. It needs only the potentials HAMS already has.
2. Mean drift per body by near-field pressure integration. This is the core of the project, because far-field cannot separate bodies in an array.
3. Difference-frequency QTFs, quadratic part. This reuses the near-field kernel.
4. Potential part of the QTF, using the indirect method.
5. Sum-frequency QTFs with free-surface integrals, only if a use case needs them.

Stages 1 and 2 deliver drift forces. Stages 1 to 4 deliver difference-frequency QTFs without free-surface integrals, an omission that Nemoh's documentation calls generally acceptable for difference frequencies.

**The three findings that shape the plan:**

| Finding | Consequence |
|---|---|
| HAMS solves a direct BIE for φ on constant panels. It has no body-surface velocity, whereas Nemoh and WAMIT both obtain velocity from a source distribution. | `∇φ` on the body is the main new numerical capability. It is the first near-field design decision (D2). |
| HAMS stores potentials and motions in different conventions: `MXPOT` holds φ/i for diffraction, and `DSPL` is in e^{+iωt} while the potentials are in e^{−iωt}. | Combining them naively gives wrong drift forces with plausible magnitudes. Phase 0 fixes the mapping and locks it with tests. |
| Nemoh is GPL-3.0-or-later and HAMS-MREL is Apache-2.0. WAMIT is proprietary. | Implement from published equations. Do not port Nemoh code (section 4.4). |

**Where HAMS-MREL can lead.** Nemoh's documentation states that its QTF computations "have not been tested for multi-body problems". Reliable per-body drift in multi-body arrays, with hydrodynamic interaction, fits the purpose of HAMS-MREL.

## 2. Scope

**Goals**

- Mean drift force and moment, per body and in total, for single and multi-body runs.
- Unidirectional and bidirectional waves.
- Finite and infinite depth, with both solvers (`ISOLV` 1 and 2), with and without irregular-frequency removal (`IRSP`), and with geometric symmetry.
- Difference and sum QTFs in later phases.
- WAMIT-format outputs that existing tools can read.
- No change in behaviour or output for runs that do not request second-order results.

**Non-goals**

- Forward speed and current.
- Near-field drift and QTFs with generalized modes. WAMIT 6.4S also excludes generalized modes from its second-order solution. Far-field drift can still support them (section 5.3).
- Fully nonlinear or time-domain second-order loads.
- Direct solution of the second-order potential φ⁽²⁾ on a free-surface mesh. Only the indirect method is planned.

**Assumptions to confirm** are listed in section 12.

## 3. The physics in brief

A second-order force has two parts:

- **Quadratic part.** Products of first-order quantities. It needs no new boundary-value problem.
- **Potential part.** The force from the second-order potential φ⁽²⁾, which has an inhomogeneous free-surface condition. It oscillates at ω₁ ± ω₂ and vanishes for the mean (ω₁ = ω₂).

So **mean drift is quadratic only**. Three methods exist:

| Method | Domain | Gives | Needs | Weakness |
|---|---|---|---|---|
| Far-field momentum (Maruo, Newman) | Large circle | F_x, F_y, M_z, total for all bodies | Kochin function from φ and ∂φ/∂n | No per-body split, no F_z, M_x, M_y |
| Near-field pressure (Pinkster, Lee) | Body surface and waterline | All 6 DOF, per body | ∇φ on the body, waterline elevation | Slow convergence, velocity singular at sharp corners |
| Control surface (momentum flux) | Surface around each body | All 6 DOF, per body | ∇φ in the fluid | User or automatic surface, more field evaluations |

## 4. Solver review

Evidence levels used below: **[source]** read directly in Nemoh's code. **[doc]** read in a manual or paper. **[derived]** worked out here and not yet run.

### 4.1 WAMIT

WAMIT is closed source, so this review covers the v7.5 manual and the v6.4S second-order chapter.

**Mean drift (v7.x)** [doc]

| Option | Method | Output | Notes |
|---|---|---|---|
| `IOPTN(7)` | Control surface | `.7` | All 6 components, per body. Recommended when pressure integration converges slowly. Needs a `.csf` file, or automatic generation. |
| `IOPTN(8)` | Far-field momentum | `.8` | I = 1, 2, 6 only. Makes no energy-conservation assumption, so it handles external damping. Uses adaptive azimuthal quadrature (`TOL` = 10⁻⁴, `MAXMIT`). |
| `IOPTN(9)` | Pressure integration | `.9` | All 6 components. Needs body velocity, so `ISOR=1` is required for the low-order method. |

Setting an option to 1 gives unidirectional results and 2 gives all heading pairs. Negative `IOPTN(6)` and `IOPTN(7)` select the source formulation in the low-order method.

Points to take from the manual:

- **Velocity needs sources.** With the low-order method, option 9 requires `ISOR=1`. With `ISOR=0`, the control-surface method gives reliable results only for F_x, F_y and M_z. This is the same problem HAMS has.
- **Cross-check by design.** The manual advises comparing option 9 with options 7 or 8.
- **Mesh guidance.** Velocity is singular at exterior corners, so panels must be finer near the free surface and corners.
- **Two control-surface alternatives** (`IALTCSF`). Alternative 1 moves the waterline integral to the control surface/free-surface intersection and is usually more accurate. Alternative 2 is for a surface very close to the body, such as a narrow gap.
- **Multi-body control surfaces.** Each body gets its own surface, which must not contain or intersect another body.
- **Output format.** `PER BETA1 BETA2 I Mod Pha Re Im`, scaled by ρgA²Lᵏ with k = 1 for forces and k = 2 for moments.
- **Bidirectional output** gives F(β₁,β₂) with F*(β₂,β₁) = F(β₁,β₂). The total mean load is ρgLᵏ[|A₁|²F(β₁,β₁) + |A₂|²F(β₂,β₂) + 2Re(A₁A₂*F(β₁,β₂))].

**Near-field formula (v7.5 eq. 15.54)** [doc]. WAMIT rotates body points with a matrix containing the quadratic term H and integrates pressure over the exact wetted surface. The result is a waterline integral, a `½|∇φ|²` integral, an `(ξ+α×x)·∇φ_t` integral, a `α×∬nφ_t` term and a hydrostatic term F_S. These equations were read from the manual's equation images, and the moment includes a `ξ×(−ρ∬nφ_t)` term, so M̄ is about a fixed point.

**QTFs (v6.4S, not in the v7 mainline)** [doc]

- `I2ND=1` in the CFG, and `ISOR=1` is **mandatory** for the low-order method.
- Generalized modes, walls and zero-thickness structures are excluded.
- Outputs: `.10s/.10d` quadratic by pressure integration, `.qcs/.qcd` quadratic by control surface, `.11s/.11d` total by the indirect method, `.12s/.12d` total by the direct method, `.13` to `.15` pressures and elevations, `.16s/.16d` second-order RAO.
- Indirect and direct methods must converge to each other, which gives another built-in check.
- When two periods are equal, the difference-frequency quadratic force equals the mean drift force.
- A `.pt2` file selects period and heading pairs, either all or a selected list.
- A `.fdf` file defines the free surface in three regions. The inner region is panelled (automatic generation is available). The intermediate annulus uses Gauss-Chebyshev azimuthally and Gauss-Legendre radially. The far field is semi-analytic with asymptotics.
- `R_part` should be O(h) in shallow water and O(λ) in deep water, where λ is the longest wavelength involved. The constant may need to be substantially larger than one in deep water.
- Omitting free-surface forcing (`RINNER<0`) is allowed as an approximation and is **not recommended for sum frequencies**.
- Scaling: F⁺/(ρgLAᵢAⱼ), F⁻/(ρgLAᵢAⱼ*), moments with L². Symmetry: F⁺ᵢⱼ = F⁺ⱼᵢ and F⁻ᵢⱼ = F⁻ⱼᵢ*.

### 4.2 Nemoh v3

Nemoh's QTF code is in `QTF/`, split into `PreProcessor`, `Solver` and `PostProcessor`, and runs after Nemoh1. I read `MQSolver.f90` (the part retrieved, about 60,000 characters, which covers all quadratic, body-force and free-surface terms and the asymptotic driver), `MQpreprocessor.f90` (about the first 40,000 characters), `MQSolverPreparation.f90` (about the first 30,000 characters), `MInfluenceMatrix.f90` and `Prepare_Waterline` in `Common/Face.f90`. I did not read `MQSolverAsymp.f90` or `MQpostproc.f90`.

**Pipeline** [source, doc]

1. `QTFpreProc` builds velocity on body panels from sources: `Velocity = MATMUL(gradS, σ)`, with total source strength σ = σ_diff + Σ(−iω ξⱼ σ_radj). The influence matrix `gradS` has extra rows for waterline points. Free-surface points are processed one row at a time to limit memory.
2. Potentials and velocities go to direct-access binary files, by frequency and heading.
3. `QTFsolver` loops over frequency and heading pairs and computes each term.
4. `QTFpostProc` sums the selected terms into `OUT_QTFM_N.dat` and `OUT_QTFP_N.dat`.

**Term structure** [source]. Option values select `DUOK`, `DUOK+HASBO` or the full `DUOK+HASBO+HASFS+ASYMP`.

| Name | Meaning | Terms |
|---|---|---|
| `DUOK` | Quadratic | Waterline (η−ζ)², body ½ρ∇φ·∇φ, body ρX·∇φ_t, R(F_I), optional translation moment and hydrostatic term |
| `HASBO` | Potential part, body forcing, indirect method | Second-order Froude-Krylov, ∂ₙφ_I⁽²⁾ψ, (Ẋ−∇φ)·R(n)ψ, two surface terms replacing second derivatives by Green's theorem, one waterline term |
| `HASFS` | Free-surface forcing on a finite mesh | Products of φ, ∇φ with ψ at ω±, in the region r < R_e |
| `ASYMP` | Free-surface forcing r > R_e | Fourier-Bessel expansion of Kochin coefficients with `N_Bessel` terms, default 30 |

Here ψ is the normalized radiation potential at the sum or difference frequency, obtained by **interpolation** from the first-order frequency list.

Implementation details worth noting:

- **Bi-harmonic helpers.** All quadratic terms go through `DOT_PRODUCT_DIFF_BIHARM` and `..._SUM_BIHARM`, then a global factor of ½. Nemoh's own comment queries that factor: `TRANSFER FUNCTION = BICHROMATIC FORCES / 2?`. Section 5.4 derives it.
- **Waterline points.** Segments are found from panel edges with both nodes near z = 0, and lid panels are skipped. The evaluation point is shifted off the body by `0.01·length` along the panel normal and placed at z = −ε·D. The documentation warns that the potential on the waterline "is rather sensitive with the epsilon value".
- **Inertia form.** R(F_I) uses F_I = F_exc + radiation − K_total ξ. The stiffness and damping include external matrices, so F_I includes external forces, not only fluid forces.
- **Waterline normal.** `ComputeNDGamma_WLINE` takes N(1:2) from the panel, sets N(3) to 0 and does **not** renormalize. On a flared hull this differs from a unit horizontal normal. Section 5.3 derives that the unit normal is correct, so do not copy this.
- **Hydrostatic term optional.** It is controlled by a flag in `Nemoh.cal` and, according to the documentation, is not always included in other software.
- **Frequency grid.** A uniform QTF grid is expected. Values are interpolated if they do not match. The manual recommends `ω_max,QTF ≤ ω_max/2`, so that sum frequencies stay in range, and first-order coverage of all difference and sum frequencies.
- **Output.** `ω1 ω2 β1 β2 DoF |Q|/ρg ∠Q Re/ρg Im/ρg`, lower triangle only.

**Stated limits** [doc]

- Multi-body QTFs are untested.
- Free-surface integrals have an unresolved issue with lid panels, so the advice is to leave lid panels out of full QTF runs. This matters because HAMS would have the same interaction with `IRSP=1`.
- A quadrilateral free-surface mesh is required for free-surface integrals.
- Compared with HydroSTAR, the authors note a sign difference in the imaginary part and need to swap the direction order for bidirectional cases, and they say this needs more investigation. So **convention errors affect even a mature code**.

**Paper** [doc]. Kurnia & Ducrozet (2023) use the near-field quadratic form of Chen (1988) and the indirect method of Molin (1979) and Chen (1988). Cases include a truncated cylinder with a heave plate (draft 10 m, radius 10 m, plate radius 15 m, 3464 panels including 262 lid panels), SOFTWIND (1872 panels) and OC4 (2196 panels). The SOFTWIND sum-frequency case uses a free-surface mesh of 8304 panels over a 50 m radius. They compare with HydroSTAR, DIODORE and WAMIT. The documentation's test cases `8b` and `9` to `11` can serve as references here.

### 4.3 Comparison

| Topic | WAMIT | Nemoh | Proposed for HAMS-MREL |
|---|---|---|---|
| First-order unknown | φ, optionally plus σ (`ISOR=1`) | σ | φ (direct BIE, existing) |
| Body velocity | From σ, or from higher-order geometry | From σ | Tangential differencing, plus optional source solve (D2) |
| Mean drift methods | Far-field, pressure, control surface | Near-field (QTF diagonal) | Far-field, pressure, then control surface |
| Multi-body | Yes | Untested for QTF | Yes, first-class |
| Generalized modes | Excluded in 6.4S | Not addressed in what I read | Far-field only |
| Irregular frequencies | Extended method | Lid, with issues for free-surface integrals | Existing `IRSP`, tested for drift |
| Quadratic rotation term | −ρα×∬nφ_t + F_S | R(F_I), optional HM term | WAMIT form (D4) |
| QTF potential part | Indirect and direct | Indirect | Indirect only |
| Free surface far field | Semi-analytic | Bessel expansion | Decide in Phase 6 |
| Second-order RAO | Yes (`.16`) | Not seen in what I read | Optional late item |
| Licence | Proprietary | GPL-3.0-or-later | Apache-2.0 |

### 4.4 What to adopt, and what to avoid

**Adopt**

- WAMIT's output file conventions, so BEMRosetta and OpenFAST-style readers can be used. Confirm each tool's parser before calling this compatible.
- WAMIT's practice of offering independent methods that must agree.
- Nemoh's term-by-term decomposition, with each term stored separately. It makes debugging far easier.
- Nemoh's one-kernel approach for sum and difference terms.
- Nemoh's progressive options (quadratic, then body forcing, then free surface).
- Row-at-a-time field evaluation for free-surface points, to limit memory.

**Avoid**

- A separate pre-processor, solver and post-processor chain of executables. HAMS has the first-order solution in memory, and one program is simpler for users.
- Non-renormalized waterline normals.
- Interpolating radiation data without a warning.
- Ambiguous normalization. Define everything in one place (section 5.4).

**Licensing.** Nemoh's source is GPL-3.0-or-later, and copying or translating it into Apache-2.0 HAMS-MREL would create a licence conflict. Implementers should work from the papers and the equations in this document, not from Nemoh's code. Reading Nemoh for the review and running it as an external reference are not the same as incorporating it. If this matters for release, get advice from TU Delft's licensing contact. I am not giving legal advice. Separately, `src/WavDynSubs.f90` is marked as part of LGPL-3.0-or-later FinGreen3D, so the repository already has mixed licences.

## 5. HAMS-MREL starting point

### 5.1 What already exists

| Need | Available | Where |
|---|---|---|
| Panel centroid, area, normal | Yes. Normals point out of the fluid (into the body). I checked this on all sampled example meshes: every panel normal points toward the body's node centroid. | `src/NormalProcess.f90` |
| Potentials on panels | Yes, per mode and per image | `MXPOT`, `MXPOT_MULTI_COMB` |
| Incident φ and ∇φ | Yes | `VINP`, `DINP` in `src/PatclVelct.f90` |
| Motions, added mass, damping, excitation | Yes | `DSPL*`, `AMAS*`, `BDMP*`, `EXFC*` |
| Field-point potential | Yes, but it recomputes influence for every mode | `CalPotentialMulti` (`src/PressureElevationMulti.f90:55`) |
| Green function and first derivatives | Yes. Derivatives are with respect to the **source** point. | `INFGREEN3D`, `FINGREEN3D` |
| Rankine panel integrals of 1/r and ∇(1/r) | Yes, on-body and off-body | `SGLINTBD_*`, `SGLINTBD_*2` |
| Second derivatives of G | **No** | Needed only for velocity from the direct BIE, which this plan avoids |
| Body velocity ∇φ | **No** | New |
| Panel adjacency | **No** | New |
| Waterline | **No** | New, and feasible: on all sampled meshes, every z = 0 edge is an open edge |
| Multi-frequency storage | **No**. `MXPOT` is overwritten for each heading and frequency. | New for QTFs |

### 5.2 Conventions [derived from code; to be locked by Phase 0 tests]

Let φ be the physical potential with time factor e^{−iωt}, and let the incident wave be φ_I = −i(gA/ω)Z(z)e^{ik(x cosβ + y sinβ)}, which gives η = Ae^{ik·x}.

| HAMS quantity | Meaning | Evidence |
|---|---|---|
| `VINP` = φ_I/i | Diffraction column stores φ/i | `src/PatclVelct.f90:68-74` |
| `MXPOT(:,6N+1,:)` | Scattered potential φ_S/i, for both `ISOL` values | `src/AssbMatxMulti.f90:214-216`, source term signs |
| `MXPOT(:,j,:)` | ψⱼ, with ∂ψⱼ/∂n = nⱼ, in e^{−iωt}, with no factor of i | A + iB/ω = ρ∬ψn from `RFORCE_MULTI` |
| `EXFC` | −iF, where F = ∬pn in e^{−iωt} | `src/PotentWavForceMulti.f90:108` |
| `DSPL` | ξ in **e^{+iωt}** | `LEFT`/`RIGHT` in `src/SolveMotionMulti.f90:83-85`, where `RIGHT` = conj(F) |

The WAMIT-format output of elevation in `WamitNondimensMulti` gives conj(η) under this mapping, and its radiation branch applies −iω and then a conjugate. Both support the table.

So the physical first-order total potential is:

```
φ_total = i (Φ_I,HAMS + Φ_S,HAMS) + Σⱼ (−iω) conj(DSPL(j)) ψⱼ
ξⱼ      = conj(DSPL(j))     (e^{−iωt})
∂φ_P/∂n = −∂φ_I/∂n + Σⱼ (−iω) ξⱼ nⱼ    (perturbation φ_P = φ_S + φ_R)
```

Internal second-order code uses e^{−iωt} throughout. Output code conjugates once, at the end, for WAMIT format.

### 5.3 Formulae [derived; to be verified in Phases 1 and 2]

Throughout, ⟨Re(a e^{−iωt}) Re(b e^{−iωt})⟩ = ½Re(a b*), and n points out of the fluid.

**Far-field.** HAMS's own field formula is 4πφ(p) = ∬(G ∂φ/∂n − φ ∂G/∂n_q)dS (`CalPotentialMulti`), with G's wave part tending to 2πi C Z(z)Z(ζ)H₀⁽¹⁾(kR). Define:

```
C = k² / ((k²−ν²)h + ν)        (deep water: C = k)
H(θ) = ∬ (∂φ_P/∂n − φ_P ∂/∂n) [Z(ζ) e^{−ik(x cosθ + y sinθ)}] dS
a(θ) = −(ω C / 2g) √(2/(πk)) e^{−iπ/4} H(θ)
```

Then the outgoing elevation is η_P ≈ a(θ)e^{ikr}/√r. A momentum and energy balance on a large circle gives, for any horizontal unit vector e:

```
F·e = (ρg/2)(C_g/C_p) [ X (e·e_β) − ∫₀^{2π} |a|² (e·e_θ) dθ ]
X   = −2 √(2π/k) Re[ A* a(β) e^{iπ/4} ]
```

Here C_g/C_p = ½ in deep water and ½(1 + 2kh/sinh 2kh) in finite depth. Checks already done by hand:

- A lossless body has X = ∫|a|²dθ, which recovers the Maruo form (ρg/2)(C_g/C_p)∫|a|²(1−cos(θ−β))dθ.
- For full reflection in the plane-wave analogue, the deep-water result is ½ρgA² per unit width, the known value.

Generalized modes need nothing beyond the modified nⱼ in the expression for ∂φ_P/∂n, so far-field can support them.

**Near-field mean force and moment** (WAMIT v7.5 eq. 15.54 to 15.57, read from the equation images on 2026-11-04). All coordinates are relative to a body origin O, taken in HAMS as `XR` of that body, and Z_o is the z coordinate of O. Angle brackets denote the time average.

```
F̄ = ½ρg ∫_WL n' ⟨ζ² − 2ζ X₃⟩ dl                                        (I)  waterline
    + ½ρg k ∫_WL [n_z/√(1−n_z²)] ⟨(ζ − X₃)²⟩ dl                        (Ia) flare, vertical only
    − ρ ∬_S n ⟨½∇φ·∇φ⟩ dS                                              (II) velocity
    − ρ ∬_S n ⟨X·∇φ_t⟩ dS                                              (III) displacement × gradient
    + ⟨α × (−ρ ∬_S n φ_t dS)⟩                                          (IV) rotation
    + F_S                                                              (V)  hydrostatic

M̄ = same with n' → x×n', k → x×k, n → x×n, IV → −ρ⟨α × ∬(x×n)φ_t dS⟩,
    plus ⟨ξ × (−ρ ∬_S n φ_t dS)⟩ and M_S instead of F_S

X₃ = ξ₃ + α₁y − α₂x,   X = ξ + α×x,   ζ = (iω/g) φ_total at the waterline
```

n' is the **unit** horizontal normal at the waterline, pointing into the body. The wetted strip has slant height dz/√(1−n_z²), so the horizontal part of n dS is n' dz dl. F_S and M_S are given in the plan (Task 7).

The waterline bracket is ζ² − 2ζX₃, not (ζ − X₃)², and term IV uses only the dynamic pressure. Together these are exactly equivalent, for the horizontal force, to the common form with (ζ − X₃)² and a rotation term that includes the first-order hydrostatic force. The difference is ∮n'X₃²dl, which Stokes' theorem turns into the cancelling hydrostatic part of the rotation term [derived]. This form needs no external loads and no inertia, so D4 follows from it directly.

**Far-field yaw moment and heading pairs** [derived]. With κ = (ρg/2)(C_g/C_p), and a₁ and a₂ per unit wave amplitude at the same frequency:

```
M_z,12 = −(κ/k) [ ∫ (1/2i)(a₂* ∂a₁ − a₁ ∂a₂*) dθ
                  + √(2π/k) ( i conj(∂a₂(β₁)e^{iπ/4}) − i ∂a₁(β₂)e^{iπ/4} ) ]
F_12   =  κ [ −√(2π/k) ( e_β₁ conj(a₂(β₁)e^{iπ/4}) + e_β₂ a₁(β₂)e^{iπ/4} ) − ∫ a₁ a₂* e_θ dθ ]
```

Total mean load is |A₁|²F₁₁ + |A₂|²F₂₂ + 2Re(A₁A₂*F₁₂). For β₁ = β₂ these reduce to the single-heading forms, and they satisfy F₂₁* = F₁₂. The moment's structure was fixed by requiring the correct change under a shift of reference point. This pins everything except a shift-invariant term, so it must also be checked against near-field M_z.

### 5.4 One kernel for mean drift and QTFs [derived]

For two wave components i and j, products of first-order signals expand as:

```
f(t) g(t) = ½ Σᵢ Σⱼ Re[ fᵢ gⱼ* e^{−i(ωᵢ−ωⱼ)t} + fᵢ gⱼ e^{−i(ωᵢ+ωⱼ)t} ]
```

For a force term c·f·g, define per unit wave amplitude:

```
T⁻ᵢⱼ = (c/2) · ½ (fᵢ gⱼ* + fⱼ* gᵢ)
T⁺ᵢⱼ = (c/2) · ½ (fᵢ gⱼ  + fⱼ  gᵢ)

F(t) = Σᵢ Σⱼ Re[ T⁻ᵢⱼ aᵢ aⱼ* e^{−i(ωᵢ−ωⱼ)t} + T⁺ᵢⱼ aᵢ aⱼ e^{−i(ωᵢ+ωⱼ)t} ]
```

The diagonal T⁻ᵢᵢ = (c/2)Re(fᵢgᵢ*) is the mean drift per unit amplitude squared. This answers the factor of ½ queried in Nemoh's comment. Writing the near-field code around this kernel from the start means mean drift, bidirectional drift and both QTFs share one tested implementation.

### 5.5 Existing issues that affect this work

Found during review. The first four can produce wrong second-order results and belong in Phase 0.

| # | Issue | Location | Effect here |
|---|---|---|---|
| 1 | `ISOL=2` drops the incident part from `EXFC`, so motions are wrong for this purpose | `EFORCE*` | Reject `ISOL≠1` when second-order output is requested |
| 2 | `VINP` uses the deep formula when `H>500`, while `DINP` uses it only when `H<0`. `SINH(k*H)` also overflows at large kH. | `src/PatclVelct.f90:64`, `:106`, `:112` | φ and ∇φ can be inconsistent, or NaN, in deep finite-depth cases. Far-field and near-field both depend on them. |
| 3 | `DISPERSION` has its own `G=9.807` while the rest of the code uses `9.80665` | `src/WavDynSubs.f90:58`, `src/WavDynMods.f90:146` | About 3×10⁻⁵ mismatch between k and ν. Small, but it will show in tight tests. |
| 4 | Generalized-mode files overwrite `DXYZ(:,dof)` | `CalNormals*` | Rigid-body displacement X = ξ + α×r no longer holds. Stop with a clear message for near-field. |
| 5 | Example meshes are not always conforming. `example5` has 29 open edges away from the waterline (59 open in total, 30 at z = 0). | Mesh data | Edge-based adjacency alone fails. D2 needs a tolerance-based node or proximity search. |
| 6 | `GRN_KIND=4` is the default, and `ACA_TOL`/`ITER_TOL` are 10⁻⁴. `developer_documentation.md` calls `GRN_KIND=8` "the default". | `src/WavDynMods.f90:47`, docs | Differencing amplifies solver noise by about 1/h, and drift then squares it. Measure in Phase 2. |
| 7 | `IRSP` normal equations use AᵀA and Aᵀb, not AᴴA and Aᴴb | `src/AssbMatx_irrMulti.f90:165`, `:309`, and the single-body file | This solves a consistent system but is not a true complex least-squares fit. It is outside this project, but check that drift near irregular frequencies is clean. |
| 8 | `Motion.4` and `ExcForce.3` scaling uses `REFL^2, ^3` for motions and no length scale for forces, which appears inconsistent with WAMIT when `REFL≠1` | `src/SolveMotionMulti.f90:123-129`, `src/PotentWavForceMulti.f90:113-115` | Second-order code must use in-memory arrays, never re-read these files. Raise it as a separate question. |
| 9 | Fixed on `origin/main`: `ReorderAmssDamp` had a hard-coded output path | `src/HAMS_Prog.f90:735` | None. New code should use `outputdir` in the same way. |
| 10 | `VINP` assigns the global `CI` on each call and prints when z > 0 | `src/PatclVelct.f90:59-60` | Called from OpenMP regions and with waterline points. Move to a pure helper for new code. |

## 6. Key design decisions

### D1. Far-field first

Far-field needs no velocity and no waterline, and uses only quantities the direct BIE gives accurately. It provides the reference for near-field work on any geometry, without another solver. **Recommended.** Its limitation is that for arrays it gives only the total force.

### D2. Body-surface velocity

| | T: tangential differencing | S: auxiliary source solve |
|---|---|---|
| Method | ∇φ = n ∂φ/∂n + least-squares surface gradient over neighbouring panels | Solve (2πI + D)σ = 4π ∂φ/∂n, with D = ∂G/∂n_p, then ∇φ = (1/4π)(self term + ∬σ∇_pG) |
| New solves | None | One extra LU per frequency, about double the factorization time |
| New data | Adjacency | D and ∇_pG from stored derivatives (with respect to p: horizontal parts flip sign, the wave and image z parts do not, the 1/r z part flips). |
| `ISOLV=2` | Works immediately | Needs its own H-matrix path |
| `IRSP=1` | Works immediately, since φ comes from the existing solve | Needs the waterplane extension to S as well |
| Weakness | One-sided at waterline and sharp edges. Sensitive to noise. Not the method in either reference code. | S and the direct BIE give slightly different φ. Velocity at panel centroids is still non-smooth. |
| Extends to the field | No | Yes: control surface and free surface need only first derivatives of G |

**Recommendation: implement T first, then S.** T gets near-field running across all solver options quickly. S is what WAMIT and Nemoh use, and it is required later for free-surface velocities. Choose the default by a measured convergence study against far-field on the same meshes (Phase 2 gate). Do not assume in advance which wins.

Avoid differentiating the direct BIE for velocity. That needs second derivatives of G and hypersingular integrals, which do not exist in the code.

**Details for T:**

- Neighbours are panels that share a node, plus a proximity search within about 1.5·`PNSZ` for non-conforming meshes.
- Exclude neighbours across a sharp edge (normal angle above a threshold, initially 45°).
- Project offsets into the panel's plane, and fit a 2D linear function with distance weighting.
- If a panel has too few neighbours, fall back to the normal component only and report a count in the log.

**Decision (2026-11-06): default T, S opt-in.** Both are implemented. The convergence study is in [mean_drift.md](../../mean_drift.md#convergence-study). Measured:

- Phase 2 gate met by T: on the free rounded cylinder at 7680 panels, near-field is within 0.24% of far-field from 0.70 to 1.80 rad/s. T is closer than S on the rounded cylinder at every frequency from 0.60 rad/s, and on the rounded box. On the bottom-mounted cylinder S is closer at 0.40 to 0.80 rad/s, with both within 0.06%.
- S converges slowly on the rounded cylinder: its error falls with refinement, local order about 0.4 to 0.6, and reaching zero is not shown.
- On the sharp cylinder S is closer up to 1.20 rad/s (ka about 1.5), and T from 1.40 rad/s. This rests on one radius, and on the sharp box S is closer at all three frequencies run, so there is no switch of method by frequency.
- S costs +98% run time at 7680 panels, cannot use `IRSP=1`, and is off by 32% to 198% at the two irregular frequencies tested.
- On the sharp box neither is reliable for F_y and M_z (+12.5% F_y with S at 7040 panels). Rounding the corners and the bilge brings T to 0.41%.

### D3. Waterline elevation

| Option | Method | Cost |
|---|---|---|
| E: extrapolate | φ_wl = φ_c + ∇φ·(x_wl − x_c) from the adjacent panel | None |
| B: BIE field point | Evaluate the field formula at the segment midpoint, offset slightly into the fluid | N_wl × N Green evaluations per frequency, once per point if the influence row is reused |

**Recommendation: implement both**, with E first since it is free. Reference solvers show this is delicate: Nemoh's documentation warns about sensitivity to its offset, and WAMIT warns about values within a panel dimension of the body. Expose the offset as a parameter and run a sensitivity study before fixing a default.

Add a `FieldInfluenceRow` routine that computes G and ∂G/∂n once per point and applies them to all modes. This avoids the repeated Green evaluations in `CalPotentialMulti`.

**Decision (2026-11-06): default E, B opt-in.** Both are implemented, with offset 0.01 of the segment length. Measured: B is first order and E second order on MacCamy-Fuchs, and E is closer on every mesh. On the free rounded cylinder at 1920 panels E is closer at all 11 frequencies. In finite depth B can jump with the mesh or the offset (term I -11.40% against about -2% at neighbouring offsets). A separate measurement on commit 36d1733, not reproduced on the drift branch, found that a direct eigenfunction sum for G takes that value to -1.96%. A fix is planned on a separate branch. B stays opt-in with a logged warning in finite depth. See [mean_drift.md](../../mean_drift.md#limits).

### D4. Rotation term

Use WAMIT's form: rotate only the dynamic pressure force, −ρ∬nφ_t, integrated from φ_total on the panels. Do **not** use −ω²Mξ or Nemoh's form with external matrices. With moorings, PTO or other external loads, the forces acting on the wetted surface are what get rotated, and the external force is not a pressure on that surface. For a free body with no external loads, the forms agree, and this gives a test.

### D5. Input

Use a new optional file, `SecondOrder.in`, with free-format `key = value` lines. If the file is absent, behaviour is unchanged.

The fixed-column `ControlFile.in` parser is fragile, and appending to it risks breaking existing inputs and test data. This also matches the existing pattern of optional `gen_mod` files detected by `INQUIRE`.

```text
! SecondOrder.in  (all keys optional unless stated)
mean_drift       = 1            ! 0 off, 1 unidirectional, 2 all heading pairs
drift_methods    = far near     ! far, near, control
velocity_method  = differencing ! differencing, source
waterline_method = extrapolate  ! extrapolate, field
waterline_offset = 0.01         ! fraction of segment length, field method only
motion_state     = free         ! free, fixed
n_azimuth        = 360          ! starting points for far-field quadrature
far_tolerance    = 1.0e-4
write_terms      = 1            ! 1: also write each near-field term
wl_z_tolerance   = 1.0e-6       ! relative to body size
! Later phases
qtf              = 0            ! 0 off, 1 diff, 2 sum, 3 both
qtf_terms        = quadratic    ! quadratic, body, full
```

`motion_state = fixed` sets ξ = 0 for second-order calculations only. It gives diffraction-only drift, which has analytic references.

### D6. Output

Write to `Output/Wamit_format/`, consistent with `AmssDamp.1` and `ExcForce.3`.

| File | Contents | Phase |
|---|---|---|
| `MeanDriftMomentum.8` | Far-field, I = 1, 2, 6, system total | 1 |
| `MeanDriftPressure.9` | Near-field, per body with index 6(b−1)+I | 2 |
| `MeanDriftControlSurface.7` | Control surface, per body | 3 |
| `MeanDriftTerms.txt` | Terms I to V per body | 2 |
| `QTF_Quadratic.10d`, `.10s` | Quadratic part | 4 |
| `QTF_Indirect.11d`, `.11s` | Total, indirect method | 5, 6 |

Format: `PER BETA1 BETA2 I Mod Pha Re Im`, scaled by ρgA²Lᵏ. Use the existing `OUFR` for the first column, and document this, since WAMIT uses period. Open files with `NEWUNIT`, and use `outputdir`. Confirm the conjugation with a WAMIT reference case. Note that readers may assume a specific form of index for multi-body data.

### D7. Single-body and multi-body paths

HAMS duplicates nearly every routine (`X` and `XMulti`). **Do not duplicate again.** Write new modules on explicit arrays with an `NBODY`-independent interface. Add one small adapter each for the single-body arrays (`XYZ_P`, `DXYZ_P`, `DS`, `MXPOT`, `MD = 7`) and the multi-body arrays (`*_MULTI_COMB`, `MD = 6N+1`). `MXPOT(…,IP)` holds the value on image IP in both, so symmetry handling is shared.

### D8. Storage

Mean drift needs only the current frequency, so it runs inside the existing loops. It stores diffraction potentials across headings within each frequency, and only when bidirectional drift is requested. QTFs need all frequencies. Store `φ`, `∇φ` and `ψ` per frequency in a direct-access file, or in memory when small.

Rough size for 5000 panels, 40 frequencies, 3 headings and 6 bodies: about 38 MB for total fields and about 460 MB for radiation fields (double-precision complex). This is modest compared with the Green's function arrays.

## 7. Architecture

New files in `src/`. `fpm` compiles everything in `src/`, but `CompileHAMS.bat` and `HAMS_MREL.vfproj` must be updated by hand (see Phase 0).

| Module | Responsibility | Depends on |
|---|---|---|
| `SecondOrderTypes.f90` | Options, first-order state view, result types | none |
| `SecondOrderInput.f90` | Read and validate `SecondOrder.in` | Types |
| `FirstOrderState.f90` | Convention conversion (section 5.2), φ, ∂φ/∂n, ξ, α, F⁽¹⁾ for a frequency and heading | Types, wave modules |
| `SurfaceGeometry.f90` | Adjacency, sharp-edge flags, waterline segments per body | Panel arrays only |
| `SurfaceVelocity.f90` | ∇φ by method T or S | Geometry, state |
| `FarFieldDrift.f90` | Kochin function, adaptive quadrature, F_x, F_y, M_z | State |
| `NearFieldDrift.f90` | Terms I to V using the bi-harmonic kernel | Geometry, velocity, state |
| `WaterlineField.f90` | Field-point influence row and BIE waterline potential | Green's function modules |
| `SecondOrderOutput.f90` | WAMIT-format writers | Types |
| `SecondOrderDriver.f90` | Single and multi-body adapters, orchestration | The above |
| `QTFDriver.f90` (later) | Frequency and heading pair loops, storage | The above |

**Integration in `HAMS_Prog.f90`.** The call sites are clear: after `SolveMotion` (line 300) and `SolveMotionMulti` (line 724), inside the heading loop, and before `MXPOT` is overwritten. Skip when `TP` flags the zero or infinite frequency limits (`SYBO=1`).

```text
read inputs → read SecondOrder.in if present → validate combinations
build surface geometry once (after BuildCombinedArrays / CalNormals)
for each frequency:
    solve radiation; for each heading: solve diffraction, motions
        state = view of current first-order solution (converted)
        far-field → accumulate
        near-field: velocity → terms → accumulate
write outputs
```

**Validation at startup.** Stop with a clear message for:

- `ISOL≠1`.
- Near-field with `gen_mod` files.
- Near-field on a surface-piercing body with no waterline found.
- `velocity_method=source` with an unsupported solver combination, until supported.
- `qtf` with a frequency grid that cannot provide the pairs.

## 8. Implementation phases

Duration estimates are rough engineering judgement for one developer who knows the code. Acceptance thresholds are proposals to be adjusted after the first results.

### Phase 0: Foundations (about 1 week)

1. Add `FirstOrderState` and lock the section 5.2 conversion with tests.
2. Fix or guard section 5.5 items 1 to 4. Agree separately whether issue 2 changes first-order behaviour, because it would change existing results for `H>500`.
3. Add `SecondOrder.in` parsing and startup validation.
4. Add new modules to `CompileHAMS.bat` and `HAMS_MREL.vfproj`, which both list source files explicitly, in dependency order where needed.
5. Set up test geometry: a cylinder mesh generator (Python), at several resolutions.

**Done when:**

- With no `SecondOrder.in`, the existing integration tests pass unchanged.
- The equation of motion residual using converted quantities is below 10⁻⁸ relative, on an example case.
- Elevation from `φ_total` reproduces `PressureElevation.6p` after conjugation.

### Phase 1: Far-field mean drift (1 to 2 weeks)

1. Kochin function for total perturbation and for each radiation mode, with symmetry images.
2. Adaptive azimuthal quadrature, with a stationary-phase cross term.
3. F_x, F_y, then M_z. Then all heading pairs.
4. `MeanDriftMomentum.8`.

**Tests**

| Test | Expected |
|---|---|
| Bottom-mounted fixed cylinder a(θ) vs MacCamy-Fuchs series | Converges with mesh, within about 1% at a moderate mesh |
| Drift on the same cylinder vs analytic (ρg/2)(C_g/C_p)∫|a|²(1−cosθ)dθ | Within about 2% |
| Radiated energy: ½ω²Bⱼⱼ = ½ρgC_g∫|aⱼ|²dθ against `BDMP` | Within about 1% |
| Haskind: excitation from H_j(β+π) against `EXFC` | Within about 1% |
| Energy balance X = ∫|a|²dθ on a free body with no external damping | Within about 1% |
| M_z about a shifted reference point | M_z changes by the lever arm × F |
| Symmetric body at symmetric heading | F_y and M_z near zero |
| F(β₁,β₂) = F*(β₂,β₁) | To round-off |

### Phase 2: Near-field mean drift, per body (3 to 4 weeks)

1. Adjacency and waterline extraction.
2. Velocity method T. Waterline method E.
3. Terms I to V, in the bi-harmonic kernel, one test per term.
4. Moments, with an explicit choice of origin convention.
5. `MeanDriftPressure.9`, `MeanDriftTerms.txt`.
6. Method B for the waterline, and the offset sensitivity study.
7. Velocity method S. Direct solver first, then `IRSP`, then `ISOLV=2` if needed.
8. Convergence study: pick default velocity and waterline methods.

**Tests**

| Test | Expected |
|---|---|
| Waterline length of a cylinder | 2πa within mesh error. Closed loop per body. Correct orientation. |
| Uniform and linear φ: surface gradient | Exact on planar neighbours |
| Incident-wave velocity from differencing vs `DINP` | First-order convergence or better |
| Fixed cylinder, near vs far | Within about 5% at a moderate mesh, and converging |
| Free truncated cylinder, near vs far | Within about 5% for F_x, away from resonances |
| Term V and D4: no external loads, form with −ω²Mξ vs fluid-force form | Equal to round-off level of the motion solve |
| With external stiffness: near vs far | Agreement holds using the fluid-force form |
| Two distant bodies | Each approaches the isolated value |
| Two bodies, sum of near-field F_x, F_y vs far-field total | Within about 5% |
| Half mesh with symmetry vs full mesh | Same to mesh-independent tolerance |
| `ISOLV` 1 vs 2, `GRN_KIND` 4 vs 8, `IRSP` 0 vs 1 | Differences quantified and documented |
| Reference run of Nemoh diagonal (cylinder) and WAMIT `.9`, if available | Agreement with sign and normalization understood |

**Gate.** Do not start Phase 4 until near-field and far-field agree and converge on at least one free body. If method T cannot meet this, S becomes the default.

### Phase 3: Control surface (optional, about 2 weeks)

Add this when near-field proves inaccurate on cases that matter, typically sharp edges or narrow gaps. WAMIT's manual recommends it for that purpose. It needs S-based field velocity, or central differences of the field potential, which are acceptable away from the body and need no second derivatives. Start with an automatic circular cylinder around each body. Multi-body surfaces must not contain another body.

Chen's middle-field formulation (2007) is a related alternative and could be reviewed at that point.

### Phase 4: QTF, quadratic part (about 3 weeks)

1. Reuse the Phase 2 kernel with ωᵢ ≠ ωⱼ.
2. Storage (D8). Option for all pairs or selected pairs, like WAMIT's `.pt2`.
3. `QTF_Quadratic.10d` and `.10s`.
4. Python post-processing: Newman's approximation from the diagonal, and slowly varying force spectra. This can be delivered earlier for users who need only slow-drift estimates.

**Tests**

- Diagonal equals mean drift to round-off.
- T⁻ᵢⱼ = T⁻ⱼᵢ* and T⁺ᵢⱼ = T⁺ⱼᵢ.
- A direct evaluation of a bichromatic quantity in time on a toy signal matches the kernel, which verifies the factor of ½.
- Comparison with Nemoh's `DUOK` output for its cylinder test case, after converting conventions. This is an external comparison only.

### Phase 5: QTF, potential part, body forcing (4 to 6 weeks)

1. Second-order incident potential and force (known analytic form; validate its free-surface condition numerically).
2. Radiation ψ and ∇ψ at Δω. This needs low-frequency solves. Add a list of "radiation-only" frequencies, so full diffraction is not repeated.
3. Body-forcing integrand, including conversion of second-derivative terms with Green's theorem.
4. `QTF_Indirect.11d`. Report clearly that the free-surface term is omitted.

**Tests**

- Result tends to the quadratic part as Δω → 0.
- Convergence with frequency spacing.
- Compare with Nemoh `DUOK+HASBO` on the cylinder and OC4 cases.
- Compare with WAMIT `.11d` with `RINNER<0`, if available.

**Risks.** Low Δω is long-wave, so the Green function must remain reliable and the first-order frequency list must cover it. Both references warn about interpolation here.

### Phase 6: Free-surface integrals and sum frequency (4 to 8 weeks, only on need)

Needs S-based velocities on free-surface points, a free-surface mesh (automatic if feasible), an inner integral and an outer integral. Sum frequency cannot omit this. WAMIT's manual advises against it, and Nemoh's documentation calls the free-surface integrals necessary. Interaction with `IRSP` must be studied, since Nemoh reports trouble with lid panels. Decide between semi-analytic asymptotics and Bessel expansion at that time.

Second-order motion RAO (WAMIT option 16) can be added here if needed.

## 9. Validation matrix

| Geometry | Depth | Case | Phases | Reference |
|---|---|---|---|---|
| Bottom-mounted cylinder, fixed | Finite | Scatter, drift | 1, 2 | MacCamy-Fuchs series (analytic) |
| Truncated cylinder, fixed and free | Deep, finite | All 6 DOF | 1, 2, 4 | Self-consistency, Nemoh `8b`, WAMIT if available |
| Hemisphere, free | Deep | Classic drift | 1, 2 | Published data: extract and verify before use. Kudou (1977) and Pinkster & van Oortmerssen (1977) are commonly cited, but I did not verify these numbers. |
| Box with sharp corners | Deep | Convergence stress test | 2, 3 | Far-field |
| Two or more cylinders | Finite | Per-body drift, gap effects | 2, 3 | Far-field total, WAMIT multi-body if available |
| `example2` array | Finite | Performance | 2 | Timing only |
| OC4 semi-submersible | 200 m | QTF | 4, 5 | Nemoh `9`, published figures |
| SOFTWIND | 200 m | Sum QTF | 6 | Nemoh `10b`, published figures |

**CI.** Add a small drift integration case to `test/data`, following `config1`, with results checked against a stored reference and against the far-field/near-field identity. Keep meshes small enough for runner time limits. Unit tests go in `test/` via `test-drive`, and no MKL-specific code should be needed for the pure geometry and kernel tests.

## 10. Risks

| Risk | Likelihood | Effect | Mitigation |
|---|---|---|---|
| Convention errors | High without Phase 0 | Wrong results that look plausible | Section 5.2, tests, far-field cross-check |
| Near-field convergence too slow with constant panels | Medium to high | Inaccurate per-body drift | Refinement guidance, S method, control surface |
| Differencing noise with `GRN_KIND=4` or GMRES | Medium | Noisy drift | Measured in Phase 2. Allow or recommend tighter settings for second-order runs. |
| Waterline sensitivity | High | Term I dominates error in many cases | Two methods, offset study, unit normal |
| Sharp edges at waterline (for example hull chines) | Medium | Local velocity errors | Sharp-edge handling, warnings like WAMIT's `VMAXOPT9` |
| Narrow-gap arrays | Medium | Near-field and CS both difficult, plus gap resonance | Document limits. Alternative 2-style CS in Phase 3. |
| Irregular frequency effects on drift | Medium | Spikes, since drift is quadratic | Test with `IRSP` 0 and 1 |
| Cost of QTF low-frequency solves | Medium | Long runs | Radiation-only frequency list |
| Licence contamination | Low if followed | Blocks release | Implement from equations, keep a record of sources |
| Two near-duplicate code paths | Certain if not designed against | Maintenance | D7 |

## 11. Effort overview

| Phase | Result | Rough effort | Depends on |
|---|---|---|---|
| 0 | Foundations | 1 week | none |
| 1 | Far-field drift, system total | 1 to 2 weeks | 0 |
| 2 | Near-field drift, per body | 3 to 4 weeks | 0, checked against 1 |
| 3 | Control surface | 2 weeks, optional | 2 (S) |
| 4 | Quadratic QTF | 3 weeks | 2 |
| 5 | Potential QTF, body forcing | 4 to 6 weeks | 4 |
| 6 | Free surface, sum frequency | 4 to 8 weeks, on need | 5 |

Drift forces are available after about 5 to 7 weeks. A difference-frequency QTF without free-surface integrals (through Phase 5, skipping Phase 3) takes about 12 to 16 weeks in total.

## 12. Open questions

1. **Scope.** Decided on 2026-11-04: drift first (Phases 0 to 2). Later phases stay in this document as direction and are not yet committed work.
2. **Target applications.** Are the main cases floating wind and WEC arrays with narrow gaps? This decides how soon the control surface is needed.
3. **External reference.** Decided on 2026-11-04: use Nemoh for now, as an external program. WAMIT can be added if a licence becomes available.
4. **Issue 2 in section 5.5.** May the `H>500` inconsistency be fixed in first-order code, even though it changes outputs for such cases? Should it be a separate pull request?
5. **Licence check.** Is a TU Delft review wanted given that Nemoh was read during this review?
6. **Output naming.** Prefer `MeanDriftPressure.9`-style names (matching `ExcForce.3`), or plain WAMIT names?
7. **Gravity and density.** Should `G` and `RHO` become inputs while this work is under way? They are hard-coded at present.

## 13. Sources

**HAMS-MREL.** Files cited above, at commit `95950e6`. Mesh checks were run on `examples/example1`, `example2` and `example5` with a throwaway script. No solver was built or run for this review, so every statement about run-time behaviour comes from reading the code.

**WAMIT**

- [WAMIT v7.5 manual, chapter 3: quantities and mean drift](https://www.wamit.com/manual7.x/html/wamit_v75manualch3.html)
- [Chapter 4: inputs](https://www.wamit.com/manual7.x/html/wamit_v75manualch4.html), [chapter 5: outputs](https://www.wamit.com/manual7.x/html/wamit_v75manualch5.html)
- [Chapter 11: control surfaces](https://www.wamit.com/manual7.x/html/wamit_v75manualch11.html), [chapter 15: theory](https://www.wamit.com/manual7.x/html/wamit_v75manualch15.html)
- [WAMIT v6.4 manual, chapter 11: second-order module](https://www.wamit.com/manual6.4/Chap11.pdf)

The v7.5 pages were read through an automated summary, except eq. 15.54 to 15.57, which were read from the page's equation images. The v6.4S chapter was read as extracted text. Chapter 3's summary places options 7 and 9 in body coordinates, but chapter 5's says forces are with respect to the inertial system. Resolve this from the manual before finalizing frames.

**Nemoh**

- [Repository](https://gitlab.com/lheea/Nemoh), read at `master` on 2026-11-03: `QTF/Solver/MQSolver.f90`, `QTF/Solver/MQSolverPreparation.f90`, `QTF/PreProcessor/MQpreprocessor.f90`, `QTF/PreProcessor/MInfluenceMatrix.f90`, `Common/Face.f90`
- Documentation: [description](https://lheea.gitlab.io/Nemoh/description.html), [input and output](https://lheea.gitlab.io/Nemoh/input-output.html), [test cases](https://lheea.gitlab.io/Nemoh/test-cases.html)
- Kurnia, R., Ducrozet, G. (2023). NEMOH: Open-source boundary element solver for computation of first- and second-order hydrodynamic loads in the frequency domain. *Computer Physics Communications* 292, 108885. [Zenodo record](https://zenodo.org/records/10362436)
- Kurnia, R., Ducrozet, G., Gilloteaux, J.-C. (2022). Second order difference- and sum-frequency wave loads in the open-source potential flow solver NEMOH. OMAE2022-79163.

**Theory** (cited from memory or from the references above, and not re-read for this review; verify details before relying on them)

- Maruo, H. (1960). The drift of a body floating on waves. *Journal of Ship Research* 4.
- Newman, J.N. (1967). The drift force and moment on ships in waves. *Journal of Ship Research* 11.
- Newman, J.N. (1974). Second-order, slowly-varying forces on vessels in irregular waves. Proc. International Symposium on Dynamics of Marine Vehicles and Structures in Waves.
- Pinkster, J.A. (1980). Low frequency second order wave exciting forces on floating structures. PhD thesis, TU Delft.
- Molin, B. (1979). Second-order diffraction loads upon three-dimensional bodies. *Applied Ocean Research* 1(4), 197-202.
- Chen, X.-B. (1988). Etudes des réponses du second ordre d'une structure soumise à une houle aléatoire. PhD thesis, Université de Nantes.
- Chen, X.-B. (2007). Middle-field formulation for the computation of wave-drift loads. *Journal of Engineering Mathematics* 59.
- MacCamy, R.C., Fuchs, R.A. (1954). Wave forces on piles: a diffraction theory. US Army Beach Erosion Board, Technical Memorandum 69.
