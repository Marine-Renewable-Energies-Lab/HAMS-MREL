# Damping lid on the free surface: design

Status: awaiting review. Base: `origin/main` at 6f4cc3b. Branch: `feature/damping-lid`.

## Goal

Let the user give a lid of any 2D planform on z = 0, anywhere on the free surface, on which the potential is zero. The lid may touch a body or lie away from every body. Results of runs without a lid do not change.

The user chose phi = 0 over a tunable dissipative lid. Where the lid sits is an input; how the equations are built stays as close to `IRSP` as the physics allows.

## What a phi = 0 lid does

Derived, not measured. With time factor e^{-iwt}, the mean power through the surface per unit area is proportional to Im(phi* dphi/dz), which is zero when phi = 0. The lid absorbs no energy. It suppresses the wave on the lid by reflecting it, as a pressure-release surface does.

The documentation must say this, although the file is called `DampingLid.pnl`. Energy balance and the Haskind relation therefore still hold with a lid, which the validation uses.

## Why the IRSP rows cannot be reused unchanged

`IRSP` adds rows and no unknowns (`AssbMatx_irr.f90`, JEL = 1..NELEM). Inside a body those rows state the null-field identity, which the true solution satisfies, so they are consistent with the hull rows. On open water the same rows state 4 pi phi(p) = 0, which the hull solution does not satisfy.

Measured with `IRSP = 1` and the waterplane mesh replaced by an open-water annulus. Binary from 0a77224 (`fix/image-term-integration`), not a build of main. Cylinder r = 5, draft 5, 180 quad panels, depth 50, 1.35 rad/s, annulus r = 6 to 10:

| Lid panels | A11 | B11 | A33 | B33 |
| ---: | ---: | ---: | ---: | ---: |
| none | 253.4 | 206.3 | 210.2 | 21.9 |
| 40 | 434.3 | -186.3 | 88.1 | 14.3 |
| 160 | -255.9 | -374.7 | 65.1 | 1327 |
| 640 | -337.2 | -102.6 | 602.9 | 239.1 |

The results do not converge with lid refinement, B11 is negative, and the runs finished with no warning.

## Formulation

Not yet tested. The wave Green function G satisfies dG/dz = nu G on z = 0, so the free-surface integral vanishes outside the lid. On the lid it leaves a single layer with the unknown

    q = dphi/dz - nu phi

The code's hull equation, with G differentiated at the source and the normal out of the fluid, gains that term. D is the double-layer and S the single-layer panel integral, j runs over hull panels and l over lid panels:

    hull row i:   2 pi phi_i + sum_j D_ij phi_j - sum_l S_il q_l = sum_j S_ij V_j
    lid row l:                 sum_j D_lj phi_j - sum_l' S_ll' q_l' = sum_j S_lj V_j + T_l

The lid row is the `IRSP` row (zero diagonal, field point on z = 0) plus one column for q. Unknowns: hull phi and lid q. Each lid panel adds one row and one unknown, so the system is square.

phi = 0 applies to the total potential:

- Radiation: phi = 0 on the lid, so T_l = 0.
- Diffraction: the code solves for the scattered potential, and phi_S = -phi_I on the lid, so T_l = 4 pi phi_I(p_l). q is the same for scattered and total, because phi_I satisfies the free-surface condition.

The incident term is what the `IRSP` rows lack. Without it the incident wave passes through the lid.

Signs come from reading the code's conventions. Check 7 below tests them against a reference that shares no code with the solver.

### Frequency limits

- nu = 0: G = 1/r + 1/r1, and the formulation holds with q = dphi/dz.
- Infinite frequency: G is zero on z = 0, so phi is already zero on the whole free surface and the lid rows vanish. The lid is skipped at this frequency and the hull system is solved alone.

## Input

- `ControlFile.in`: a new optional line `Damping_lid 1` after `Solver_type`, read the same backward-compatible way. Absent or 0 means no lid. The lid is off unless asked for, so a stray file cannot change a result.
- `DampingLid.pnl`, in the same format as `WaterPlaneMesh.pnl`: triangles or quads, all nodes on z = 0. One file for any number of bodies, in global coordinates, since a lid belongs to no body.
- Holes and separate patches need nothing special: the lid is a list of panels. Vertex order must not matter, since no normal is used; a test checks this on the panel integrals.

The run stops with a message when:

- the line is 1 and the file is missing (`VerifyInputFilesExist`)
- a node is off z = 0
- `ISX` or `ISY` differ from the hull's, in the lid file (it must not overwrite them, as the `WaterPlaneMesh.pnl` read does)
- `Solver_type = 2`
- a lid panel centroid lies inside a hull. Test: the sum of the hull double-layer Rankine integrals and their images is a solid angle at z = 0, about 4 pi in magnitude inside a closed hull and 0 outside. An interior lid is what `IRSP` is for.

## Panel integrals

Lid data needs its own arrays. The existing `IRR = 2` and `4` branches of `SGLINTBD_TRI2` and `SGLINTBD_QUAD2` read the waterplane arrays `iXYZ` and `iNCON`, which hold the `IRSP` lid when both are on. Add a new `IRR` value that reads the lid arrays, and the same in `SingularIntgrMulti.f90`. Hull field points and lid field points both use this branch.

`S` for a lid source is the sum of:

- **Rankine and image:** exact panel integrals. For a field point on z = 0 the image integral equals the direct one. `IMAGE_TERM_EXACT` cannot be used on the lid, since it divides by R1 = 0 at a panel's own centroid.
- **Wave term, panels apart:** one point at the centroid, as elsewhere in the code.
- **Wave term, self panel and field points inside a panel:** `HavelockGF` divides by the zero distance and takes the log of zero there. Split the panel into triangles at the field point and integrate the log singularity with a Duffy transform.

### Finite depth

The sibling measurement of `Limes` (`limes-measurement-result.md`) found the finite-depth G unreliable at small R/h, especially near the surface and with both points at the same depth. Scaled error in G: 68 at R/h = 0.0005, 2.0 at 0.005, 0.14 at 0.01, 7e-4 at 0.03. That is where lid-to-lid quadrature would sample. So for the self term:

- Integrate the deep-water wave term by the Duffy rule. It has no `Limes`.
- Add the remainder G_finite - G_deep, which has no log singularity if the two log coefficients match, as a smooth term sampled at R/h of at least 0.03.

The matching coefficient is an assumption until checked: the plan starts by evaluating the remainder near R = 0 against the direct series reference from that measurement.

Lid panels with centroid spacing below about 0.01 h stay exposed to the same scatter as hull panels, through the one-point term. The documentation will state this. Sibling work may change `FINGREEN3D`; this design does not depend on that.

## Solver

- New code only runs when a lid is on, so the no-lid paths are not edited in behaviour. New routines live in `src/DampingLid.f90`, plus a multi-body counterpart, and `ReadControlFile`, `VerifyInputFilesExist` and `HAMS_Prog.f90` gain small branches. Check at implementation whether `HAMS_MREL.vfproj` and `CompileHAMS.bat` list sources by name.
- Direct LU on the square system. Hull phi is copied to `MXPOT`, and lid q is stored in a new array with the same mode and symmetry layout.
- With `IRSP = 1` as well: interior rows gain the same q column, the system is rectangular, and it goes through normal equations as `IRSP` does. Those rows stay consistent, because a point inside the body is still outside the fluid. Rows of the first kind make squaring riskier, so check 8 measures it.
- Symmetry: `ISX` or `ISY` follow the existing `RXY` pattern. The half mesh must match the full mesh.
- GMRES is out of scope. The H-matrix assumes hull geometry, and lid rows would need their own preconditioner.

## Outputs

- `RFORCE` and `EFORCE` integrate `MXPOT` over the hull and need no change.
- `CalPotential` gains the lid term, so pressure and elevation at field points include the lid's effect, in `PressureElevation.f90` and `PressureElevationMulti.f90`. Field points inside a lid panel use the Duffy rule above.
- `ErrorCheck.txt` and the console report the lid panel count.
- Mean drift (`second-order-drift`) is out of scope.

## Scope and order

1. Single body, finite and infinite depth, symmetry, with `IRSP`.
2. Multi-body on the same branch, through the combined-mesh arrays. A lid in a gap between bodies is likely the main use. `example4` and `HAMS_Prog.f90` show a multi-body `IRSP` path, although a code comment says single body only; stage 2 relies on the combined hull path and reuses the `IRSP` path only after checking it.

No push or PR unless asked.

## Known numerical limits

Expected, to be measured in check 4:

- q is singular at the lid edge, roughly as 1/sqrt(distance), as for an electrified disc. Constant panels will converge slowly near the edge.
- Lid rows are of the first kind, so conditioning worsens with refinement.

## Validation

Measured values go in a fixtures README, with limits set from them, as in `test/data/image_term/README.md`.

| # | Check | Pass means |
| --- | --- | --- |
| 1 | No lid | Existing unit and integration tests pass with stored outputs unchanged |
| 2 | Self term at nu = 0, against the analytic electrified disc | The lid-to-lid operator reproduces the disc's known total charge, converging with refinement |
| 3 | Duffy rule and finite-depth remainder | Agree with refined quadrature and with the direct series reference |
| 4 | Convergence: the annulus above at 40, 160, 640 panels, opposite of the measured table | A, B, X converge, B >= 0; the rate is reported, along with the condition number |
| 5 | Lid condition | Total elevation about 0 at lid centroids; off-centroid values are reported |
| 6 | Lossless checks | A and B symmetric; Haskind holds, using the helper behind `haskind_heave` |
| 7 | Independent reference | Bottom-mounted cylinder with an annular phi = 0 lid, solved by eigenfunction matching in Python, with no BEM or Green function. Inside the lid the radial functions are I and K with vertical modes cos((n - 1/2) pi (z + h)/h); outside, the usual free-surface modes. Compare surge force and elevation. As b -> a it must reduce to MacCamy-Fuchs, which also checks the hull mesh on the sea bed. |
| 8 | Combinations | Half mesh matches full mesh; with `IRSP` matches without it away from irregular frequencies; shrinking lid tends to no lid |
| 9 | Stops | Each stop condition above exits non-zero with its message |
| 10 | Multi-body | Lid with a far second body tends to the single-body result; two-body A and B symmetric |

Tests stay small enough for the laptop: ifort 2022.1 and fpm, per `docs/developer_documentation.md`. Any larger convergence study goes to Snellius, with job sizes and cost shown first.

## Documentation

- `docs/input_files.md`: the `Damping_lid` line and `DampingLid.pnl`, with the reflecting behaviour, the GMRES stop and the small-panel caveat.
- `docs/developer_documentation.md`: the formulation and the measured results.
- `readme.md`: one line on the input.

## Literature to check before citing

Lid methods for gap resonance appear in Buchner et al. (2001), Huijsmans et al. (2001), Newman, and Chen (2004 or 2005; sources differ). These were found through search and are not read at the source. This design does not rely on them, but any comparison in a manuscript must read them first.
