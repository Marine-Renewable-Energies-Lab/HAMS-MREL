# Dissipative (epsilon) damping lid: design

Status: approved in conversation on 2026-11-06 (four sections). Branch: `feature/damping-lid`, on top of 2854242.
Replaces the phi = 0 lid of `2026-11-05-damping-lid-design.md`.

## Goal

A damping lid of any 2D planform on z = 0 that damps gap and moonpool resonance by absorbing energy, as the WAMIT damping
lid does (WAMIT v7.5 manual, Section 12.8 and eqs. 15.91 to 15.93, IDAMPER = 1). The first-order formulation of HAMS-MREL
is not changed: the hull boundary integral equation, the Green functions, the panel integrals, symmetry and the solvers stay
as they are. The lid only adds terms.

## Why the phi = 0 lid is replaced

Derived. phi = 0 makes the mean energy flux through the lid, proportional to Im(phi* dphi/dz), zero: the lid reflects
and absorbs nothing. It also lets water cross a lid that does not move (eta = 0 while dphi/dz is not zero). The WAMIT lid
puts a pressure proportional to the vertical velocity on the free surface, p/rho = epsilon dphi/dz, which absorbs energy.

## Formulation

Time factor e^{-i omega t}, nu = omega^2/g. WAMIT uses e^{+i omega t}; the equations below are the complex conjugates of
WAMIT's. The dynamic condition with p/rho = epsilon phi_z on the lid gives

    phi_z = beta phi,    beta = nu g / (g - i omega epsilon)

The wave Green function satisfies dG/dz = nu G on z = 0, so the free-surface integral of Green's identity leaves only the
lid, with q = phi_z - nu phi = alpha phi and

    alpha = beta - nu = i omega epsilon nu / (g - i omega epsilon)

Unknowns: phi on the hull panels and phi on the lid panels. D and S are the double- and single-layer panel integrals as now,
j runs over hull panels, l over lid panels, alpha_l uses the epsilon of panel l:

    hull row i:  2 pi phi_i + sum_j D_ij phi_j - sum_l alpha_l S_il phi_l = sum_j S_ij V_j
    lid row l:   4 pi phi_l + sum_j D_lj phi_j - sum_l' alpha_l' S_ll' phi_l' = sum_j S_lj V_j

The lid row has the 4 pi free term of WAMIT eq. 15.93. Interior rows of IRSP (field points inside the hull) have no free
term, as now.

Diffraction is solved for the scattered potential. The incident potential satisfies phi_I,z = nu phi_I on z = 0, so the
scattered potential has q_S = alpha (phi_S + phi_I) on the lid, and every row (hull, lid, interior) gains

    + sum_l alpha_l S_xl phi_I,l     on the right-hand side

The 4 pi phi_I lid term of the phi = 0 lid is removed.

Field points:

    4 pi phi(x) = sum_j (S_xj V_j - D_xj phi_j) + sum_l alpha_l S_xl (phi_l + phi_I,l)

with phi_I,l only for the diffraction problem.

Elevation on a lid panel l, from the total potential phi_T: eta = i omega phi_T / (g - i omega epsilon_l). Elsewhere it is
i omega phi_T / g as now. The kinematic condition gives the equivalent phi_z = -i omega eta.

Absorbed power, used by the tests only: P = 1/2 rho sum_l epsilon_l |phi_z,l|^2 A_l = 1/2 rho omega^2 sum_l epsilon_l
|eta_l|^2 A_l >= 0.

### Limits that follow

- epsilon = 0: alpha = 0. The lid rows only evaluate the representation and do not feed the hull rows, so the hull
  solution equals the one without a lid.
- epsilon to infinity: alpha to -nu, phi_z to 0 and eta to 0: a rigid lid.
- Infinite frequency: the lid is skipped, as now.
- Zero frequency: the run stops when the lid is on, as WAMIT does.
- A and B stay symmetric by reciprocity, since the lid condition is the same homogeneous condition in every problem. B is
  positive semi-definite.
- The Haskind relation over the hull no longer holds for epsilon > 0: the lid adds sum_l alpha_l phi_j phi_I A_l.

### Symmetry

The lid columns of symmetry class IP are, as now, sum_IS RXY(IS,IP) S(x reflected for IS, l) times alpha_l. The incident
term on the right uses the class components of phi_I on the lid in the scaling of the unknowns, U_I^IP(l) = sum_J
RXY(IP,J) phi_I(centroid l reflected for J), so the term is sum_l alpha_l LS^IP(x,l) U_I^IP(l), where LS^IP is the lid
column block before it is scaled by -alpha and factorised. LS is stored by the assembly for this purpose.

## Input

`Damping_lid 1` in `ControlFile.in` stays the switch. `DampingLid.pnl` keeps its header, nodes and panels and gains one
block after the panels, in the layout of the WAMIT `.dmp` file:

    # Start Definition of Damping
          2                       ! NDAMPERS
          1    40    0.20         ! IP1  IP2  EDAMPER (m/s)
         41    80    0.10
    # End Definition of Damping

- IP1 and IP2 are panel numbers of `DampingLid.pnl`, 1 to LNELEM. With ISX or ISY the images take the epsilon of their
  source panel.
- epsilon is the dimensional WAMIT coefficient in m/s, constant over frequency. There is no default.

Stops, each naming the line or panel: block missing; NDAMPERS < 1; IP1 > IP2 or outside 1 to LNELEM; overlapping ranges; a
lid panel in no range; epsilon < 0; an unreadable value. epsilon = 0 is allowed. A zero frequency in the frequency list
with the lid on stops the run. The existing stops stay (GMRES, IRSP with several bodies, lid inside a hull, file and mesh
checks).

The documentation gives the block, the WAMIT and Gourlay (2019) starting values (0.1 to 1.0 in the WAMIT tests, 0.2 for
side-by-side vessels) with the statement that epsilon must be calibrated against model tests or CFD, and Gourlay's
small-epsilon conversion from Chen's form, epsilon ~ g epsilon_Chen / omega, for information only.

## Outputs

No file format changes, with one exception: a single-body lid run also writes the radiation elevation rows of field
points at z = 0 to `PressureElevation.6p`, which a no-lid run does not write (no-lid output is unchanged). A, B,
excitation and RAOs are integrated over the hull only, as now and as in WAMIT. Field points
on z = 0 inside a lid panel use the lid elevation; others as now. Pressure in the fluid is unchanged in form. One log line
per damper at start: range, number of panels, epsilon. No absorbed-power output.

`MXQ` is renamed `MXLPOT` (lid potential), so the old meaning cannot be used by mistake.

## Code

- `DampingLid.f90`: read the damping block (`READ_LID_MESH`), store `LEPS(LNELEM)`; `LID_ASSEMBLE` stores the lid column
  block and fills -alpha_l times it, plus 4 pi on the lid diagonal; `LID_DIFFRACTION` adds the incident lid term instead of
  4 pi phi_I; a function for alpha; a function that returns the epsilon at a point on z = 0 (0 when on no lid panel),
  testing the symmetry images of the point.
- `DampingLidMulti.f90`: the same for several bodies.
- `PressureElevation.f90`, `PressureElevationMulti.f90`: the field lid term and the lid elevation.
- `HAMS_Prog.f90`: the zero-frequency stop, the damper log lines.
- `WavDynMods.f90` (`Lid_mod`): `MXLPOT`, `LEPS`.
- Not touched: CALGREEN, the hull D and S, the panel integrals, finite-depth handling, symmetry routines, ZGETRF, ZGETRS,
  ZGELS, `FILL_HULL_LEFT`.

## Testing

Fixtures stay at 640 lid panels or fewer and run on the laptop. A convergence study with finer lids may run on Snellius
(SSH key login only, job size and cost shown before submitting).

1. epsilon = 0 gives no lid: hull outputs of `lid_ann_160`, `lid_gap` and `lid_two_far` with epsilon = 0 equal the runs
   without a lid to 1e-12 relative.
2. Large epsilon gives a rigid lid: epsilon = 1e6 against the eigenfunction reference with phi_z = 0 on the lid, and
   |eta| on the lid centroids small against the incident amplitude.
3. Infinite frequency: identical to no lid.
4. `annular_lid_reference.py` takes epsilon: in the lid region the dispersion relation k tanh(kh) = nu g/(g - i omega
   epsilon), complex roots by Newton iteration from the open-water roots, matching as now.
5. `mm_cyl_lid` a, b, c (2, 4, 8 rings) for epsilon = 0.2, 1, 1e6: A, B, X and eta at the field points against the
   reference; error falling about as 1/N; target within 1 % at 8 rings, set from the first measurement.
6. Energy balance in the reference: 1/2 B_jj |U|^2 = far-field power + lid power, to 1e-6.
7. Lid power from HAMS-MREL via field points at the lid centroids, P >= 0, against the reference lid power within the
   tolerance of 5.
8. Reciprocity for `lid_gap` with epsilon = 0.2: A and B symmetric, B positive semi-definite.
9. `lid_ann_160_half` equals the full mesh.
10. `lid_gap`: the peak gap elevation falls for epsilon 0, 0.2, 1, near the epsilon = 0 resonance frequency.
11. One fixture per new input stop and for the zero frequency.
12. All runs without a lid bit-identical to before; lid fixtures get new references, recorded in
    `test/data/damping_lid/README.md`.

Not covered: no WAMIT or Gourlay geometry is available, so there is no cross-code comparison.
