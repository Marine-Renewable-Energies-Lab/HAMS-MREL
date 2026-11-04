# Mean drift integration fixtures

Used by `test/test_integration_drift.f90`. Generated once with `pre_postprocessing_python/SecondOrder/make_mesh.py` and committed, so the tests do not call Python.

| Case | Purpose | Mesh | Frequencies (rad/s) | Headings | `SecondOrder.in` |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `drift_cyl` | Physics checks | see below | 0.70, 0.95, 1.10 | 0, 30 | `mean_drift=1`, free, `write_terms=2` |
| `drift_cyl_sybo` | `SYBO=1`, heading pairs | as `drift_cyl` | same, plus the limits | 0, 30 | `mean_drift=2`, free |
| `drift_cyl_gmres` | `ISOLV=2` | as `drift_cyl` | same | 0, 30 | `mean_drift=1`, free |
| `drift_offset`, `drift_offset_refl2` | `REFL` scaling on non-zero M_z | a=5, d=5, 96 panels, centred at (6, 8) | 1.2 | 30 | `mean_drift=1`, fixed |
| `drift_two_cyl` | Multi-body path, reusable by later tasks | two of a=5, d=5, 180 panels each, 30 m apart in y | 1.35, 1.55 | 0, 30 | `mean_drift=2`, free, `write_terms=2` |

## drift_cyl

    python make_mesh.py cylinder --radius 10 --draft 10 --ntheta 32 --nz 8 --nr 6 --graded --out HullMesh.pnl

Cylinder a = 10 m, draft 10 m, depth 100 m, 448 panels, no external loads, rotation centre (0, 0, 0). `field_panels.txt` lists the hull panel whose centroid is each field point in `ControlFile.in`.

Hydrostatics, from the 32-sided waterplane polygon, for a uniform solid: V = 3121.445 m3, Awp = 312.1445 m2, zG = zB = -5 m, M = 3.121445e6 kg, I44 = I55 = 1.815845e8 kg m2, I66 = 1.550726e8 kg m2, C33 = 3.061092e6 N/m, C44 = C55 = 7.603715e7 N m, and M15 = -M24 = M zG.

Natural frequencies, from a first-order sweep with the frequency-dependent added mass: heave 0.80 rad/s, pitch coupled with surge 0.76 rad/s. The frequencies straddle them without being at them. In a 0-degree wave, the motion per unit amplitude is about 1.8 m in heave and 0.2 rad in pitch at 0.70 rad/s, and 0.5 m in heave at 0.95 rad/s.

The irregular frequencies of the interior problem start at about 1.55 rad/s (J0(ka) = 0), above the frequencies used.

The accuracy of the energy and Haskind checks is limited by the panel method and improves with the mesh, at about first order in panel size. On this mesh the 1% limit does not hold beyond the frequencies used: at 0.65 rad/s the energy balance is off by 1.3% (the free body follows the long wave and radiates almost nothing), and at 1.20 and 1.30 rad/s the Haskind heave force is off by 1.3% and 2.2% (the heave force has decayed to about a tenth of the surge force).
