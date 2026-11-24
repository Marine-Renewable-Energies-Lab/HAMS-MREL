# J14 and J15: OC3 low-row checks, cost (nothing submitted)

Two checks of the OC3 Hywind difference-QTF rows ω₁ ≤ 0.25 rad/s, where HAMS-MREL differs from HydroStar by up to 65 % of
max|QTF|. The rows above that agree to a median of 0.25–1.7 %. Both checks use 20 frequencies from 0.0628 to 1.256 rad/s
(step 0.0628, the grid of HydroStar and J5C), so they cover the rows ω₁ ≤ 0.25 against every ω₂ up to 1.256 rad/s.
Inputs come from `oc3_edge_study.py`.

- **J14, Nemoh v3.0.3 on its own test case 11.**
  - Mesh: 4842 panels with Nemoh's lid. Terms: DUOK + HASBO, no free-surface integral, as shipped.
  - Only the frequency list is changed. Every difference is on the list, as HASBO needs.
  - One thread: the build has no OpenMP.
  - `results/QTF` keeps the DUOK and HASBO parts apart, so the run also splits Nemoh's QTF into its quadratic and potential parts.
- **J15, HAMS-MREL on three structured meshes of the hull.**
  - Each mesh is refined only next to the bottom edge. Wall rows and bottom rings are graded by 1.25 per row from the
    edge size up to the uniform size.
    - e0: uniform (1.5 m rows, 5138 panels).
    - e1: edge rows 0.375 m (5306 panels).
    - e2: edge rows 0.094 m (5866 panels).
  - Settings, lid and Hydrostatic.in are those of J5C.
  - The control surface is held at J5C's (`cs_margin`, `cs_bottom_margin` 4 m, `cs_panel_size` 1 m). The job checks
    that C has J5C's 7370 panels.
  - Binary: `hamsmrel_50fc55e_J5C`, md5 checked.

fat_genoa bills at least 24 CPUs per job, at 1.5 SBU per core-hour.

| Job | Expected wall (h) | Core-h | SBU | Limit (h) | SBU at limit |
|---|---:|---:|---:|---:|---:|
| J14 Nemoh 3.0.3, OC3 test case (4842 panels), 1 thread | 4.0 (range 2.6–5.1) | 96 | 144 | 12 | 432 |
| J15 oc3_e0, uniform (5138 hull panels) | 6.0 | 143 | 214 | 9 | 324 |
| J15 oc3_e1, edge row 0.375 m (5306) | 6.2 | 149 | 223 | 10 | 360 |
| J15 oc3_e2, edge row 0.094 m (5866) | 7.0 | 168 | 252 | 11 | 396 |
| **All four** | | 556 | **834** | | **1512** |

How the times are estimated:
- **J15.** From J5C's stage times on the same binary.
  - The B2 field passes of the first 20 frequencies took 16 640 s. They are scaled by the number of hull panels.
  - The rest of the run took 127 s per frequency. It is scaled by (panels)^2.5.
  - The low frequencies carry most of the field time: 1100–1530 s each up to 0.25 rad/s. Cutting the list to 14
    frequencies (ω₂ ≤ 0.88) would save only about 20 %.
- **J14.** From the Nemoh runs J8 and J9 (deep water, LU solver): about 73 s per frequency at 4700 panels.
  - A finite-depth factor of 3–6 is assumed; it has not been measured.
  - QTFpreProc is taken to cost as much again.
  - This is the least certain number in the table, which is why its limit is 12 h.
- **Saving.** Nemoh uses one of the 24 billed CPUs. Run inside the e0 job (HAMS on 23 threads), J14 would add almost
  nothing to the bill, about 144 SBU less, but it would tie the two runs together.

Submit, from `<run_root>/run_nemohpaper`:

```
sbatch run_J14_nemoh_oc3.sbatch
sbatch -J np_J15_oc3_e0 -t 09:00:00 --export=ALL,CASE=oc3_e0 run_J15_oc3_edge.sbatch
sbatch -J np_J15_oc3_e1 -t 10:00:00 --export=ALL,CASE=oc3_e1 run_J15_oc3_edge.sbatch
sbatch -J np_J15_oc3_e2 -t 11:00:00 --export=ALL,CASE=oc3_e2 run_J15_oc3_edge.sbatch
```

What would settle the question:
- **J14 (Nemoh).**
  - If Nemoh's rows ω₁ ≤ 0.25 agree with HAMS-MREL to within a few % of max|QTF|, while both differ from HydroStar,
    then the difference is on HydroStar's side or is a convention we do not share.
  - If Nemoh agrees with HydroStar, the error is in HAMS-MREL. The DUOK/HASBO split then tells whether it sits in the
    quadratic part (QTF_Quadratic.10d) or the potential part (QTF_PotentialTerms.txt).
  - As a control, the rows ω₁ ≥ 0.31 must agree between the three codes as they do now.
- **J15 (edge refinement).**
  - The edge singularity is the cause if the low rows move toward HydroStar from e0 to e2 and settle (e1 to e2
    smaller than e0 to e1), with term III and the heave mean drift at 0.188 rad/s moving most.
  - The edge is ruled out if e0, e1 and e2 differ by less than about 2 % of max|QTF| in those rows while the gap to
    HydroStar stays at 20–65 %.
  - In both cases e0 must match J5C on the rows ω₁ ≥ 0.31, which checks that the structured mesh is equivalent to the
    Nemoh mesh.
