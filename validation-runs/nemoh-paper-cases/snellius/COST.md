# Nemoh-paper cases: cost (nothing submitted)

fat_genoa, 24 CPUs billed per job (the minimum), about 1.5 SBU per core-hour. No build job: the jobs use a copy of the
binary built from `development_2` at 8f27716 (the branch these cases were written from). Estimates from
`cost_estimate.py` (its docstring gives the basis: the measured log of the 9208-panel OC4 run of the Roland benchmark).

| Job | Content | Expected wall (h) | Core-h | SBU | Limit (h) | Core-h at limit | SBU at limit |
|---|---|---:|---:|---:|---:|---:|---:|
| J1 small | heave-plate cylinder m1 and m2, OC4 at the towing draught, cylinder 8b, SOFTWIND first order; Nemoh on the heave-plate m1 mesh beside | 1.9 | 45 | 68 | 5 | 120 | 180 |
| J2 sw_diff | SOFTWIND QTF⁻, 50 frequencies to 3.14 rad/s, no free-surface integral | 5.2 | 124 | 186 | 13 | 312 | 468 |
| J3 sw_full | SOFTWIND QTF⁻ and QTF⁺ with the free-surface integral, 32 frequencies to 2.01 rad/s | 10.0 | 241 | 361 | 26 | 624 | 936 |
| J4 oc4_bi | OC4 QTF⁻, headings 0 and 30, every heading pair, 50 frequencies to 3.14 rad/s | 0.6 | 14 | 20 | 3 | 72 | 108 |
| J5 oc3_diff | OC3 QTF⁻, 50 frequencies to 3.14 rad/s (Nemoh hydrosCal first) | 8.8 | 210 | 315 | 22 | 528 | 792 |
| J5b (option) | OC3 as J5, 32 frequencies to 2.01 rad/s | 5.7 | 136 | 204 | 15 | 360 | 540 |
| J2b (option) | SOFTWIND QTF⁻ as J2, 32 frequencies to 2.01 rad/s | 3.4 | 81 | 121 | 9 | 216 | 324 |
| **J1 to J5** | | | **633** | **950** | | **1656** | **2484** |
| J1, J2b, J3, J4, J5b | | | 516 | 774 | | 1392 | 2088 |

## Why J2, J3 and J5 cost most

The B2 term on the control surface C and the free surface F (`qtf_b2 = control`) takes one field pass per frequency, of
about N_hull × N_points Green-function rows. SOFTWIND (draft 91.4 m in 200 m) and OC3 (120 m in 320 m) are deep-draft
hulls: depth / draft is 2.2 and 2.7, below `field_deep_draft = 10`, so every pass keeps the finite-depth Green function
(65 times slower per row than the deep-water rows on the OC4 log). That is also the check of the guard asked for in
case D: the log line "the finite-depth Green function is kept" should appear. The generated C of these spars is large
(7328 and 9101 points, the default margins of 4 median panels). OC4 (depth / draft 10) gets the deep-water rows from
0.70 rad/s and is cheap.

The options J2b and J5b stop at 2.01 rad/s instead of 3.14 rad/s. The paper plots the QTFs to π rad/s; most of the
energy of the slow-drift excitation of these spars lies below 2 rad/s.

## Grids

The step is 0.0628 rad/s (2π/100, as in the HydroStar files of the Nemoh test cases), so every HAMS pair is a HydroStar
pair. J3 uses the grid of HydroStar's sum-frequency file (`SPG*.rao`: ω to 2.0096 rad/s); the largest difference there,
1.95 rad/s, stays below the pairs where `.11d` FS_outer diverged in the Roland full-range run (one frequency ≤ 0.30 and
the other ≥ 2.25 rad/s).
