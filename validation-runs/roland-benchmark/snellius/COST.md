# VolturnUS-S run: cost (nothing submitted)

fat_genoa, account <account>, 24 CPUs (the billed minimum), 32G, one job. No build job: the run uses the existing binary
`run_fsspeed/hamsmrel_grn4` (fs-speedups branch, source of 3af6cf8) with `qtf_fs_circle = panels`, which equals the
branch tip 80eaf3d (the later source changes are log lines and that default).

| Option | Job | Expected wall | Expected core-h | Limit | Max core-h |
|---|---|---|---|---|---|
| A (recommended) | indirect and direct, free surface on (`qtf_direct = 1`) | 5.1 h | 122 | 18 h | 432 |
| B | indirect only, free surface on | 3.5 h | 84 | 14 h | 336 |

About 1.5 SBU per core-hour on fat_genoa (the fs-speedups estimate: 366 core-h about 550 SBU), so A is about 180 SBU
expected and 650 at most.

## Basis of the estimate

The case is close in size to the OC4 WAMIT model: 5076 hull panels against 5022, the same 34 frequencies as the OC4
difference runs, a 4412-triangle inner free surface (OC4: about 4044 inner points).

| Part | Basis | Expected | Upper bound used for the limit |
|---|---|---|---|
| First order and the rest | OC4 sum run with the new defaults: 0.48 h for 21 frequencies; times 34/21 | 0.8 h | 0.8 h |
| Body part of the difference QTF | OC4 difference run of the full hull, 5022 panels, older code: 5.08 h in all, 4.3 h of it second order; the body part of the sum run fell 2.94 to 0.76 h with the new code | 1.2 h | 4.3 h (no saving) |
| Free-surface part | OC4 sum run, new defaults: 0.83 h; times 34/21 frequencies and 1.1 for the inner mesh | 1.5 h | 7.7 h (old-code rate) |
| Direct method (option A only) | 33 difference frequencies, one assembly, factorisation and free-surface forcing each, about 3 min each | 1.6 h | 3 h |

The limit is the upper-bound sum (15.8 h for A, 12.8 h for B) rounded up. A wrong input stops the run in minutes (the log
check at the end of the job), so a failed start costs little.

## Why option A

WAMIT's second-order run (the reference) used its direct method (".12d ... DIRECT METHOD" in its log). Option A gives the
HAMS direct result (`QTF_Direct.12d`) as well as the indirect one (`QTF_Indirect.11d`), so both are compared with WAMIT,
like with like, and the two HAMS methods against each other on a new hull.

## Extension to the full range of the paper (prepared, not submitted)

Roland et al. (2026) run 0.30 to 2.50 rad/s for both platforms (their section 3.3; Figs. 5 and 7 reach 2.50, the OC4 WAMIT
file 3.0). HAMS stops at 1.90 because `velocity_method = source` cannot remove irregular frequencies (IRSP = 1 is rejected).
The first ones are about 1.98 rad/s for the OC4 offset columns (radius 6 m) and 1.94 rad/s for the VolturnUS-S outer
columns (6.25 m), with more below 2.55 rad/s: J_m(ka) = 0 gives 2.50 (OC4, j_11) and 2.17, 2.45 (VolturnUS-S, central
column j_01 and outer j_11) [derived: k = j_mn / a, omega^2 = g k for these drafts]. So the extension needs
`velocity_method = differencing` with IRSP = 1 and a lid, and a full rerun: a QTF pair needs the first-order states of both
its frequencies in one run, so the missing frequencies cannot be added to the existing runs. The overlap 0.25 to 1.90 rad/s
then also compares differencing with the source velocity of the existing runs.

Wavelength at 2.55 rad/s: 9.5 m (deep water). Longest panel edges: OC4 9208 panels median 0.44 m, 90 % 1.87 m (lambda/5),
largest in the top 3 m 2.45 m (lambda/3.9); OC4 5022 panels 0.69 / 2.24 / 2.29 m; VolturnUS-S 1.62 / 1.62 m, 1.40 m at the
waterline (lambda/7 at 2.50). The OC4 meshes resolve 2.55 rad/s with most panels but not with the largest near-surface ones;
the VolturnUS-S mesh does. The VolturnUS-S inner free surface of 2 m is lambda/4.9 at 2.50 rad/s; 1.25 m (about 11,000
triangles) keeps lambda/8.

| Job | Frequencies | Expected | Limit (max billed) |
|---|---|---|---|
| OC4 full hull, 9208 panels (recommended) | 47, 0.25 to 2.55 | 8.4 h, 202 core-h | 22 h, 528 core-h |
| OC4 full hull, 5022 panels | 47, 0.25 to 2.55 | 3.1 h, 74 core-h | 9 h, 216 core-h |
| VolturnUS-S option A, inner FS 1.25 m (recommended) | 46, 0.25 to 2.50 | 9.9 h, 238 core-h | 30 h, 720 core-h |
| VolturnUS-S option A, inner FS 2 m | 46, 0.25 to 2.50 | 6.7 h, 161 core-h | 18 h, 432 core-h |
| VolturnUS-S option A as approved | 34, 0.25 to 1.90 | 5.1 h, 122 core-h | 18 h, 432 core-h |

Basis. OC4: the logs of the existing runs. The control-surface field passes took 4.03 h of 5.08 h (5022 panels) and 9.27 h
of 12.72 h (9208 panels), with 67 passes for 34 frequencies and 4 radiation-only solves; the source solve 0.47 and 1.48 h.
With 47 frequencies the passes are 98 (times 1.46). The fs-speedups code cut these passes to 0.36 of their time on OC4
(4538 s to 1613 s, deep-water rows above k0 h = 10); the upper bound takes no saving. The rest scales with the frequencies
(times 1.38) and the lid (times 1.3, assumed); the source solve drops. VolturnUS-S: the parts of option A times 46/34 for
the passes and the free-surface part, times 45/33 for the direct solves, times 2.6 for the finer inner mesh, and the lid for
the first order. The 9208-panel OC4 run is enough for the comparison: the 5022-panel run adds only the mesh trend, which the
existing runs already show up to 1.90 rad/s.

Chosen jobs (cases written by `make_oc4_extension.py` and `make_volturnus_case.py --fmax 2.50 --fs-panel 1.25 --lid`;
job files `run_oc4_ext.sbatch` and `run_volturnus_ext.sbatch`): OC4 9208 panels, expected 202 core-h, limit 528;
VolturnUS-S 0.25-2.50 with the 1.25 m inner mesh, expected 238, limit 720. Total expected about 440 core-h (about 660 SBU),
at most 1248 core-h. The lids: OC4 1674 panels on the 7 waterline loops (4 columns, 3 brace crossings), nearest centroid
0.44 waterline row heights from the wall; VolturnUS-S 654 panels on the 4 columns, 0.50 row heights.
