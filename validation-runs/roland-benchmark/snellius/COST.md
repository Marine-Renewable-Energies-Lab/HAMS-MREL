# VolturnUS-S run: cost (nothing submitted)

fat_genoa, account tdsei3290, 24 CPUs (the billed minimum), 32G, one job. No build job: the run uses the existing binary
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
