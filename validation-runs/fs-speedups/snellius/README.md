# fs-speedups on Snellius: G1 item 3 and the OC4 comparison (cluster-side steps)

Nothing here has been submitted. Submit only after the cost is approved (COST.md: about 366 core-hours expected, 852 at
most). Login by key only; in Git Bash every ssh/scp line needs the `MSYS_NO_PATHCONV=1` prefix (paths start with /).

- Run root: `/home/vraghavan/HAMS_drift_study/run_fsspeed` (new).
- Repository: `/home/vraghavan/HAMS_drift_study/repo_fsspeed` (new), a clone of a git bundle of feature/qtf-fs-speedups.
  The job scripts are run from the clone (`repo_fsspeed/validation-runs/fs-speedups/snellius/...`), so they arrive with
  Unix line endings and need no scp of their own. build.sbatch stops unless the source equals that of 3af6cf8.
- Inputs of the OC4 run: those of the old OC4 run (code at 22086da), already on Snellius in `run_oc4wsum/q_oc4wfull_sym_indsum_free_L0/input`
  (checked against `oc4_compare/md5sums_input.txt` by the job). The outputs of the old run there are the reference.

| Script | Runs | Limit | CPUs | Memory | Depends on |
|---|---|---|---|---|---|
| build.sbatch | `run_fsspeed/hamsmrel_grn4` | 0:30 | 24 | 32G | - |
| g1_item3/run_sharp.sbatch | array 1-4 of sharp.txt: q_sharp_directsum_fixed_m2 and _free_m2, samples and panels | 6:00 each | 24 | 32G | build |
| g1_item3/run_round.sbatch | array 1-2 of round.txt: q_round_directsum_free_m2, samples and panels | 1:00 each | 24 | 32G | build |
| oc4_compare/run_oc4.sbatch | q_oc4wfull_sym_indsum_free_L0_fsspeed (the inputs of the old run + qtf_fs_circle = panels) | 9:00 | 24 | 32G | build |

The G1 cases are the 7c-2 cases of drift_convergence.py (qtf = 3, so sum and difference; qtf_direct = 1; qtf_terms = full;
IRSP = 1), written by g1_item3/make_case.py with 24 threads and `qtf_fs_circle = samples` or `= panels` appended; every
other key at the branch default (field_deep_kh = 10, qtf_fs_inner = stencil). Both routes run on the same build, so the
difference between them is the circle route: the panels route itself, plus the deep-water rows that the samples route uses
for its circle samples where k0 h >= 10 (about 1 rad/s and up in depth 100 m), while the panels route is always finite
depth (pw_from_panels). Against the exact references that is the comparison G1 wants: each route as it would run by
default. The OC4 run adds only `qtf_fs_circle = panels` to the inputs of the old OC4 run; the
other items of the branch are on by default (field_deep_kh = 10, psi reuse, the threaded mode check).

## 1. Bundle (laptop, Git Bash)

```
cd /c/Users/vraghavan1/Downloads/MREL_HAMS-MREL
mkdir -p validation-runs/fs-speedups/xfer
git -C HAMS-MREL-fs-speed bundle create ../validation-runs/fs-speedups/xfer/fsspeed.bundle feature/qtf-fs-speedups
md5sum validation-runs/fs-speedups/xfer/fsspeed.bundle
```

## 2. Folders, bundle and clone (Snellius)

```
MSYS_NO_PATHCONV=1 ssh -o BatchMode=yes -o PasswordAuthentication=no vraghavan@snellius.surf.nl \
    'mkdir -p /home/vraghavan/HAMS_drift_study/run_fsspeed/logs /home/vraghavan/HAMS_drift_study/xfer_fsspeed'
MSYS_NO_PATHCONV=1 scp validation-runs/fs-speedups/xfer/fsspeed.bundle \
    vraghavan@snellius.surf.nl:/home/vraghavan/HAMS_drift_study/xfer_fsspeed/
MSYS_NO_PATHCONV=1 ssh -o BatchMode=yes -o PasswordAuthentication=no vraghavan@snellius.surf.nl 'cd /home/vraghavan/HAMS_drift_study && \
    md5sum xfer_fsspeed/fsspeed.bundle && \
    git clone -b feature/qtf-fs-speedups xfer_fsspeed/fsspeed.bundle repo_fsspeed && git -C repo_fsspeed log -2 --oneline && \
    ls run_oc4wsum/q_oc4wfull_sym_indsum_free_L0/input run_oc4wsum/q_oc4wfull_sym_indsum_free_L0/result.json'
```

`git log -2` must show the docs commit of Task 7 on top of `3af6cf8 test: circle radius check ...`.

## 3. Submit (only after the approval)

```
cd /home/vraghavan/HAMS_drift_study/repo_fsspeed/validation-runs/fs-speedups/snellius
B=$(sbatch --parsable build.sbatch)
sbatch --dependency=afterok:$B g1_item3/run_sharp.sbatch
sbatch --dependency=afterok:$B g1_item3/run_round.sbatch
sbatch --dependency=afterok:$B oc4_compare/run_oc4.sbatch
```

Logs in `run_fsspeed/logs` (build_<job>.log, g1_sharp_<array>_<task>.log, g1_round_..., oc4_<job>.log). Each run writes
`log.txt`, `time.txt`, `result.json` and `output/` in its case folder: `run_fsspeed/g1_<route>/<case>/` and
`run_fsspeed/q_oc4wfull_sym_indsum_free_L0_fsspeed/`. What to look for early in a log.txt: no line "qtf_fs_circle and
qtf_fs_inner are ignored", and no "the partial waves miss the field on the circle" warning (the jobs grep for both at the
end). Cost check afterwards: `sacct -j <ids> -X --format=JobID,JobName,Elapsed,NCPUS,CPUTimeRAW,MaxRSS,State`.

## 4. Afterwards (laptop)

Copy back into `validation-runs/fs-speedups/snellius/results/` (outside git: the outputs are large):

```
MSYS_NO_PATHCONV=1 ssh -o BatchMode=yes -o PasswordAuthentication=no vraghavan@snellius.surf.nl 'cd /home/vraghavan/HAMS_drift_study/run_fsspeed && \
    tar -czf fsspeed_back.tgz g1_samples/*/log.txt g1_samples/*/result.json g1_samples/*/time.txt g1_samples/*/output/Wamit_format \
    g1_panels/*/log.txt g1_panels/*/result.json g1_panels/*/time.txt g1_panels/*/output/Wamit_format \
    q_oc4wfull_sym_indsum_free_L0_fsspeed/log.txt q_oc4wfull_sym_indsum_free_L0_fsspeed/result.json \
    q_oc4wfull_sym_indsum_free_L0_fsspeed/time.txt q_oc4wfull_sym_indsum_free_L0_fsspeed/output/Wamit_format \
    q_oc4wfull_sym_indsum_free_L0_fsspeed/compare_oc4.txt'
MSYS_NO_PATHCONV=1 scp vraghavan@snellius.surf.nl:/home/vraghavan/HAMS_drift_study/run_fsspeed/fsspeed_back.tgz \
    /c/Users/vraghavan1/Downloads/MREL_HAMS-MREL/validation-runs/fs-speedups/xfer/
cd /c/Users/vraghavan1/Downloads/MREL_HAMS-MREL
mkdir -p validation-runs/fs-speedups/results && tar -xzf validation-runs/fs-speedups/xfer/fsspeed_back.tgz -C validation-runs/fs-speedups/results
V=validation-runs; S=HAMS-MREL-fs-speed/validation-runs/fs-speedups/snellius
python -I $S/g1_item3/g1_compare.py $V/fs-speedups/results HAMS-MREL-fs-speed \
    --sum-ref $V/snellius-qtf-directsum/qtfdirectsum_reference.csv --diff-ref $V/snellius-qtf-fs/qtffs_reference.csv
python -I $S/oc4_compare/compare_oc4.py $V/fs-speedups/results/q_oc4wfull_sym_indsum_free_L0_fsspeed \
    $V/oc4wfull-sum/results/q_oc4wfull_sym_indsum_free_L0
```

What the exact references hold (checked 2026-11-20 on the CSV headers and read_potential_reference of qtf_compare.py):
both `qtfdirectsum_reference.csv` (stages 1-2 and 5; kernel sum, fixed and free, 0.40-1.40 by 0.05, every pair
omega_i <= omega_j) and `qtffs_reference.csv` (the 6a reference; kernel diff, fixed and free, 0.40-1.40 by 0.01,
omega_j - omega_i <= 0.3) have the quadratic part, the five body-forcing terms and the free-surface part FS. So each FS part
of HAMS is compared with FS, each total without its FS part with the total without FS, and each total with FS with the
total with FS: like with like. There is no exact reference for the round hull (in-run check only).

g1_compare.py writes `g1_compare.txt` and `g1_compare.png` in the results folder (qtf_compare.py is imported for its
readers, not run: its --direct mode compares the .12d/.11d of these qtf_terms = full cases, FS included, with the
reference without FS):
1. sharp m2, both kernels, both routes against the exact cylinder: the FS parts (difference: FS_inner + FS_outer of
   QTF_PotentialTerms.txt and D_FS of QTF_DirectTerms.txt; sum: P_FS and D_FS of QTF_DirectTermsSum.txt), the totals minus
   their FS part, the body parts, and the totals with FS; free body also with 0.75-0.80 rad/s left out. Then the
   panels-minus-samples change and the G1 verdict per quantity on the max and the median;
2. the in-run direct-indirect check of both kernels (round hull included);
3. every output file of panels against samples; 4. the core-hours.
It exits 1 and lists PROBLEMS when a case, route, output file, row set or result.json is missing (cases from sharp.txt
and round.txt). compare_oc4.py writes `compare_oc4.txt` and `compare_oc4.png` beside the OC4 run: every QTF file against
the old run (max and median), the wall time and core-hours of both, and the logged passes; it too exits 1 on any problem.

Today's route at 47d5911 (stage 5 outputs, results_wl), through g1_compare.py with the stage-5 run as both routes:
FS parts against the exact FS, max / median %, difference kernel FS_inner + FS_outer: fixed 1.82 / 0.03 (F_x), 1.91 / 0.01
(F_z), 1.41 / 0.01 (M_y); free, 0.75-0.80 left out, 2.34 / 0.00, 1.21 / 0.00, 1.97 / 0.00. Sum kernel D_FS: fixed 1.74 /
0.28, 1.05 / 0.31, 1.50 / 0.36; free, 0.75-0.80 left out, 2.03 / 0.52, 0.99 / 0.21, 1.83 / 0.46 (as compare_stage5*.txt).
The difference totals stay far off (fixed F_z 15.8 / 11.4 %) because of the quadratic part .10d (20.8 / 15.0 % on its own),
which is the same in both routes.

G1 rule (spec G1): item 3 becomes the default only if no quantity of the panels runs is worse than the samples runs by
more than the larger of 0.5 percentage points and 10 % of the samples error, against the exact reference. Read the FS rows
and the change column first: the totals carry the route-independent quadratic error, whose large values widen the 10 %
margin.

Known limits of the compare scripts: `qtf_files_compare.py` reports a file that exists in one folder only as a problem (exit 1), so a legitimate extra file in a run would fail the check although the outputs are written; `g1_compare.py` indexes the terms files by column and raises IndexError (non-zero exit) on a terms file from a run without `qtf_terms = full`; all G1 cases use `qtf_terms = full`.
