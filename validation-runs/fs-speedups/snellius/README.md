# fs-speedups on Snellius: G1 item 3 and the OC4 comparison (cluster-side steps)

Nothing here has been submitted. Submit only after the cost is approved (COST.md: about 363 core-hours expected, 852 at
most). Login by key only; in Git Bash every ssh/scp line needs the `MSYS_NO_PATHCONV=1` prefix (paths start with /).

- Run root: `/home/vraghavan/HAMS_drift_study/run_fsspeed` (new).
- Repository: `/home/vraghavan/HAMS_drift_study/repo_fsspeed` (new), a clone of a git bundle of feature/qtf-fs-speedups.
  The job scripts are run from the clone (`repo_fsspeed/validation-runs/fs-speedups/snellius/...`), so they arrive with
  Unix line endings and need no scp of their own. build.sbatch stops unless the source equals that of 3af6cf8.
- Inputs of the OC4 run: those of 29453291, already on Snellius in `run_oc4wsum/q_oc4wfull_sym_indsum_free_L0/input`
  (checked against `oc4_compare/md5sums_input.txt` by the job). The 29453291 outputs there are the reference.

| Script | Runs | Limit | CPUs | Memory | Depends on |
|---|---|---|---|---|---|
| build.sbatch | `run_fsspeed/hamsmrel_grn4` | 0:30 | 24 | 32G | - |
| g1_item3/run_sharp.sbatch | array 1-4 of sharp.txt: q_sharp_directsum_fixed_m2 and _free_m2, samples and panels | 6:00 each | 24 | 32G | build |
| g1_item3/run_round.sbatch | array 1-2 of round.txt: q_round_directsum_free_m2, samples and panels | 1:00 each | 24 | 32G | build |
| oc4_compare/run_oc4.sbatch | q_oc4wfull_sym_indsum_free_L0_fsspeed (29453291 + qtf_fs_circle = panels) | 9:00 | 24 | 32G | build |

The G1 cases are the 7c-2 cases of drift_convergence.py (qtf = 3, so sum and difference; qtf_direct = 1; qtf_terms = full;
IRSP = 1), written by g1_item3/make_case.py with 24 threads and `qtf_fs_circle = samples` or `= panels` appended; every
other key at the branch default (field_deep_kh = 10, qtf_fs_inner = stencil). Both routes run on the same build, so the
difference between them is item 3 alone. The OC4 run adds only `qtf_fs_circle = panels` to the inputs of 29453291; the
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

g1_compare.py writes `g1_compare.txt` and `g1_compare.png` in the results folder: the sum QTF (.12s, .12s - .10s, P, D_FS)
and the difference QTF (.12d, .11d) of both routes against the exact sharp cylinder (qtf_compare.py --direct-sum and
--direct), the in-run direct-indirect check (round hull included), every output file of panels against samples, and the
core-hours. compare_oc4.py writes `compare_oc4.txt` and `compare_oc4.png` beside the OC4 run: every QTF file against
29453291 (max and median), the wall time and core-hours of both, and the logged passes. The stage-5 numbers of today's
route at 47d5911 are in `validation-runs/snellius-qtf-directsum/compare_stage5*.txt`.

G1 rule (spec G1): item 3 becomes the default only if no quantity of the panels runs is worse than the samples runs by
more than the larger of 0.5 percentage points and 10 % of the samples error, against the exact reference.
