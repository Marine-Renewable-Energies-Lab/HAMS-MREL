# VolturnUS-S on Snellius (cluster-side steps)

Nothing here has been submitted. Submit only after the user's approval of the cost (COST.md). Login by key only; in Git
Bash every ssh/scp line needs the `MSYS_NO_PATHCONV=1` prefix.

- Run root: `/home/vraghavan/HAMS_drift_study/run_roland` (new).
- Binary: `/home/vraghavan/HAMS_drift_study/run_fsspeed/hamsmrel_grn4` (existing; see run_volturnus.sbatch).
- Inputs: written on the laptop by `make_volturnus_case.py` into
  `C:/Users/vraghavan1/Downloads/MREL_HAMS-MREL/validation-runs/roland-benchmark/cases/q_volturnus_free_full/input`
  (outside git: derived from the IEA mesh) and checked on the cluster against `md5sums_input.txt`.

## 1. Copy (laptop, Git Bash)

A Windows checkout may give the job file and the checksum list CRLF line ends; the last line strips them on the cluster.

```
C=/c/Users/vraghavan1/Downloads/MREL_HAMS-MREL/validation-runs/roland-benchmark/cases/q_volturnus_free_full
J=validation-runs/roland-benchmark/snellius
MSYS_NO_PATHCONV=1 ssh -o BatchMode=yes -o PasswordAuthentication=no vraghavan@snellius.surf.nl \
    'mkdir -p /home/vraghavan/HAMS_drift_study/run_roland/logs /home/vraghavan/HAMS_drift_study/run_roland/q_volturnus_free_full/input'
MSYS_NO_PATHCONV=1 scp $C/input/* vraghavan@snellius.surf.nl:/home/vraghavan/HAMS_drift_study/run_roland/q_volturnus_free_full/input/
MSYS_NO_PATHCONV=1 scp $J/md5sums_input.txt $J/run_volturnus.sbatch vraghavan@snellius.surf.nl:/home/vraghavan/HAMS_drift_study/run_roland/
MSYS_NO_PATHCONV=1 ssh -o BatchMode=yes -o PasswordAuthentication=no vraghavan@snellius.surf.nl \
    "cd /home/vraghavan/HAMS_drift_study/run_roland && sed -i 's/\r\$//' md5sums_input.txt run_volturnus.sbatch"
```

For option B of COST.md (indirect only), delete the line `qtf_direct = 1` from `SecondOrder.in` before the copy, rewrite
`md5sums_input.txt`, and set `#SBATCH -t 14:00:00`.

## 2. Submit (only after the approval)

```
cd /home/vraghavan/HAMS_drift_study/run_roland && sbatch run_volturnus.sbatch
```

Log in `run_roland/logs/volturnus_<job>.log`; the case writes `log.txt`, `time.txt`, `result.json` and `output/`. Look
early in `log.txt` for the free-surface lines (inner points, circle R, modes) and for any stop message. Cost afterwards:
`sacct -j <id> -X --format=JobID,JobName,Elapsed,NCPUS,CPUTimeRAW,MaxRSS,State`.

## 3. Back (laptop)

```
MSYS_NO_PATHCONV=1 ssh -o BatchMode=yes -o PasswordAuthentication=no vraghavan@snellius.surf.nl 'cd /home/vraghavan/HAMS_drift_study/run_roland && \
    tar -czf volturnus_back.tgz q_volturnus_free_full/log.txt q_volturnus_free_full/time.txt q_volturnus_free_full/result.json \
    q_volturnus_free_full/output/Wamit_format'
MSYS_NO_PATHCONV=1 scp vraghavan@snellius.surf.nl:/home/vraghavan/HAMS_drift_study/run_roland/volturnus_back.tgz \
    /c/Users/vraghavan1/Downloads/MREL_HAMS-MREL/validation-runs/roland-benchmark/results/
```

Then `compare_diagonals.py` (figure 7) and `compare_maps.py` against `IEA-15-240-RWT-UMaineSemi.12d` (see ../README.md).

## Full-range jobs (not submitted; cost in COST.md)

The same steps for `q_oc4wfull_sym_free_L1_ext` (job file `run_oc4_ext.sbatch`, inputs from `make_oc4_extension.py`)
and `q_volturnus_free_full_ext` (`run_volturnus_ext.sbatch`, inputs from `make_volturnus_case.py --fmax 2.50 --fs-panel
1.25 --lid`): copy each case's `input/` to `run_roland/<case>/input/`, copy `md5sums_<case>.txt` and the job file to
`run_roland/`, strip CR as above, then `sbatch` the job file. Both reuse the binary of section 2.

Each case runs in its own folder, `run_roland/<case>/` (log.txt, time.txt, result.json, output/), with its Slurm log in
`run_roland/logs/` (`oc4_ext_<job>.log`, `volturnus_ext_<job>.log`). fat_genoa allows up to 5 days, so the limits of
22 h and 30 h fit. The commands (laptop, Git Bash, from the repository root):

```
C=/c/Users/vraghavan1/Downloads/MREL_HAMS-MREL/validation-runs/roland-benchmark/cases
J=validation-runs/roland-benchmark/snellius
R=/home/vraghavan/HAMS_drift_study/run_roland
for c in q_oc4wfull_sym_free_L1_ext q_volturnus_free_full_ext; do
  MSYS_NO_PATHCONV=1 ssh -o BatchMode=yes -o PasswordAuthentication=no vraghavan@snellius.surf.nl "mkdir -p $R/logs $R/$c/input"
  MSYS_NO_PATHCONV=1 scp $C/$c/input/* vraghavan@snellius.surf.nl:$R/$c/input/
  MSYS_NO_PATHCONV=1 scp $J/md5sums_$c.txt vraghavan@snellius.surf.nl:$R/
done
MSYS_NO_PATHCONV=1 scp $J/run_oc4_ext.sbatch $J/run_volturnus_ext.sbatch vraghavan@snellius.surf.nl:$R/
MSYS_NO_PATHCONV=1 ssh -o BatchMode=yes -o PasswordAuthentication=no vraghavan@snellius.surf.nl \
    "cd $R && sed -i 's/\r\$//' md5sums_q_*_ext.txt run_*_ext.sbatch && \
     (cd q_oc4wfull_sym_free_L1_ext/input && md5sum -c ../../md5sums_q_oc4wfull_sym_free_L1_ext.txt) && \
     (cd q_volturnus_free_full_ext/input && md5sum -c ../../md5sums_q_volturnus_free_full_ext.txt)"
```

Submit (only after the approval): `cd /home/vraghavan/HAMS_drift_study/run_roland && sbatch run_oc4_ext.sbatch && sbatch
run_volturnus_ext.sbatch`.
