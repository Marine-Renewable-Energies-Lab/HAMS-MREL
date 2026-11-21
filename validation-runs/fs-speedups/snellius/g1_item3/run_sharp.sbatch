#!/bin/bash
#SBATCH -J fsspeed_g1_sharp
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=32G
#SBATCH -t 4:00:00
#SBATCH -o <run root>/run_fsfix/logs/g1_sharp_%A_%a.log
#SBATCH --array=2,4
# Re-run with the fix (C/F on finite-depth rows): the panels lines of sharp.txt only (2 and 4); expected 2.21 h each
# (6667 s of the first panels runs plus 1291 s of C/F passes on finite-depth rows), limit 4 h.
#SBATCH -A <account>
# G1 item 3 (qtf_fs_circle = panels against samples), sharp m2 cases of Phase 7c-2 (qtf = 3: sum and difference,
# qtf_direct = 1, qtf_terms = full, IRSP = 1), one task per line of sharp.txt (route, case). Expected 3.2 h each: in stage 5
# (stage 5, 47d5911) HAMS took 11436 s fixed and 11394 s free (the "code 0, ... s" lines of drift_convergence.py in
# snellius-qtf-directsum/logs/29405797_1.log and _2.log). Limit 6:00:00. Same modules, threads and stack as
# run_qtfdirectsum/jobs/stage5_wlfit.sbatch.
module load 2024
module load intel/2024a
ulimit -s unlimited
export OMP_NUM_THREADS=24
export OMP_STACKSIZE=1G
JOBS="<run root>/repo_fsfix/validation-runs/fs-speedups/snellius"
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $JOBS/g1_item3/sharp.txt)
test -n "$LINE"
exec bash $JOBS/g1_item3/run_one.sh $LINE
