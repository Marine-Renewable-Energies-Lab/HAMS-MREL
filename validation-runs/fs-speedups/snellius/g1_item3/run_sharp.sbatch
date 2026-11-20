#!/bin/bash
#SBATCH -J fsspeed_g1_sharp
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=32G
#SBATCH -t 6:00:00
#SBATCH -o /home/vraghavan/HAMS_drift_study/run_fsspeed/logs/g1_sharp_%A_%a.log
#SBATCH --array=1-4%4
#SBATCH -A tdsei3290
# G1 item 3 (qtf_fs_circle = panels against samples), sharp m2 cases of Phase 7c-2 (qtf = 3: sum and difference,
# qtf_direct = 1, qtf_terms = full, IRSP = 1), one task per line of sharp.txt (route, case). Expected 3.2 h each: stage 5 (job 29405797, 47d5911) took 11436 s fixed and 11394 s free with the same case.
# Limit 6:00:00. Same modules, threads and stack as the stage-5 job run_qtfdirectsum/jobs/stage5_wlfit.sbatch.
module load 2024
module load intel/2024a
ulimit -s unlimited
export OMP_NUM_THREADS=24
export OMP_STACKSIZE=1G
JOBS=/home/vraghavan/HAMS_drift_study/repo_fsspeed/validation-runs/fs-speedups/snellius
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $JOBS/g1_item3/sharp.txt)
test -n "$LINE"
exec bash $JOBS/g1_item3/run_one.sh $LINE
