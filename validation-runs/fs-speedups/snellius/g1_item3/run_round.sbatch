#!/bin/bash
#SBATCH -J fsspeed_g1_round
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=32G
#SBATCH -t 1:00:00
#SBATCH -o /home/vraghavan/HAMS_drift_study/run_fsspeed/logs/g1_round_%A_%a.log
#SBATCH --array=1-2%2
#SBATCH -A tdsei3290
# G1 item 3 (qtf_fs_circle = panels against samples), round m2 cases of Phase 7c-2 (qtf = 3: sum and difference,
# qtf_direct = 1, qtf_terms = full, IRSP = 1), one task per line of round.txt (route, case). Expected 0.25 h each: stage 2 (job 29336170_2, ed27bd3) took 768 s for this case, before the near-waterline fit.
# Limit 1:00:00. Same modules, threads and stack as the stage-5 job run_qtfdirectsum/jobs/stage5_wlfit.sbatch.
module load 2024
module load intel/2024a
ulimit -s unlimited
export OMP_NUM_THREADS=24
export OMP_STACKSIZE=1G
JOBS=/home/vraghavan/HAMS_drift_study/repo_fsspeed/validation-runs/fs-speedups/snellius
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $JOBS/g1_item3/round.txt)
test -n "$LINE"
exec bash $JOBS/g1_item3/run_one.sh $LINE
