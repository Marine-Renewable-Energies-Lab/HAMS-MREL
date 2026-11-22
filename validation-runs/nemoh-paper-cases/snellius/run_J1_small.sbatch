#!/bin/bash
#SBATCH -J np_J1_small
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -o <run_root>/run_nemohpaper/logs/J1_small_%j.log
#SBATCH -A <account>
# J1: the small cases in sequence on 23 threads (heave-plate cylinder m1 and m2, cylinder of test case 8b after
# Nemoh's hydrosCal, SOFTWIND first order), with Nemoh v3.0.3 on the heave-plate m1 mesh on 1 thread beside them. Expected 1.3 h, limit 4 h (COST.md).
export RUN_ROOT=<run_root>/run_nemohpaper
export NEMOH_BIN=<nemoh_root>/nemoh_v303/bin
export HYDROSCAL=<nemoh_root>/nemoh_v303/build/Mesh/hydrosCal
export NEMOH_TC=<nemoh_root>/Nemoh-3.0.3/TestCases
source $RUN_ROOT/common.sh
fail=0
for c in hp_m1 hp_m2 sw_fo; do check_inputs $c || exit 1; done
check_inputs cyl8b_diff || exit 1
(run_nemoh nemoh_hp_m1 > $RUN/nemoh_hp_m1.out 2>&1) &
NEM=$!
for c in hp_m1 hp_m2 sw_fo; do run_hams $c 23 || fail=1; done
nemoh_mechanics 8b_QTF_Cylinder/mesh634_floating cyl8b_diff && run_hams cyl8b_diff 23 || fail=1
wait $NEM || fail=1
cat $RUN/nemoh_hp_m1.out
exit $fail
