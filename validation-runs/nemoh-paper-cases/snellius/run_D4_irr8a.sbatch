#!/bin/bash
#SBATCH -J np_D4_irr8a
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH -o <run_root>/run_nemohpaper/logs/D4_irr8a_%j.log
#SBATCH -A <account>
# D4: the irregular-frequency removal (IRSP = 1, lid) of HAMS-MREL on Nemoh's test case 8a (cylinder R 6 m, T 14 m, 502
# hull panels, deep water), whose HydroStar reference has IRR to 6.28 rad/s: with the test case's lid (132 panels), our ring
# lid (96) and no lid, 100 frequencies 0.0628-6.28 rad/s. D3 showed the lid adding a broad error that grows above 3 rad/s
# on the heave-plate cylinder; this checks the lid on a body with a reference. Expected a few minutes.
export RUN_ROOT=<run_root>/run_nemohpaper
export NEMOH_BIN=unused HYDROSCAL=unused NEMOH_TC=unused
source $RUN_ROOT/common.sh
for c in d4_lid_nemoh d4_lid_ours d4_nolid; do check_inputs $c || exit 1; done
fail=0
for c in d4_lid_nemoh d4_lid_ours d4_nolid; do run_hams $c 24 || fail=1; done
exit $fail
