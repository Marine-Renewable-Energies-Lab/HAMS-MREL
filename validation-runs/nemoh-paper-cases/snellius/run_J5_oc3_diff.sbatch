#!/bin/bash
#SBATCH -J np_J5_oc3_diff
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=32G
#SBATCH -t 22:00:00
#SBATCH -o <run_root>/run_nemohpaper/logs/J5_oc3_diff_%j.log
#SBATCH -A <account>
# J5: OC3 Hywind (Nemoh test case 11, 4710 hull panels, 96-panel lid), difference QTF, 50 frequencies
# 0.0628-3.14 rad/s, depth 320 m; Hydrostatic.in from Nemoh's hydrosCal (the test case has no Mechanics). Expected 8.8 h,
# limit 22 h (COST.md). Option J5b: 32 frequencies to 2.01 rad/s (edit ControlFile.in on the laptop, new md5sums).
export RUN_ROOT=<run_root>/run_nemohpaper
export NEMOH_BIN=<nemoh_root>/nemoh_v303/bin
export HYDROSCAL=<nemoh_root>/nemoh_v303/build/Mesh/hydrosCal
export NEMOH_TC=<nemoh_root>/Nemoh-3.0.3/TestCases
source $RUN_ROOT/common.sh
check_inputs oc3_diff || exit 1
nemoh_mechanics 11_QTF_OC3_Hywind/mesh4842_Floating oc3_diff || exit 1
run_hams oc3_diff 24
