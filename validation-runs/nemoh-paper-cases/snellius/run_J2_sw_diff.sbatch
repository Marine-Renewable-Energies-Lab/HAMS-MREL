#!/bin/bash
#SBATCH -J np_J2_sw_diff
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=32G
#SBATCH -t 13:00:00
#SBATCH -o <run_root>/run_nemohpaper/logs/J2_sw_diff_%j.log
#SBATCH -A <account>
# J2: SOFTWIND (Nemoh test case 10, 1830-panel half hull, ISY = 1, 56-panel lid), difference QTF with qtf_terms =
# body, 50 frequencies 0.0628-3.14 rad/s, heading 0, depth 200 m. Expected 5.2 h, limit 13 h (COST.md).
export RUN_ROOT=<run_root>/run_nemohpaper
export NEMOH_BIN=<nemoh_root>/nemoh_v303/bin
export HYDROSCAL=<nemoh_root>/nemoh_v303/build/Mesh/hydrosCal
export NEMOH_TC=<nemoh_root>/Nemoh-3.0.3/TestCases
source $RUN_ROOT/common.sh
check_inputs sw_diff || exit 1
run_hams sw_diff 24
