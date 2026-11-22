#!/bin/bash
#SBATCH -J np_J4_oc4_bi
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=32G
#SBATCH -t 03:00:00
#SBATCH -o <run_root>/run_nemohpaper/logs/J4_oc4_bi_%j.log
#SBATCH -A <account>
# J4: OC4 (Nemoh test case 9, 2068 hull panels, 112-panel lid, C and F from make_mesh), difference QTF of every
# heading pair of 0 and 30 degrees (mean_drift = 2), 50 frequencies 0.0628-3.14 rad/s. Expected 0.6 h, limit 3 h.
export RUN_ROOT=<run_root>/run_nemohpaper
export NEMOH_BIN=<nemoh_root>/nemoh_v303/bin
export HYDROSCAL=<nemoh_root>/nemoh_v303/build/Mesh/hydrosCal
export NEMOH_TC=<nemoh_root>/Nemoh-3.0.3/TestCases
source $RUN_ROOT/common.sh
check_inputs oc4_bi || exit 1
run_hams oc4_bi 24
