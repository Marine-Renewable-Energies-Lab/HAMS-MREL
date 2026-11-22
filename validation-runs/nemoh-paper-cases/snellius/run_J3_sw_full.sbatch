#!/bin/bash
#SBATCH -J np_J3_sw_full
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=64G
#SBATCH -t 26:00:00
#SBATCH -o <run_root>/run_nemohpaper/logs/J3_sw_full_%j.log
#SBATCH -A <account>
# J3: SOFTWIND, difference and sum QTF with the free-surface integral (qtf = 3, qtf_terms = full), 32 frequencies
# 0.0628-2.0096 rad/s (the grid of HydroStar's SPG files), heading 0. Expected 10 h, limit 26 h (COST.md).
export RUN_ROOT=<run_root>/run_nemohpaper
export NEMOH_BIN=<nemoh_root>/nemoh_v303/bin
export HYDROSCAL=<nemoh_root>/nemoh_v303/build/Mesh/hydrosCal
export NEMOH_TC=<nemoh_root>/Nemoh-3.0.3/TestCases
source $RUN_ROOT/common.sh
check_inputs sw_full || exit 1
run_hams sw_full 24
