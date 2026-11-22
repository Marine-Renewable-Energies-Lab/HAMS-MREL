#!/bin/bash
#SBATCH -J np_J6_sw_fo_nolid
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH -o <run_root>/run_nemohpaper/logs/J6_sw_fo_nolid_%j.log
#SBATCH -A <account>
# J6: SOFTWIND first order without a lid (IRSP = 0), 20 frequencies 0.0628-1.256 rad/s, heading 0: the check of the
export RUN_ROOT=<run_root>/run_nemohpaper
export NEMOH_BIN=unused HYDROSCAL=unused NEMOH_TC=unused
source $RUN_ROOT/common.sh
check_inputs sw_fo_nolid || exit 1
run_hams sw_fo_nolid 24
