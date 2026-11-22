#!/bin/bash
#SBATCH -J np_J7_hp_m3
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH -o <run_root>/run_nemohpaper/logs/J7_hp_m3_%j.log
#SBATCH -A <account>
# J7: heave-plate cylinder on the mesh refined on the plate (hp_m3, 11400 hull panels, plate-edge rows 0.2 m and
export RUN_ROOT=<run_root>/run_nemohpaper
export NEMOH_BIN=unused HYDROSCAL=unused NEMOH_TC=unused
source $RUN_ROOT/common.sh
check_inputs hp_m3 || exit 1
run_hams hp_m3 24
