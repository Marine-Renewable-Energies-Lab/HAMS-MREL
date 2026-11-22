#!/bin/bash
#SBATCH -J np_J5B_oc3_diff
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=48G
#SBATCH -t 22:00:00
#SBATCH -o <run_root>/run_nemohpaper/logs/J5B_oc3_diff_%j.log
#SBATCH -A <account>
# J5B: rerun of J5 (OC3 Hywind, Nemoh test case 11, difference QTF without the free surface, 50 frequencies 0.0628-3.14
# rad/s, depth 320 m), which ran all 50 frequencies on the binary of 8f27716 and stopped with SIGSEGV in the automatic
# deallocation of control_fields (addr2line) before writing a QTF, as J2 did. Binary of 1515e89 with -traceback, built by
# J2B (md5 checked). Inputs: those of the laptop and the Hydrostatic.in that J5 wrote from Nemoh's hydrosCal. Output to
# oc3_diff_B, apart from the failed oc3_diff. Expected 12.6 h (J5's run), limit 22 h.
export RUN_ROOT=<run_root>/run_nemohpaper
export NEMOH_BIN=unused HYDROSCAL=unused NEMOH_TC=unused
source $RUN_ROOT/common.sh
md5sum -c $RUN/md5sums_bin_1515e89_tb.txt --quiet || { echo "binary differs from the one J2B built"; exit 1; }
check_inputs oc3_diff_B || exit 1
BIN=$RUN/hamsmrel_1515e89_tb
run_hams oc3_diff_B 24
