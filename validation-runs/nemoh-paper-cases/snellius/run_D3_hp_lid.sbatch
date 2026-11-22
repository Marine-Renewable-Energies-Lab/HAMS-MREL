#!/bin/bash
#SBATCH -J np_D3_hp_lid
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH -o <run_root>/run_nemohpaper/logs/D3_hp_lid_%j.log
#SBATCH -A <account>
# D3: the fall-off of A_33, A_55 and B_55 of the heave-plate cylinder above 3.3 rad/s (all of m1, m2, m3). The mesh m1 at
# 2.0-4.0 rad/s by 0.05 (41 frequencies) with (a) its lid (744 panels, 1 m rings), (b) no lid (IRSP = 0), (c) a lid of
# 0.35 m panels (2602). Same binary as hp_m1 (8f27716; first order is the same in 1515e89). Expected about 16 min.
export RUN_ROOT=<run_root>/run_nemohpaper
export NEMOH_BIN=unused HYDROSCAL=unused NEMOH_TC=unused
source $RUN_ROOT/common.sh
for c in d3_lid_m1 d3_nolid d3_lid035; do check_inputs $c || exit 1; done
fail=0
for c in d3_lid_m1 d3_nolid d3_lid035; do run_hams $c 24 || fail=1; done
exit $fail
