#!/bin/bash
#SBATCH -J np_J8_paper
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH -o <run_root>/run_nemohpaper/logs/J8_paper_%j.log
#SBATCH -A <account>
# J8: the heave-plate cylinder on a mesh of the paper's coarseness (p0: near-square panels of 0.75 m, hull 3233, lid 262;
# the paper: 3202 + 262), 0.2-4.0 rad/s by 0.05 (77 frequencies). HAMS-MREL on 17 threads and, beside it, Nemoh v3.0.3 on
# the same mesh in 7 single-thread runs of 11 frequencies each (nemoh_p0_c1 .. c7). Expected about 7 min, limit 1 h.
export RUN_ROOT=<run_root>/run_nemohpaper
export NEMOH_BIN=<nemoh_root>/nemoh_v303/bin
export HYDROSCAL=<nemoh_root>/nemoh_v303/build/Mesh/hydrosCal
export NEMOH_TC=unused
source $RUN_ROOT/common.sh
check_inputs p0 || exit 1
for k in 1 2 3 4 5 6 7; do (cd $RUN/nemoh_p0_c$k && md5sum -c $RUN/md5sums_nemoh_p0_c$k.txt --quiet) || { echo "nemoh_p0_c$k: inputs differ"; exit 1; }; done
fail=0
pids=""
for k in 1 2 3 4 5 6 7; do (run_nemoh nemoh_p0_c$k > $RUN/nemoh_p0_c$k.out 2>&1) & pids="$pids $!"; done
run_hams p0 17 || fail=1
for p in $pids; do wait $p || fail=1; done
for k in 1 2 3 4 5 6 7; do cat $RUN/nemoh_p0_c$k.out; done
exit $fail
