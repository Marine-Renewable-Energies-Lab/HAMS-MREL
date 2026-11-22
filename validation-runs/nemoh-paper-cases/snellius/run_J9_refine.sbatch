#!/bin/bash
#SBATCH -J np_J9_refine
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=96G
#SBATCH -t 10:00:00
#SBATCH -o <run_root>/run_nemohpaper/logs/J9_refine_%j.log
#SBATCH -A <account>
# J9: one-part-at-a-time refinements of the heave-plate mesh p0 (paper's coarseness), 0.2-4.0 rad/s by 0.05 (77
# frequencies). After D3 (HAMS's lid spoils A, B above 3 rad/s on this body, more with more lid rows) and D4 (the lid works
# on the plain cylinder 8a): HAMS-MREL on each hull refinement with the paper's lid and without one (_nolid, IRSP = 0),
# p0 without a lid, lids of 0.375 m (p_lid) and 2.2 m (p_lid2); Nemoh v3.0.3 on each refinement with the lid, 7 chunks
# of 11 frequencies each, 7 single-thread runs at a time beside HAMS (17 threads while Nemoh runs, 24 after; the largest
# cases last). Expected about 4.5 h, limit 10 h.
export RUN_ROOT=<run_root>/run_nemohpaper
export NEMOH_BIN=<nemoh_root>/nemoh_v303/bin
export HYDROSCAL=<nemoh_root>/nemoh_v303/build/Mesh/hydrosCal
export NEMOH_TC=unused
source $RUN_ROOT/common.sh
HAMS_CASES="p_offset_nolid p0_nolid p_lid2 p_offset p_wl_nolid p_wl p_edge_nolid p_edge p_lid p_wall_nolid p_wall p_plate"
NEMOH_CASES=""
for v in p_plate p_wall p_lid p_edge p_wl p_offset; do for k in 1 2 3 4 5 6 7; do NEMOH_CASES="$NEMOH_CASES nemoh_${v}_c$k"; done; done
for c in $HAMS_CASES; do check_inputs $c || exit 1; done
for c in $NEMOH_CASES; do (cd $RUN/$c && md5sum -c $RUN/md5sums_$c.txt --quiet) || { echo "$c: inputs differ"; exit 1; }; done
export -f run_nemoh
export RUN NEMOH_BIN HYDROSCAL
# Nemoh: 7 at a time, in the order of NEMOH_CASES (the largest meshes first)
(echo $NEMOH_CASES | tr ' ' '\n' | grep . | \
   xargs -P 7 -I{} bash -c 'run_nemoh {} > $RUN/{}.out 2>&1 || echo "{} failed"') > $RUN/J9_nemoh.out 2>&1 &
NEM=$!
fail=0
for c in $HAMS_CASES; do
  if kill -0 $NEM 2>/dev/null; then t=17; else t=24; fi
  run_hams $c $t || fail=1
done
wait $NEM
cat $RUN/J9_nemoh.out
grep -l "failed\|returned" $RUN/nemoh_p_*.out && fail=1
for c in $NEMOH_CASES; do cat $RUN/$c/result.json 2>/dev/null || { echo "$c: no result"; fail=1; }; done
exit $fail
