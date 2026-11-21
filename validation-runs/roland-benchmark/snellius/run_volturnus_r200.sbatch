#!/bin/bash
#SBATCH -J roland_volturnus_r200
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=32G
#SBATCH -t 8:00:00
#SBATCH -o <run_root>/run_roland/logs/volturnus_r200_%j.log
#SBATCH -A <account>
# VolturnUS-S, the 200 m cut-off check of the off-diagonal gap against WAMIT (docs/second_order_validation.md, section 8). As
# q_volturnus_free_full_ext (same hull, lid, C, F, matrices, 46 frequencies 0.25-2.50 rad/s, differencing), but indirect
# only (no qtf_direct) and with the inner free surface to R = 200 m (WAMIT's RINNER) in 10 m panels: FSInnerMesh.pnl with
# 2640 triangles about (-12.946, 0), written by make_volturnus_r200.py. Only FS_outer of QTF_PotentialTerms.txt is
# used: the free-surface integral beyond 200 m, whose partial waves come from the hull integral (qtf_fs_circle = panels)
# and do not depend on the inner mesh. .11d of q_volturnus_free_full_ext minus it is HAMS cut at 200 m.
# Binary: run_fsspeed/hamsmrel_grn4, as q_volturnus_free_full_ext. Expected about 3.3 h, limit 8 h (COST.md); the run took 3.2 h.
module load 2024
module load intel/2024a
ulimit -s unlimited
export OMP_NUM_THREADS=24
export OMP_STACKSIZE=1G
R=<run_root>
RUN=$R/run_roland
CASE=q_volturnus_r200
BIN=$R/run_fsspeed/hamsmrel_grn4
cd $RUN/$CASE/input || exit 1
md5sum -c $RUN/md5sums_$CASE.txt --quiet || { echo "input checksums differ from the laptop"; exit 1; }
cd $RUN/$CASE
START=$(date +%s)
/usr/bin/time -v -o time.txt $BIN $RUN/$CASE/input/ $RUN/$CASE/output/ > log.txt 2>&1
CODE=$?
END=$(date +%s)
# The solver stops on a bad input with return code 0, so its last message is the check (drift_convergence.run_hams)
grep -q "completes successfully" log.txt || { [ $CODE -eq 0 ] && CODE=1; }
grep -n "QTF free surface:\|partial waves miss the field\|overflow\|warning" log.txt | head -40
PEAK=$(awk -F': ' '/Maximum resident set size/ {print $2 * 1024}' time.txt 2>/dev/null)
printf '{\n "code": %d,\n "wall_s": %d,\n "peak_bytes": %s,\n "name": "%s",\n "panels": 5076,\n "threads": 24,\n "host": "%s"\n}\n' \
  $CODE $((END - START)) "${PEAK:-null}" $CASE "$(hostname)" > result.json
cat result.json
exit $CODE
