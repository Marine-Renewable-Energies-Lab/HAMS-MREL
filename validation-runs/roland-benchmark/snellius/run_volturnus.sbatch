#!/bin/bash
#SBATCH -J roland_volturnus
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=32G
#SBATCH -t 18:00:00
#SBATCH -o <run_root>/run_roland/logs/volturnus_%j.log
#SBATCH -A <account>
# VolturnUS-S (IEA 15 MW) matched to the second-order WAMIT run of IEA-15-240-RWT v1.1.18: difference QTF with the
# free-surface integral (qtf_terms = full), indirect and direct methods (qtf_direct = 1), free body, heading 0, depth 200 m,
# 34 frequencies 0.25-1.90 rad/s, 5076 hull panels (the IEA low-order mesh without its lid, mirrored), FSInnerMesh.pnl with
# 4412 triangles. Inputs: make_volturnus_case.py on the laptop, copied to $RUN/$CASE/input by scp and checked against
# md5sums_input.txt. Binary: run_fsspeed/hamsmrel_grn4 of the fs-speedups branch (built from 4cbd625; its source is that of
# 3af6cf8). The later commits of that branch up to 80eaf3d change log lines and make qtf_fs_circle = panels the default; the
# job sets that key, so the run is that of 80eaf3d. Expected about 5 h, limit 18 h (COST.md).
module load 2024
module load intel/2024a
ulimit -s unlimited
export OMP_NUM_THREADS=24
export OMP_STACKSIZE=1G
R=<run_root>
RUN=$R/run_roland
CASE=q_volturnus_free_full
BIN=$R/run_fsspeed/hamsmrel_grn4
cd $RUN/$CASE/input || exit 1
md5sum -c $RUN/md5sums_input.txt --quiet || { echo "input checksums differ from the laptop"; exit 1; }
grep -q "qtf_fs_" SecondOrder.in && { echo "SecondOrder.in already sets a qtf_fs_ key"; exit 1; }
printf 'qtf_fs_circle = panels\n' >> SecondOrder.in
cd $RUN/$CASE
START=$(date +%s)
/usr/bin/time -v -o time.txt $BIN $RUN/$CASE/input/ $RUN/$CASE/output/ > log.txt 2>&1
CODE=$?
END=$(date +%s)
# The solver stops on a bad input with return code 0, so its last message is the check (drift_convergence.run_hams)
grep -q "completes successfully" log.txt || { [ $CODE -eq 0 ] && CODE=1; }
grep -n "qtf_fs_inner are ignored\|partial waves miss the field\|warning" log.txt | head -40
PEAK=$(awk -F': ' '/Maximum resident set size/ {print $2 * 1024}' time.txt 2>/dev/null)
printf '{\n "code": %d,\n "wall_s": %d,\n "peak_bytes": %s,\n "name": "%s",\n "panels": 5076,\n "threads": 24,\n "host": "%s"\n}\n' \
  $CODE $((END - START)) "${PEAK:-null}" $CASE "$(hostname)" > result.json
cat result.json
exit $CODE
