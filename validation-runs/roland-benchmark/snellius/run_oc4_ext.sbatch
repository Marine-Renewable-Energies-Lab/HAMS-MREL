#!/bin/bash
#SBATCH -J roland_oc4_ext
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=32G
#SBATCH -t 22:00:00
#SBATCH -o <run_root>/run_roland/logs/oc4_ext_%j.log
#SBATCH -A <account>
# OC4 full WAMIT hull, 9208 panels, over the range of Roland et al. (2026): difference QTF, qtf_terms = body (as the
# existing 0.25-1.90 run), free body, heading 0, depth 200 m, 47 frequencies 0.25-2.55 rad/s, IRSP = 1 with a 1674-panel
# lid, velocity_method = differencing. Inputs: make_oc4_extension.py on the laptop (meshes and matrices of the existing
# case), copied to $RUN/$CASE/input by scp and checked against md5sums_$CASE.txt. Binary as run_volturnus.sbatch (the
# qtf_fs_circle key it appends has no effect without qtf_terms = full). Expected about 8.4 h, limit 22 h (COST.md).
module load 2024
module load intel/2024a
ulimit -s unlimited
export OMP_NUM_THREADS=24
export OMP_STACKSIZE=1G
R=<run_root>
RUN=$R/run_roland
CASE=q_oc4wfull_sym_free_L1_ext
BIN=$R/run_fsspeed/hamsmrel_grn4
cd $RUN/$CASE/input || exit 1
md5sum -c $RUN/md5sums_$CASE.txt --quiet || { echo "input checksums differ from the laptop"; exit 1; }
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
printf '{\n "code": %d,\n "wall_s": %d,\n "peak_bytes": %s,\n "name": "%s",\n "panels": 9208,\n "threads": 24,\n "host": "%s"\n}\n' \
  $CODE $((END - START)) "${PEAK:-null}" $CASE "$(hostname)" > result.json
cat result.json
exit $CODE
