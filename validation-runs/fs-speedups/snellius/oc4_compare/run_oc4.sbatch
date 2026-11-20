#!/bin/bash
#SBATCH -J fsspeed_oc4
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=32G
#SBATCH -t 9:00:00
#SBATCH -o <run root>/run_fsspeed/logs/oc4_%j.log
#SBATCH -A <account>
# The OC4 WAMIT-model sum case of the old run (q_oc4wfull_sym_indsum_free_L0: qtf = 2, qtf_terms = full, 5022 hull panels,
# FSInnerMesh.pnl with 4044 inner points, depth 200 m) with every item of the fs-speedups branch at its intended default:
# field_deep_kh at its default (10, no key), psi reuse and the threaded mode check (always on), qtf_fs_circle = panels
# (the one key added), qtf_fs_inner left at stencil. The inputs are those of the old run, copied from run_oc4wsum and checked
# against md5sums_input.txt (the laptop copy, validation-runs/oc4wfull-sum/cases). Modules, threads, stack and the run
# lines as run_oc4wsum/jobs/run_hams.sbatch. At the end compare_oc4.py compares every QTF file with the old run and writes
# compare_oc4.txt beside the run (python3, standard library). Expected about 2 h (findings_speedups.md of oc4wfull-sum:
# 7.25 h -> 2.0 h with deep water from k0 h = 10, psi reuse and the threaded fit); limit 9 h (the old run took 7.25 h).
module load 2024
module load intel/2024a
ulimit -s unlimited
export OMP_NUM_THREADS=24
export OMP_STACKSIZE=1G
R="<run root>"
RUN=$R/run_fsspeed
REPO=$R/repo_fsspeed
JOBS=$REPO/validation-runs/fs-speedups/snellius
REF=$R/run_oc4wsum/q_oc4wfull_sym_indsum_free_L0
CASE=q_oc4wfull_sym_indsum_free_L0_fsspeed
BIN=$RUN/hamsmrel_grn4
mkdir -p $RUN/$CASE/input
cp $REF/input/* $RUN/$CASE/input/
cd $RUN/$CASE/input
md5sum -c $JOBS/oc4_compare/md5sums_input.txt --quiet || { echo "input checksums differ from the old run"; exit 1; }
grep -q "qtf_fs_" SecondOrder.in && { echo "SecondOrder.in of the old run already sets a qtf_fs_ key"; exit 1; }
printf 'qtf_fs_circle = panels\n' >> SecondOrder.in
cd $RUN/$CASE
git -C $REPO log -1 --format="repository at %H %s"
START=$(date +%s)
if [ -x /usr/bin/time ]; then
  /usr/bin/time -v -o time.txt $BIN $RUN/$CASE/input/ $RUN/$CASE/output/ > log.txt 2>&1
else
  $BIN $RUN/$CASE/input/ $RUN/$CASE/output/ > log.txt 2>&1
fi
CODE=$?
END=$(date +%s)
# The solver stops on a bad input with return code 0, so its last message is the check (drift_convergence.run_hams)
grep -q "completes successfully" log.txt || { [ $CODE -eq 0 ] && CODE=1; }
grep -n "qtf_fs_inner are ignored\|partial waves miss the field" log.txt && echo "check the log lines above"
PEAK=$(awk -F': ' '/Maximum resident set size/ {print $2 * 1024}' time.txt 2>/dev/null)
printf '{\n "code": %d,\n "wall_s": %d,\n "peak_bytes": %s,\n "name": "%s",\n "panels": 5022,\n "threads": 24,\n "host": "%s"\n}\n' \
  $CODE $((END - START)) "${PEAK:-null}" $CASE "$(hostname)" > result.json
cat result.json
python3 $JOBS/oc4_compare/compare_oc4.py $RUN/$CASE $REF --out $RUN/$CASE || echo "compare_oc4.py failed; run it on the laptop"
exit $CODE
