#!/bin/bash
# One G1 item-3 run: run_one.sh <samples|panels> <case>. Writes the case under run_fsspeed/g1_<route>/<case> (make_case.py),
# runs HAMS as run_oc4wsum/jobs/run_hams.sbatch does (input and output folders with a trailing slash, log.txt, time.txt,
# result.json), and checks that the route reached the log. Called by run_sharp.sbatch and run_round.sbatch.
set -u
R=/home/vraghavan/HAMS_drift_study
RUN=$R/run_fsspeed
REPO=$R/repo_fsspeed
JOBS=$REPO/validation-runs/fs-speedups/snellius
BIN=$RUN/hamsmrel_grn4
ROUTE=$1
CASE=$2
ROOT=$RUN/g1_$ROUTE
mkdir -p $ROOT
IN=$(python3 $JOBS/g1_item3/make_case.py $REPO $ROOT $CASE $ROUTE --threads 24) || { echo "make_case.py failed"; exit 1; }
grep -q "^qtf_fs_circle = $ROUTE\$" $IN/SecondOrder.in || { echo "qtf_fs_circle not set"; exit 1; }
cd $ROOT/$CASE
git -C $REPO log -1 --format="repository at %H %s"
START=$(date +%s)
if [ -x /usr/bin/time ]; then
  /usr/bin/time -v -o time.txt $BIN $IN/ $ROOT/$CASE/output/ > log.txt 2>&1
else
  $BIN $IN/ $ROOT/$CASE/output/ > log.txt 2>&1
fi
CODE=$?
END=$(date +%s)
# The solver stops on a bad input with return code 0, so its last message is the check (drift_convergence.run_hams)
grep -q "completes successfully" log.txt || { [ $CODE -eq 0 ] && CODE=1; }
# The panels route must not be reported as ignored, and its mode check must not warn
grep -n "qtf_fs_inner are ignored\|partial waves miss the field" log.txt && echo "check the log lines above"
PEAK=$(awk -F': ' '/Maximum resident set size/ {print $2 * 1024}' time.txt 2>/dev/null)
printf '{\n "code": %d,\n "wall_s": %d,\n "peak_bytes": %s,\n "name": "%s",\n "route": "%s",\n "threads": 24,\n "host": "%s"\n}\n' \
  $CODE $((END - START)) "${PEAK:-null}" $CASE $ROUTE "$(hostname)" > result.json
cat result.json
exit $CODE
