# Shared steps of the Nemoh-paper jobs (sourced by the run_*.sbatch files). Cluster paths and the account are placeholders
# filled in at staging (see README.md): RUN_ROOT is the run folder on the cluster.
module load 2024
module load intel/2024a
ulimit -s unlimited
export OMP_STACKSIZE=1G
RUN=${RUN_ROOT:?set RUN_ROOT}
BIN=$RUN/hamsmrel                       # HAMS-MREL built from development_2 at 8f27716 (md5 in md5sums_bin.txt)
NEMOH_BIN=${NEMOH_BIN:?set NEMOH_BIN}   # Nemoh v3.0.3 programs (preProc, solver, postProc)
HYDROSCAL=${HYDROSCAL:?set HYDROSCAL}   # Nemoh v3.0.3 hydrosCal
NEMOH_TC=${NEMOH_TC:?set NEMOH_TC}      # Nemoh v3.0.3 TestCases folder

check_inputs() {
  # $1 case: the inputs must be those written on the laptop (md5sums_<case>.txt)
  (cd $RUN/$1/input && md5sum -c $RUN/md5sums_$1.txt --quiet) || { echo "$1: input checksums differ from the laptop"; return 1; }
}

run_hams() {
  # $1 case, $2 threads. The solver stops on a bad input with return code 0, so its last message is the check.
  local c=$1 n=$2 start end code
  cd $RUN/$c || return 1
  export OMP_NUM_THREADS=$n
  sed -i "s/^    Number of threads *[0-9]*/    Number of threads       $n/" input/ControlFile.in
  start=$(date +%s)
  /usr/bin/time -v -o time.txt $BIN $RUN/$c/input/ $RUN/$c/output/ > log.txt 2>&1
  code=$?
  end=$(date +%s)
  grep -q "completes successfully" log.txt || { [ $code -eq 0 ] && code=1; }
  grep -n "warning\|Warning\|error\|finite-depth Green function is kept" log.txt | head -20
  printf '{"case": "%s", "code": %d, "wall_s": %d, "threads": %d, "host": "%s"}\n' $c $code $((end - start)) $n "$(hostname)" \
    > result.json
  cat result.json
  return $code
}

nemoh_mechanics() {
  # $1 Nemoh test case folder (relative to NEMOH_TC), $2 HAMS case: Nemoh's preProc and hydrosCal on the test case as
  # shipped (it has no Mechanics), then Hydrostatic.in for HAMS from its Inertia.dat and Kh.dat
  local src=$NEMOH_TC/$1 dst=$RUN/nemoh_mech_$2
  rm -rf $dst && mkdir -p $dst && cp -r $src/. $dst/ && cd $dst || return 1
  $NEMOH_BIN/preProc . > log_preProc.txt 2>&1 || { echo "$2: preProc failed"; return 1; }
  $HYDROSCAL . > log_hydrosCal.txt 2>&1 || { echo "$2: hydrosCal failed"; return 1; }
  python3 $RUN/mechanics_to_hydrostatic.py $dst $RUN/$2/input/Hydrostatic.in || return 1
}

run_nemoh() {
  # $1 Nemoh case folder under RUN: preProc, hydrosCal, solver, postProc on one thread
  local c=$1 p start end
  cd $RUN/$c || return 1
  export OMP_NUM_THREADS=1
  start=$(date +%s)
  for p in preProc hydrosCal solver postProc; do
    if [ $p = hydrosCal ]; then $HYDROSCAL . > log_$p.txt 2>&1; else $NEMOH_BIN/$p . > log_$p.txt 2>&1; fi \
      || { echo "$c: $p returned $?"; return 1; }
  done
  end=$(date +%s)
  printf '{"case": "%s", "wall_s": %d}\n' $c $((end - start)) > result.json
  cat result.json
}
