#!/bin/bash
#SBATCH -J np_J14_nemoh_oc3
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH -o <run_root>/run_nemohpaper/logs/J14_nemoh_oc3_%j.log
#SBATCH -A <account>
# J14: Nemoh v3.0.3 on its own OC3 Hywind test case (11_QTF_OC3_Hywind/mesh4842_Floating: 4842 panels with its lid,
# h = 320 m), difference QTF DUOK + HASBO (Contrib 2, no free surface, as shipped), as a third opinion on the rows
# omega_1 <= 0.25 rad/s where HAMS-MREL (J5B, J5C) and HydroStar differ by up to 65 % of max|QTF|. As shipped except the
# frequency list: 20 frequencies 0.0628-1.256 rad/s for the first order and the QTF (the shipped list is 65 to 4.082),
# so every difference is on the list; the files in LF (oc3_edge_study.py). Nemoh runs on one thread (no OpenMP in the
# build); fat_genoa bills 24 CPUs at least. Inputs staged as nemoh_oc3_input (md5 checked), copied to a fresh nemoh_oc3.
# Expected 3-5 h (scaled from the Nemoh runs J8/J9, deep water, 11 frequencies in 6-13 min at 3400-4700 panels, times a
# finite-depth factor), limit 12 h. Writes stages.txt with the time of each program.
export RUN_ROOT=<run_root>/run_nemohpaper
export NEMOH_BIN=<nemoh_root>/nemoh_v303/bin
export HYDROSCAL=<nemoh_root>/nemoh_v303/build/Mesh/hydrosCal
export NEMOH_TC=unused
source $RUN_ROOT/common.sh
C=nemoh_oc3
(cd $RUN/${C}_input && md5sum -c $RUN/md5sums_$C.txt --quiet) || { echo "$C: inputs differ from the laptop"; exit 1; }
rm -rf $RUN/$C && cp -r $RUN/${C}_input $RUN/$C && cd $RUN/$C || exit 1
export OMP_NUM_THREADS=1
: > stages.txt
for p in preProc hydrosCal solver postProc QTFpreProc QTFsolver QTFpostProc; do
  s=$(date +%s)
  if [ $p = hydrosCal ]; then $HYDROSCAL . > log_$p.txt 2>&1; else $NEMOH_BIN/$p . > log_$p.txt 2>&1; fi
  code=$?
  e=$(date +%s)
  echo "$p $((e - s)) s, code $code" | tee -a stages.txt
  [ $code -eq 0 ] || { tail -20 log_$p.txt; exit 1; }
done
# Nemoh can return 0 on a failure: the difference QTF file is the check
[ -s results/QTF/OUT_QTFM_N.dat ] || { echo "$C: no results/QTF/OUT_QTFM_N.dat"; ls -l results/QTF; exit 1; }
ls -l results/QTF results/RAO.tec
printf '{"case": "%s", "wall_s": %d}\n' $C $(awk '{s += $2} END {print s}' stages.txt) > result.json
cat result.json
