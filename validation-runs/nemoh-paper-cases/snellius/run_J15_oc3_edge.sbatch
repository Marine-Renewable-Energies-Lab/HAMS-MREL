#!/bin/bash
#SBATCH -J np_J15_oc3_edge
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=48G
#SBATCH -t 11:00:00
#SBATCH -o <run_root>/run_nemohpaper/logs/J15_%x_%j.log
#SBATCH -A <account>
# J15: OC3 Hywind bottom-edge mesh study, one mesh per job: submit once per case with
#   sbatch -J np_J15_oc3_e0 -t 09:00:00 --export=ALL,CASE=oc3_e0 run_J15_oc3_edge.sbatch
#   sbatch -J np_J15_oc3_e1 -t 10:00:00 --export=ALL,CASE=oc3_e1 run_J15_oc3_edge.sbatch
#   sbatch -J np_J15_oc3_e2 -t 11:00:00 --export=ALL,CASE=oc3_e2 run_J15_oc3_edge.sbatch
# Structured meshes of the OC3 hull (oc3_edge_study.py): e0 uniform (5138 hull panels), e1 and e2 with the rows next to
# the bottom edge graded from 0.375 m and 0.094 m (5306 and 5866 panels); the rest of the hull the same. 20 frequencies
# 0.0628-1.256 rad/s, so the rows omega_1 <= 0.25 against omega_2 up to 1.256. Settings of J5C (difference QTF with the
# body forcing, qtf_quad_terms = 1, lid with IRSP = 1, Hydrostatic.in of J5C), with the control surface held at J5C's
# (cs_margin = cs_bottom_margin = 4 m, cs_panel_size = 1 m: the defaults on the Nemoh mesh). Binary: J5C's
# hamsmrel_50fc55e_J5C (md5 checked). Expected 6.0 / 6.2 / 7.0 h (J5C stage times scaled), limits 9 / 10 / 11 h.
export RUN_ROOT=<run_root>/run_nemohpaper
export NEMOH_BIN=unused HYDROSCAL=unused NEMOH_TC=unused
source $RUN_ROOT/common.sh
case $CASE in oc3_e0|oc3_e1|oc3_e2) ;; *) echo "CASE must be oc3_e0, oc3_e1 or oc3_e2 (got '$CASE')"; exit 1 ;; esac
md5sum -c $RUN/md5sums_bin_50fc55e_J5C.txt --quiet || { echo "binary differs from the one J5C built"; exit 1; }
check_inputs $CASE || exit 1
BIN=$RUN/hamsmrel_50fc55e_J5C
run_hams $CASE 24
code=$?
# C must be J5C's (7370 panels); the direct-sum note once (1.256 rad/s is above nu*h = 50); no scatter warning
grep -q "QTF B2 control surface: 7370 panels on C" $RUN/$CASE/log.txt || { echo "$CASE: C is not J5C's:"; grep "QTF B2 control surface:" $RUN/$CASE/log.txt; code=1; }
grep -c "sum the eigenfunction series directly" $RUN/$CASE/log.txt
grep -c "may scatter" $RUN/$CASE/log.txt
grep "QTF B2 control surface fields" $RUN/$CASE/log.txt | awk '{s += $8} END {print "B2 field passes", s, "s"}'
ls -l $RUN/$CASE/output/Wamit_format/QTF_QuadraticTerms.txt
exit $code
