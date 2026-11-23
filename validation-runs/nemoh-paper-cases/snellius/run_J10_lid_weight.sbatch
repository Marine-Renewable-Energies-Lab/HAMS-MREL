#!/bin/bash
#SBATCH -J np_J10_lid_weight
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH -o <run_root>/run_nemohpaper/logs/J10_lid_weight_%j.log
#SBATCH -A <account>
# J10: prototype of the lid-row weighting (optional ControlFile line "Lid_row_weight ALPHA": each hull and lid row of the
# IRSP = 1 least-squares solve weighted by the square root of its panel area, the lid rows also by sqrt(ALPHA)).
# Builds HAMS-MREL from the staged source of the branch fix/field-limes-lid-weight at the lid commit (SRC_SHA = 62e46e2, after the direct-sum commit 9d32939; md5 checked),
# with the flags of the earlier binaries, then runs:
#   - heave plate p0 (paper's lid, 1.09 m, 262 rows) and p_lid (0.375 m lid, 2265 rows) at ALPHA = 1, 0.1 and 0.01;
#   - Nemoh test case 8a with the test case's lid (d4_lid_nemoh) and our ring lid (d4_lid_ours) at the same ALPHA, to check
#     that the irregular frequencies are still removed;
#   - p0 without the line (unweighted) on the new binary: must reproduce p0 of the 8f27716 binary.
# The reference results already exist: p0_nolid, p0, p_lid (unweighted, J9) and d4_* (D4). Each run gets its own folder
# <case>_w<tag>, with the inputs of <case> and the line added after "Number of threads" (CRLF kept). Expected about
# 1.1 h (p0 400 s, p_lid 650 s at 24 threads, x 3 weights, + the check run, + build), limit 3 h.
export RUN_ROOT=<run_root>/run_nemohpaper
export NEMOH_BIN=unused HYDROSCAL=unused NEMOH_TC=unused
source $RUN_ROOT/common.sh
SRC_SHA=62e46e2
(cd $RUN/src_$SRC_SHA && md5sum -c md5sums_src.txt --quiet) || { echo "source differs from the staged $SRC_SHA"; exit 1; }
grep -q "SUBROUTINE SET_ROW_WEIGHTS" $RUN/src_$SRC_SHA/src/AssbMatx_irr.f90 || { echo "source has not the lid weighting"; exit 1; }
B=$RUN/build_${SRC_SHA}_J10
rm -rf $B && mkdir -p $B && cp -r $RUN/src_$SRC_SHA/src/. $B/ && cd $B || exit 1
FILES="WavDynMods.f90 InputFiles.f90 WavDynSubs.f90 SecondOrderTypes.f90 SecondOrderInput.f90 FirstOrderState.f90 SurfaceGeometry.f90 SharpEdge.f90 PatclVelct.f90 FinGrnExtSubs.f90 InfGreen_Appr.f90 FinGreen3D.f90 SingularIntgr.f90 SingularIntgrMulti.f90 CalGreenFunc.f90 CalGreenFuncMulti.f90 BodyIntgr_irr.f90 BodyIntgr.f90 BodyIntgr_irrMulti.f90 BodyIntgrMulti.f90 HMatrix.f90 IterativeSolvers.f90 AssbMatx_irr.f90 AssbMatx.f90 AssbMatx_irrMulti.f90 AssbMatxMulti.f90 HydroStatic.f90 HydroStaticMulti.f90 NormalProcess.f90 ReadPanelMesh.f90 ReadPanelMeshMulti.f90 ImplementSubs.f90 PotentWavForce.f90 PotentWavForceMulti.f90 PressureElevation.f90 PressureElevationMulti.f90 SolveMotion.f90 SolveMotionMulti.f90 PrintOutput.f90 FarFieldDrift.f90 NearFieldDrift.f90 SecondOrderIncident.f90 PartialWaves.f90 ControlSurfaceDrift.f90 WaterlineField.f90 SurfaceVelocity.f90 FieldVelocity.f90 ControlSurfaceMesh.f90 SecondOrderOutput.f90 QTFPotential.f90 QTFField.f90 QTFFreeSurface.f90 QTFDirect.f90 QTFDirectFS.f90 QTFDirectFSBlock.f90 QTFDirectField.f90 QTFMotion.f90 QTFDriver.f90 SecondOrderDriver.f90 HAMS_Prog.f90"
ifx -qmkl -O3 -heap-arrays 0 -qopenmp -traceback -o $RUN/hamsmrel_${SRC_SHA}_J10 $FILES > build.log 2>&1 || { tail -30 build.log; exit 1; }
md5sum $RUN/hamsmrel_${SRC_SHA}_J10 | tee $RUN/md5sums_bin_${SRC_SHA}_J10.txt
BIN=$RUN/hamsmrel_${SRC_SHA}_J10

make_weighted() {
  # $1 base case, $2 new case, $3 ALPHA ("none": no line). Copies the checked inputs and adds the line after the
  # "Number of threads" line, keeping the CRLF line ends of the laptop files; stops if the line is not there once.
  local base=$1 new=$2 a=$3
  check_inputs $base || return 1
  rm -rf $RUN/$new && mkdir -p $RUN/$new/output && cp -r $RUN/$base/input $RUN/$new/ || return 1
  [ "$a" = none ] && return 0
  awk -v a="$a" '{print} /^ *Number of threads/ {cr = ($0 ~ /\r$/) ? "\r" : ""; printf "    Lid_row_weight          %s%s\n", a, cr}' \
    $RUN/$base/input/ControlFile.in > $RUN/$new/input/ControlFile.in || return 1
  # grep has no \r escape: the CR is passed as a real character ($'\r'), the dot of the weight escaped
  [ $(grep -cE "^ *Lid_row_weight +${a//./\\.}"$'\r?$' $RUN/$new/input/ControlFile.in) -eq 1 ] \
    || { echo "$new: the Lid_row_weight line was not added once"; return 1; }
  [ $(grep -c $'\r$' $RUN/$new/input/ControlFile.in) -eq $(grep -c '' $RUN/$new/input/ControlFile.in) ] \
    || { echo "$new: mixed line ends after the edit"; return 1; }
}

fail=0
CASES=""
make_weighted p0 p0_wnone none && CASES="$CASES p0_wnone" || fail=1
for w in 1 0.1 0.01; do
  t=${w/./}
  for c in d4_lid_nemoh d4_lid_ours p0 p_lid; do
    make_weighted $c ${c}_w$t $w && CASES="$CASES ${c}_w$t" || fail=1
  done
done
[ $fail -eq 0 ] || exit 1
for c in $CASES; do
  run_hams $c 24 || fail=1
  grep -n "Lid_row_weight" $RUN/$c/log.txt | head -2
done
exit $fail
