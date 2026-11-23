#!/bin/bash
#SBATCH -J np_J5C_oc3_diff
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=48G
#SBATCH -t 16:00:00
#SBATCH -o <run_root>/run_nemohpaper/logs/J5C_oc3_diff_%j.log
#SBATCH -A <account>
# J5C: full rerun of J5B (OC3 Hywind, Nemoh test case 11, difference QTF without the free surface, 50 frequencies
# 0.0628-3.14 rad/s, depth 320 m) on the branch fix/field-limes-lid-weight at SRC_SHA (50fc55e: the direct sum 9d32939, the opt-in lid weight 62e46e2, off here, and qtf_quad_terms): the field rows of C sum the
# eigenfunction series directly above nu*h = 50 (from 1.256 rad/s here, where J5B scattered through the Limes of PIDCOCK),
# and qtf_quad_terms = 1 writes the six near-field terms of every row of QTF_Quadratic.10d (QTF_QuadraticTerms.txt),
# for the desk check of the rows omega_1 <= 0.19 rad/s. Inputs: those of J5B (oc3_diff_B, md5 checked), with the one
# line "qtf_quad_terms = 1" added to SecondOrder.in (CRLF kept, checked). Builds HAMS-MREL from the staged source (md5
# checked) with the flags of the earlier binaries. J5B took 12.65 h; the direct sum costs 8.4 us per near pair against
# 26 us for PIDCOCK, so expected 12.7 h or less, limit 16 h.
export RUN_ROOT=<run_root>/run_nemohpaper
export NEMOH_BIN=unused HYDROSCAL=unused NEMOH_TC=unused
source $RUN_ROOT/common.sh
SRC_SHA=50fc55e
(cd $RUN/src_$SRC_SHA && md5sum -c md5sums_src.txt --quiet) || { echo "source differs from the staged $SRC_SHA"; exit 1; }
grep -q "SUBROUTINE EIGEN_FIELD_DIRECT" $RUN/src_$SRC_SHA/src/FinGreen3D.f90 || { echo "source has not the direct field sum"; exit 1; }
grep -q "QTF_QuadraticTerms.txt" $RUN/src_$SRC_SHA/src/QTFDriver.f90 || { echo "source has not qtf_quad_terms"; exit 1; }

make_quad_case() {
  # oc3_diff_C: the checked inputs of oc3_diff_B, with "qtf_quad_terms = 1" added at the end of SecondOrder.in, the
  # CR of its CRLF lines kept; stops unless the line is there once and every line still ends in CR
  local f=$RUN/oc3_diff_C/input/SecondOrder.in
  check_inputs oc3_diff_B || return 1
  rm -rf $RUN/oc3_diff_C && mkdir -p $RUN/oc3_diff_C/output && cp -r $RUN/oc3_diff_B/input $RUN/oc3_diff_C/ || return 1
  grep -q "qtf_quad_terms" $f && { echo "oc3_diff_C: SecondOrder.in has qtf_quad_terms already"; return 1; }
  printf 'qtf_quad_terms = 1\r\n' >> $f || return 1
  # grep has no \r escape: the CR is passed as a real character ($'\r')
  [ $(grep -c "^qtf_quad_terms = 1"$'\r$' $f) -eq 1 ] || { echo "oc3_diff_C: the qtf_quad_terms line was not added once"; return 1; }
  [ $(grep -c $'\r$' $f) -eq $(grep -c '' $f) ] || { echo "oc3_diff_C: mixed line ends after the edit"; return 1; }
  # the other four files are those of oc3_diff_B
  (cd $RUN/oc3_diff_C/input && grep -v SecondOrder.in $RUN/md5sums_oc3_diff_B.txt | md5sum -c --quiet) \
    || { echo "oc3_diff_C: inputs differ from oc3_diff_B"; return 1; }
}

make_quad_case || exit 1
B=$RUN/build_${SRC_SHA}_J5C
rm -rf $B && mkdir -p $B && cp -r $RUN/src_$SRC_SHA/src/. $B/ && cd $B || exit 1
FILES="WavDynMods.f90 InputFiles.f90 WavDynSubs.f90 SecondOrderTypes.f90 SecondOrderInput.f90 FirstOrderState.f90 SurfaceGeometry.f90 SharpEdge.f90 PatclVelct.f90 FinGrnExtSubs.f90 InfGreen_Appr.f90 FinGreen3D.f90 SingularIntgr.f90 SingularIntgrMulti.f90 CalGreenFunc.f90 CalGreenFuncMulti.f90 BodyIntgr_irr.f90 BodyIntgr.f90 BodyIntgr_irrMulti.f90 BodyIntgrMulti.f90 HMatrix.f90 IterativeSolvers.f90 AssbMatx_irr.f90 AssbMatx.f90 AssbMatx_irrMulti.f90 AssbMatxMulti.f90 HydroStatic.f90 HydroStaticMulti.f90 NormalProcess.f90 ReadPanelMesh.f90 ReadPanelMeshMulti.f90 ImplementSubs.f90 PotentWavForce.f90 PotentWavForceMulti.f90 PressureElevation.f90 PressureElevationMulti.f90 SolveMotion.f90 SolveMotionMulti.f90 PrintOutput.f90 FarFieldDrift.f90 NearFieldDrift.f90 SecondOrderIncident.f90 PartialWaves.f90 ControlSurfaceDrift.f90 WaterlineField.f90 SurfaceVelocity.f90 FieldVelocity.f90 ControlSurfaceMesh.f90 SecondOrderOutput.f90 QTFPotential.f90 QTFField.f90 QTFFreeSurface.f90 QTFDirect.f90 QTFDirectFS.f90 QTFDirectFSBlock.f90 QTFDirectField.f90 QTFMotion.f90 QTFDriver.f90 SecondOrderDriver.f90 HAMS_Prog.f90"
ifx -qmkl -O3 -heap-arrays 0 -qopenmp -traceback -o $RUN/hamsmrel_${SRC_SHA}_J5C $FILES > build.log 2>&1 || { tail -30 build.log; exit 1; }
md5sum $RUN/hamsmrel_${SRC_SHA}_J5C | tee $RUN/md5sums_bin_${SRC_SHA}_J5C.txt
BIN=$RUN/hamsmrel_${SRC_SHA}_J5C
run_hams oc3_diff_C 24
code=$?
grep -c "sum the eigenfunction series directly" $RUN/oc3_diff_C/log.txt
grep -c "may scatter" $RUN/oc3_diff_C/log.txt
ls -l $RUN/oc3_diff_C/output/Wamit_format/QTF_QuadraticTerms.txt
exit $code
