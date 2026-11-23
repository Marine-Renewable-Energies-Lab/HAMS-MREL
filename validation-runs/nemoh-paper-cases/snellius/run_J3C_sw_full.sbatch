#!/bin/bash
#SBATCH -J np_J3C_sw_full
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=64G
#SBATCH -t 22:00:00
#SBATCH -o <run_root>/run_nemohpaper/logs/J3C_sw_full_%j.log
#SBATCH -A <account>
# J3C: rerun of J3B (SOFTWIND, QTF- and QTF+ with the free-surface integral, 32 frequencies to 2.01 rad/s, the inputs of
# sw_full_B: qtf_fs_panel = 0.76 m) on the branch fix/field-limes-lid-weight at SRC_SHA (62e46e2: the direct sum, 9d32939, then the lid weighting, off without its line), which sums the eigenfunction
# series directly in the field rows above nu*h = 50 (1.57 rad/s here, where J3B's QTFs broke down through the Limes of
# PIDCOCK). The build has -traceback and the flags of the earlier binaries. Output to sw_full_C; the log must name the
# direct field sum once ("sum the eigenfunction series directly") and the depth guard ("finite-depth Green function is
# kept"). J3B took 17.4 h; the new rows cost 8 us per pair against 26 us for PIDCOCK (desk timing), so expected 16.5 to
# 17.4 h, limit 22 h.
export RUN_ROOT=<run_root>/run_nemohpaper
export NEMOH_BIN=unused HYDROSCAL=unused NEMOH_TC=unused
source $RUN_ROOT/common.sh
SRC_SHA=62e46e2
(cd $RUN/src_$SRC_SHA && md5sum -c md5sums_src.txt --quiet) || { echo "source differs from the staged $SRC_SHA"; exit 1; }
grep -q "SUBROUTINE EIGEN_FIELD_DIRECT" $RUN/src_$SRC_SHA/src/FinGreen3D.f90 || { echo "source has not the direct field sum"; exit 1; }
grep -q "join_states(states, ps, sts)" $RUN/src_$SRC_SHA/src/SecondOrderDriver.f90 || { echo "source has not the join_states fix"; exit 1; }
check_inputs sw_full_B || exit 1
rm -rf $RUN/sw_full_C && mkdir -p $RUN/sw_full_C/output && cp -r $RUN/sw_full_B/input $RUN/sw_full_C/ || exit 1
B=$RUN/build_${SRC_SHA}_J3C
rm -rf $B && mkdir -p $B && cp -r $RUN/src_$SRC_SHA/src/. $B/ && cd $B || exit 1
FILES="WavDynMods.f90 InputFiles.f90 WavDynSubs.f90 SecondOrderTypes.f90 SecondOrderInput.f90 FirstOrderState.f90 SurfaceGeometry.f90 SharpEdge.f90 PatclVelct.f90 FinGrnExtSubs.f90 InfGreen_Appr.f90 FinGreen3D.f90 SingularIntgr.f90 SingularIntgrMulti.f90 CalGreenFunc.f90 CalGreenFuncMulti.f90 BodyIntgr_irr.f90 BodyIntgr.f90 BodyIntgr_irrMulti.f90 BodyIntgrMulti.f90 HMatrix.f90 IterativeSolvers.f90 AssbMatx_irr.f90 AssbMatx.f90 AssbMatx_irrMulti.f90 AssbMatxMulti.f90 HydroStatic.f90 HydroStaticMulti.f90 NormalProcess.f90 ReadPanelMesh.f90 ReadPanelMeshMulti.f90 ImplementSubs.f90 PotentWavForce.f90 PotentWavForceMulti.f90 PressureElevation.f90 PressureElevationMulti.f90 SolveMotion.f90 SolveMotionMulti.f90 PrintOutput.f90 FarFieldDrift.f90 NearFieldDrift.f90 SecondOrderIncident.f90 PartialWaves.f90 ControlSurfaceDrift.f90 WaterlineField.f90 SurfaceVelocity.f90 FieldVelocity.f90 ControlSurfaceMesh.f90 SecondOrderOutput.f90 QTFPotential.f90 QTFField.f90 QTFFreeSurface.f90 QTFDirect.f90 QTFDirectFS.f90 QTFDirectFSBlock.f90 QTFDirectField.f90 QTFMotion.f90 QTFDriver.f90 SecondOrderDriver.f90 HAMS_Prog.f90"
ifx -qmkl -O3 -heap-arrays 0 -qopenmp -traceback -o $RUN/hamsmrel_${SRC_SHA}_J3C $FILES > build.log 2>&1 || { tail -30 build.log; exit 1; }
md5sum $RUN/hamsmrel_${SRC_SHA}_J3C | tee $RUN/md5sums_bin_${SRC_SHA}_J3C.txt
BIN=$RUN/hamsmrel_${SRC_SHA}_J3C
run_hams sw_full_C 24
code=$?
grep -c "sum the eigenfunction series directly" $RUN/sw_full_C/log.txt
grep -c "may scatter" $RUN/sw_full_C/log.txt
exit $code
