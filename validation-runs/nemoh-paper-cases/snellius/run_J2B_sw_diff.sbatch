#!/bin/bash
#SBATCH -J np_J2B_sw_diff
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=32G
#SBATCH -t 13:00:00
#SBATCH -o <run_root>/run_nemohpaper/logs/J2B_sw_diff_%j.log
#SBATCH -A <account>
# J2B: rerun of J2 (SOFTWIND difference QTF without the free surface, 50 frequencies to 3.14 rad/s), which ran all 50
# frequencies on the binary of 8f27716 and stopped with SIGSEGV before writing a QTF. D2 ran the same inputs on 1515e89 (the
# join_states fix) with bounds and pointer checks through the last frequency, the radiation-only pass and the QTF output
# without an error. Build 1515e89 with the flags of J3B plus -traceback (no run-time cost), so a crash would name its line.
# Output to sw_diff_B, apart from the failed sw_diff. Expected 6.5 h, limit 13 h.
export RUN_ROOT=<run_root>/run_nemohpaper
export NEMOH_BIN=unused HYDROSCAL=unused NEMOH_TC=unused
source $RUN_ROOT/common.sh
(cd $RUN/src_1515e89 && md5sum -c md5sums_src.txt --quiet) || { echo "source differs from the staged 1515e89"; exit 1; }
grep -q "join_states(states, ps, sts)" $RUN/src_1515e89/src/SecondOrderDriver.f90 || { echo "source has not the fix"; exit 1; }
check_inputs sw_diff_B || exit 1
B=$RUN/build_1515e89_tb
rm -rf $B && mkdir -p $B && cp -r $RUN/src_1515e89/src/. $B/ && cd $B || exit 1
FILES="WavDynMods.f90 InputFiles.f90 WavDynSubs.f90 SecondOrderTypes.f90 SecondOrderInput.f90 FirstOrderState.f90 SurfaceGeometry.f90 SharpEdge.f90 PatclVelct.f90 FinGrnExtSubs.f90 InfGreen_Appr.f90 FinGreen3D.f90 SingularIntgr.f90 SingularIntgrMulti.f90 CalGreenFunc.f90 CalGreenFuncMulti.f90 BodyIntgr_irr.f90 BodyIntgr.f90 BodyIntgr_irrMulti.f90 BodyIntgrMulti.f90 HMatrix.f90 IterativeSolvers.f90 AssbMatx_irr.f90 AssbMatx.f90 AssbMatx_irrMulti.f90 AssbMatxMulti.f90 HydroStatic.f90 HydroStaticMulti.f90 NormalProcess.f90 ReadPanelMesh.f90 ReadPanelMeshMulti.f90 ImplementSubs.f90 PotentWavForce.f90 PotentWavForceMulti.f90 PressureElevation.f90 PressureElevationMulti.f90 SolveMotion.f90 SolveMotionMulti.f90 PrintOutput.f90 FarFieldDrift.f90 NearFieldDrift.f90 SecondOrderIncident.f90 PartialWaves.f90 ControlSurfaceDrift.f90 WaterlineField.f90 SurfaceVelocity.f90 FieldVelocity.f90 ControlSurfaceMesh.f90 SecondOrderOutput.f90 QTFPotential.f90 QTFField.f90 QTFFreeSurface.f90 QTFDirect.f90 QTFDirectFS.f90 QTFDirectFSBlock.f90 QTFDirectField.f90 QTFMotion.f90 QTFDriver.f90 SecondOrderDriver.f90 HAMS_Prog.f90"
ifx -qmkl -O3 -traceback -heap-arrays 0 -qopenmp -o $RUN/hamsmrel_1515e89_tb $FILES > build.log 2>&1 || { tail -30 build.log; exit 1; }
md5sum $RUN/hamsmrel_1515e89_tb | tee $RUN/md5sums_bin_1515e89_tb.txt
BIN=$RUN/hamsmrel_1515e89_tb
run_hams sw_diff_B 24
