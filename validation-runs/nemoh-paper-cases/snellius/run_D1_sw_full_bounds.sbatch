#!/bin/bash
#SBATCH -J np_D1_bounds
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH -o <run_root>/run_nemohpaper/logs/D1_bounds_%j.log
#SBATCH -A <account>
# D1: the SIGSEGV of J3 (sw_full: SOFTWIND half mesh ISY = 1, IRSP = 1, qtf = 3, qtf_terms = full) in the deallocation at the
# end of fs_frequency (SecondOrderDriver), at the first frequency, i.e. a write out of bounds in the free-surface pass.
# Build the same source (development_2 at 8f27716) with run-time bounds and pointer checks and a traceback, then run the
# same inputs until the check stops it; the run goes no further than the first frequencies. Expected about 0.5 h.
export RUN_ROOT=<run_root>/run_nemohpaper
export NEMOH_BIN=unused HYDROSCAL=unused NEMOH_TC=unused
source $RUN_ROOT/common.sh
SRC=$RUN/src_8f27716/src    # git archive of 8f27716 (staging), checked against its md5 list
(cd $RUN/src_8f27716 && md5sum -c md5sums_src.txt --quiet) || { echo "source differs from the staged 8f27716"; exit 1; }
B=$RUN/dbg_build
rm -rf $B && mkdir -p $B && cp -r $SRC/. $B/ && cd $B || exit 1
FILES="WavDynMods.f90 InputFiles.f90 WavDynSubs.f90 SecondOrderTypes.f90 SecondOrderInput.f90 FirstOrderState.f90 SurfaceGeometry.f90 SharpEdge.f90 PatclVelct.f90 FinGrnExtSubs.f90 InfGreen_Appr.f90 FinGreen3D.f90 SingularIntgr.f90 SingularIntgrMulti.f90 CalGreenFunc.f90 CalGreenFuncMulti.f90 BodyIntgr_irr.f90 BodyIntgr.f90 BodyIntgr_irrMulti.f90 BodyIntgrMulti.f90 HMatrix.f90 IterativeSolvers.f90 AssbMatx_irr.f90 AssbMatx.f90 AssbMatx_irrMulti.f90 AssbMatxMulti.f90 HydroStatic.f90 HydroStaticMulti.f90 NormalProcess.f90 ReadPanelMesh.f90 ReadPanelMeshMulti.f90 ImplementSubs.f90 PotentWavForce.f90 PotentWavForceMulti.f90 PressureElevation.f90 PressureElevationMulti.f90 SolveMotion.f90 SolveMotionMulti.f90 PrintOutput.f90 FarFieldDrift.f90 NearFieldDrift.f90 SecondOrderIncident.f90 PartialWaves.f90 ControlSurfaceDrift.f90 WaterlineField.f90 SurfaceVelocity.f90 FieldVelocity.f90 ControlSurfaceMesh.f90 SecondOrderOutput.f90 QTFPotential.f90 QTFField.f90 QTFFreeSurface.f90 QTFDirect.f90 QTFDirectFS.f90 QTFDirectFSBlock.f90 QTFDirectField.f90 QTFMotion.f90 QTFDriver.f90 SecondOrderDriver.f90 HAMS_Prog.f90"
ifx -qmkl -O2 -g -traceback -check bounds,pointers -heap-arrays 0 -qopenmp -o $RUN/hamsmrel_dbg $FILES > build.log 2>&1 \
  || { tail -30 build.log; exit 1; }
D=$RUN/sw_full_dbg
rm -rf $D && mkdir -p $D/input && cp $RUN/sw_full/input/* $D/input/ || exit 1
cd $D
export OMP_NUM_THREADS=24
/usr/bin/time -v -o time.txt $RUN/hamsmrel_dbg $D/input/ $D/output/ > log.txt 2>&1
echo "exit $?"
grep -n -B2 -A25 "forrtl" log.txt | head -60
