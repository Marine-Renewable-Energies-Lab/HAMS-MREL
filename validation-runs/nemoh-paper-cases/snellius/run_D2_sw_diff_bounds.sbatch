#!/bin/bash
#SBATCH -J np_D2_bounds
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH -o <run_root>/run_nemohpaper/logs/D2_bounds_%j.log
#SBATCH -A <account>
# D2: J2 (sw_diff: SOFTWIND half mesh ISY = 1, IRSP = 1, qtf_terms = body, qtf_b2 = control, 6384 points on C) ran all
# 50 frequencies and stopped with SIGSEGV at the last one (3.14 rad/s), before writing any QTF; no traceback, as the binary
# had none. The last frequency takes the other branch of control_fields (no psi states) and is followed by the quadratic
# part and the QTF output. Build the fixed source (1515e89, as J3B) with run-time bounds and pointer checks and a
# traceback, then run the inputs of sw_diff on its last two frequencies only (3.0772 and 3.14 rad/s): the first takes the
# psi branch, the second the last-frequency branch, then the QTF stages. Expected about 1 h, limit 3 h.
export RUN_ROOT=<run_root>/run_nemohpaper
export NEMOH_BIN=unused HYDROSCAL=unused NEMOH_TC=unused
source $RUN_ROOT/common.sh
(cd $RUN/src_1515e89 && md5sum -c md5sums_src.txt --quiet) || { echo "source differs from the staged 1515e89"; exit 1; }
check_inputs sw_diff || exit 1
B=$RUN/dbg_build_1515e89
rm -rf $B && mkdir -p $B && cp -r $RUN/src_1515e89/src/. $B/ && cd $B || exit 1
FILES="WavDynMods.f90 InputFiles.f90 WavDynSubs.f90 SecondOrderTypes.f90 SecondOrderInput.f90 FirstOrderState.f90 SurfaceGeometry.f90 SharpEdge.f90 PatclVelct.f90 FinGrnExtSubs.f90 InfGreen_Appr.f90 FinGreen3D.f90 SingularIntgr.f90 SingularIntgrMulti.f90 CalGreenFunc.f90 CalGreenFuncMulti.f90 BodyIntgr_irr.f90 BodyIntgr.f90 BodyIntgr_irrMulti.f90 BodyIntgrMulti.f90 HMatrix.f90 IterativeSolvers.f90 AssbMatx_irr.f90 AssbMatx.f90 AssbMatx_irrMulti.f90 AssbMatxMulti.f90 HydroStatic.f90 HydroStaticMulti.f90 NormalProcess.f90 ReadPanelMesh.f90 ReadPanelMeshMulti.f90 ImplementSubs.f90 PotentWavForce.f90 PotentWavForceMulti.f90 PressureElevation.f90 PressureElevationMulti.f90 SolveMotion.f90 SolveMotionMulti.f90 PrintOutput.f90 FarFieldDrift.f90 NearFieldDrift.f90 SecondOrderIncident.f90 PartialWaves.f90 ControlSurfaceDrift.f90 WaterlineField.f90 SurfaceVelocity.f90 FieldVelocity.f90 ControlSurfaceMesh.f90 SecondOrderOutput.f90 QTFPotential.f90 QTFField.f90 QTFFreeSurface.f90 QTFDirect.f90 QTFDirectFS.f90 QTFDirectFSBlock.f90 QTFDirectField.f90 QTFMotion.f90 QTFDriver.f90 SecondOrderDriver.f90 HAMS_Prog.f90"
ifx -qmkl -O2 -g -traceback -check bounds,pointers -heap-arrays 0 -qopenmp -o $RUN/hamsmrel_dbg_1515e89 $FILES > build.log 2>&1 \
  || { tail -30 build.log; exit 1; }
D=$RUN/sw_diff_dbg
rm -rf $D && mkdir -p $D/input && cp $RUN/sw_diff/input/* $D/input/ || exit 1
# the last two frequencies of the list: the count keeps its columns (26x, i16), the list line is read free-format;
# the input files have CRLF line ends, hence \r? in the match, and the guard stops the job if the count is not 2
awk '/Number_of_frequencies/ {sub(/ 50\r?$/, "  2"); print; getline; print " 3.0772 3.1400"; next} {print}' \
  $RUN/sw_diff/input/ControlFile.in > $D/input/ControlFile.in
grep -A1 Number_of_frequencies $D/input/ControlFile.in; grep -q "Number_of_frequencies  *2" $D/input/ControlFile.in || { echo "frequency count not set to 2"; exit 1; }
cd $D
export OMP_NUM_THREADS=24
/usr/bin/time -v -o time.txt $RUN/hamsmrel_dbg_1515e89 $D/input/ $D/output/ > log.txt 2>&1
echo "exit $?"
grep -n -B2 -A25 "forrtl" log.txt | head -60
tail -5 log.txt
