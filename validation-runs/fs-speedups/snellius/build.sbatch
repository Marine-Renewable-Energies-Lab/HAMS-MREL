#!/bin/bash
#SBATCH -J fsspeed_build
#SBATCH -p fat_genoa
#SBATCH -N 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem=32G
#SBATCH -t 0:30:00
#SBATCH -o <run root>/run_fsspeed/logs/build_%j.log
#SBATCH -A <account>
# The HAMS binary of the fs-speedups branch (feature/qtf-fs-speedups) for the G1 item-3 runs and the OC4 comparison.
# repo_fsspeed is a clone of the git bundle of the branch (README.md, section 2). The source must be that of 3af6cf8 (the
# commits after it change docs and validation-runs only): the job stops otherwise. Modules, flags and file list as
# run_oc4wsum/jobs/build.sbatch (the build of the old OC4 run); the list is every src/*.f90 of 3af6cf8 in
# dependency order (checked on the laptop: the same files as at 22086da, no module used before it is compiled). The job
# checks that every file of the list is in src/ and that src/ has no .f90 file the list leaves out, naming the file.
# GRN_KIND = 4 only: the cases of both job sets use the GRN_KIND = 4 binary (drift_convergence.py Case.grn = 4).
module load 2024
module load intel/2024a
ulimit -s unlimited
set -e
R="<run root>"
REPO=$R/repo_fsspeed
RUN=$R/run_fsspeed
C1=3af6cf85b743a1c7e8d15eacadc0f1f302ac6ca2
git -C $REPO merge-base --is-ancestor $C1 HEAD
git -C $REPO diff --quiet $C1 HEAD -- src
FLAGS="-qmkl -qopenmp"
LIST="WavDynMods.f90 PatclVelct.f90 BodyIntgr.f90 HMatrix.f90 IterativeSolvers.f90 AssbMatx.f90 BodyIntgrMulti.f90 AssbMatxMulti.f90 BodyIntgr_irr.f90 AssbMatx_irr.f90 BodyIntgr_irrMulti.f90 AssbMatx_irrMulti.f90 FinGrnExtSubs.f90 FinGreen3D.f90 InfGreen_Appr.f90 SingularIntgr.f90 CalGreenFunc.f90 SingularIntgrMulti.f90 CalGreenFuncMulti.f90 SecondOrderTypes.f90 FirstOrderState.f90 NearFieldDrift.f90 ControlSurfaceDrift.f90 SurfaceGeometry.f90 ControlSurfaceMesh.f90 FarFieldDrift.f90 WaterlineField.f90 SurfaceVelocity.f90 FieldVelocity.f90 HydroStatic.f90 HydroStaticMulti.f90 ImplementSubs.f90 InputFiles.f90 NormalProcess.f90 PartialWaves.f90 PotentWavForce.f90 PotentWavForceMulti.f90 PressureElevation.f90 PressureElevationMulti.f90 PrintOutput.f90 SecondOrderIncident.f90 SecondOrderOutput.f90 QTFField.f90 QTFFreeSurface.f90 QTFPotential.f90 SharpEdge.f90 QTFDirect.f90 QTFDirectFS.f90 QTFDirectFSBlock.f90 QTFDirectField.f90 QTFDriver.f90 ReadPanelMesh.f90 ReadPanelMeshMulti.f90 SecondOrderDriver.f90 SecondOrderInput.f90 SolveMotion.f90 SolveMotionMulti.f90 WavDynSubs.f90 HAMS_Prog.f90"
d=$RUN/build_grn4; rm -rf $d; mkdir -p $d; cp $REPO/src/*.f90 $d/; cd $d
for f in $LIST; do [ -f "$f" ] || { echo "build: $f of the file list is not in src/"; exit 1; }; done
for f in *.f90; do case " $LIST " in *" $f "*) ;; *) echo "build: src/$f is not in the file list"; exit 1;; esac; done
grep -q "GRN_KIND = 4" WavDynMods.f90
ifx $FLAGS -o $RUN/hamsmrel_grn4 $LIST
git -C $REPO log -1 --format="built from %H %s"
