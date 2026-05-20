%ComSpec% /E:ON /K ""C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 vs2019"
rem  HAMS-MREL build flags (Phase 3.5):
rem    /Qmkl       — link Intel MKL (provides multi-threaded ZGETRF/ZGETRS)
rem    /O3         — high-level optimisation
rem    /QxHost     — generate the best instruction set for the build machine (AVX/AVX2/AVX-512 as available)
rem    /Qipo       — inter-procedural optimisation across files
rem    /Qopenmp    — OpenMP support (replaces the legacy /openmp spelling)
rem    /align:array64byte — align arrays to 64 bytes for AVX-512 cache lines
rem  Notes:
rem    Do NOT enable /fp:fast — it relaxes complex-arithmetic IEEE rules and breaks bit-identity on the Green's-function kernels.
rem    /QxHost ties the binary to the build machine's ISA; for distribution builds, drop it and use /arch:AVX2 (broad compatibility) instead.
ifort /Qmkl /O3 /QxHost /Qipo /align:array64byte /exe:HAMS_ifort_Win.exe WavDynMods.f90 WavDynSubs.f90 PatclVelct.f90 FinGrnExtSubs.f90 InfGreen_Appr.f90 FinGreen3D.f90 SingularIntgr.f90 SingularIntgrMulti.f90 CalGreenFunc.f90 CalGreenFuncMulti.f90 BodyIntgr_irr.f90 BodyIntgr.f90 BodyIntgr_irrMulti.f90 BodyIntgrMulti.f90 AssbMatx_irr.f90 AssbMatx.f90 AssbMatx_irrMulti.f90 AssbMatxMulti.f90 HydroStatic.f90 HydroStaticMulti.f90 NormalProcess.f90 ReadPanelMesh.f90 ReadPanelMeshMulti.f90 ImplementSubs.f90 PotentWavForce.f90 PotentWavForceMulti.f90 PressureElevation.f90 PressureElevationMulti.f90 SolveMotion.f90 PrintOutput.f90 InputFiles.f90 HAMS_Prog.f90 /Qopenmp
pause