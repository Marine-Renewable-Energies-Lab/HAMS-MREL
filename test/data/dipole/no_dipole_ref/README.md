# Outputs without dipole panels, base commit

The WAMIT-format outputs (`Wamit_format/`) of the cases `image_cos` (IRSP = 0) and `image_cos_lid` (IRSP = 1) of
`test/data/image_term`. They were written by the integration suite of commit 8f27716 on `development_2`, the base of the
dipole-panel branch, built with ifort on Windows (`/Qmkl /heap-arrays0 /fpp /Qopenmp`). The integration test
`no_dipole_bit_identical` (`test/test_integration_dipole.f90`) reruns both cases and requires every line to be
identical. A case without `DipoleMesh.pnl` must not be changed by the dipole-panel code.

If a later change of the solver changes these outputs on purpose, regenerate them from a build of that change, and say
so in its commit.
