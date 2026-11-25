# Outputs without dipole panels, base commit

The WAMIT-format outputs (`Wamit_format/`) of the cases `image_cos` (IRSP = 0) and `image_cos_lid` (IRSP = 1) of
`test/data/image_term`. They were written by the integration suite of commit 8f27716 on `development_2`, the base of the
dipole-panel branch, built with ifort on Windows (`/Qmkl /heap-arrays0 /fpp /Qopenmp`). The integration test
`no_dipole_bit_identical` (`test/test_integration_dipole.f90`) reruns both cases and requires every line to be
identical, or, where another compiler prints a different last digit (ifx in CI), every number within 1e-5 of the
largest number on its line. A case without dipole panels must not be changed by the dipole-panel code.

If a later change of the solver changes these outputs on purpose, regenerate them from a build of that change, and say
so in its commit.

Regenerated at the merge of main's single g into `development_2` (DISPERSION from g = 9.807 to `Const_mod`'s 9.80665):
the values moved by 2e-5 to 5e-5 of the largest value of each kind (forces or moments; AmssDamp.1 and
PressureElevation.6p 2.2e-5 and 2.5e-5), the size of the change of g (3.5e-5).
