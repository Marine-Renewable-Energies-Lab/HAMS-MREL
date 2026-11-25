# Column plate, earlier mesh layout

The WAMIT-format outputs (`Wamit_format/`) of the cases `nolid` (IRSP = 0) and `lid` (IRSP = 1) of
`test/data/dipole/column_plate`, written before the dipole panels moved into `HullMesh.pnl`. In that layout the plate
was a second mesh file with its own node list and the same 96 panels; the hull mesh was the 192 panels of the column.
They were written by the integration suite of commit 0876a90 on `development_2`, built with ifort on Windows
(`/Qmkl /heap-arrays0 /fpp /Qopenmp`).

The fixtures were then converted once: the plate's nodes were appended to the hull's node list (nodes 745 to 1128), its
panels written after the hull panels under `# Start Definition of Dipole Panels`, and `NDIPOLE = 96` added to line 4.
The integration test `column_plate_one_file_layout` (`test/test_integration_dipole.f90`) runs the converted cases and
requires these outputs, line by line as text or, where another compiler prints a different last digit, every number
within 1e-5 of the largest number on its line.

If a later change of the solver changes these outputs on purpose, regenerate them from a build of that change, and say
so in its commit.
