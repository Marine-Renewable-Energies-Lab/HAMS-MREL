# Input Files

HAMS-MREL requires a set of input files to run. The number files depend on whether you are analyzing a single body or multiple bodies:
- `ControlFile.in`
- `Hydrostatic[_i].in`
- `HullMesh[_i].pnl`
- `WaterPlaneMesh[_i].pnl` (only if analyzing irregular frequencies)
- `gen_mod_[i]_[dof].txt` (optional, for generalized modes)

For multi-body cases, `_i` denotes the body number (`i = 1 … numberofbodies`). For a single body, no index is used (`HullMesh.pnl`, `Hydrostatic.in` and `WaterPlaneMesh.pnl`).

The program will check that all required files exist before running. If any file is missing, execution will stop.

### Mesh Files

As described by the [HAMS user manual](https://github.com/YingyiLiu/HAMS/blob/master/Manual/A%20brief%20manual.pdf), to obtain the water and body meshes in a HAMS-compatible format, you can start by exporting the hydrodynamic CAD model from Rhinoceros to WAMIT’s *.gdf format. Next, use the built-in tool WAMIT_MeshTran.exe to convert the *.gdf mesh into the HAMS mesh format. To do this, run RunWAMIT_MeshTran.bat and enter the filename of the *.gdf file when prompted. Ensure that the *.gdf file includes meshes for both the waterplane and the submerged bodies. The tool will automatically separate them into two files: WaterplaneMesh.pnl and HullMesh.pnl. Note that the WAMIT_MeshTran tool is currently only available on Windows.

### Fortran Format Descriptors

HAMS-MREL uses Fortran format descriptors to read values from the file. These are indicated in the Format column in the tables below:
- `x`: skip spaces. For example, `14x` means skip 14 characters.
- `iN`: integer with width N. Example: `i16` reads a 16-character wide integer.
- `fW.D`: floating-point number with width W and D digits after the decimal point. Example: `f30.15` reads a number in 30 characters with 15 digits after the decimal.
- `A100`: alphanumeric string of length 100.

**Example:** `20x,f30.15` means skip the first 20 characters of the line, then read a floating-point number using 30-character width and 15 digits after the decimal point.

**Example:** `FORMAT(6(2x,E12.5))` reads 6 numbers per line, each in 12-character wide scientific notation with 5 digits after the decimal, preceded by 2 spaces.

## ControlFile.in

The general settings of the simulation are specified in `ControlFile.in`. The file is read by the subroutine `ReadControlFile` in [InputFiles.f90](../src/InputFiles.f90). The format for the file is as follows.

**Water Depth and Wave Frequencies**

| Line | Format | Variable | Description |
| :--- | :--- | :--- | :--- |
| 1 | free | - | Comment or blank line (ignored) |
| 2 | free | - | Comment or blank line (ignored) |
| 3 | `14x,f30.15` | `H` | Water depth (m). <0 for infinite depth, >0 for finite depth |
| 4 | free | - | Comment or blank line |
| 5 | free | - | Comment or blank line |
| 6 | `27x,i16` | `SYB0` | Switch: 0 or 1. 1 to include zero and infinite frequency added mass calculations; 0 to exclude them.  |
| 7 | `25x,i16` | `INFT` | Input frequency type: 1 for deepwater wave number, 2 for finite depth wave number, 3 for wave frequency, 4 for wave period and 5 for wave length. 25 spaces before the integer |
| 8 | `25x,i17` | `OUFT` | Output frequency type; same options as input frequency type |
| 9 | `26x,i16` | `NPET` | Number of frequencies |
| if `NPET > 0` | space-separated list of length `NPET` (single line) | `WVNB` | Discrete wave frequencies, periods, numbers or lengths, depending on `Input_frequency_type` (see [example](https://github.com/YingyiLiu/HAMS/blob/master/CertTest/Moonpool/Input/ControlFile.in))
| if `NPET < 0` | `27x,f30.15` | `WK1` | Minimum frequency
| if `NPET < 0` | `19x,f30.15` | `DWK` | Frequency step

**Body Information**

As shown in the table below, the first two lines after the last value from the previous section are ignored.

| Line | Format | Variable | Description |
| :--- | :--- | :--- | :--- |
| next | free | - | Comment or blank line (ignored) |
| next | free | - | Comment or blank line (ignored) |
| next | `A100` | `BODY_CHECK` | String indicating whether multiple bodies are used. Include `"multi"` if multiple bodies. |
| if `BODY_CHECK` contains "multi": | `24x,i16` | `NBODY` | Number of bodies; must be ≥1 |
| next `NBODY` lines | `26x,4f12.3` |`LCS_MULTI` | origin coordinates of the local coordinate system (LCS), per body mesh

**Wave Headings**

Wave incident angles with respect to the x-direction. All wave heading input values are expressed in degrees.

As shown in the table below, the first three lines after the last value from the previous section are ignored.

| Line | Format | Variable | Description |
| :--- | :--- | :--- | :--- |
| next | free | - | Comment or blank line (ignored) |
| next | free | - | Comment or blank line (ignored) |
| next | free | - | Comment or blank line (ignored) |
| next | `23x,i16` | `NBETA` | Number of wave headings; a negative sign indicates values are generated from a minimum and a step size rather than listed explicitly
| if `NBETA > 0` | space-separated list of length `NBETA` | `WVHD` | `NBETA` wave headings |
| if `NBETA < 0` | `20x,f30.15` | `BETA1` | minimum heading
| if `NBETA < 0` | `17x,f30.15` | `DBETA` | heading step

**Rotation Centers**

As shown in the table below, the first two lines after the last value from the previous section are ignored.

| Line | Format | Variable | Description |
| :--- | :--- | :--- | :--- |
| next | free | - | Comment or blank line (ignored) |
| next | free | - | Comment or blank line (ignored) |
| if `NBODY=1` | `28x,3f12.3` | `XR(1:3)` | Rotation center coordinates |
| if `NBODY>1` | `28x,3f12.3` | `XR_MULTI(NBODY,1:3)` | Each body's rotation center coordinates; one set of coordinates per line |

**Solver Options**

This section begines on the next line after the previous one. 

| Line | Format | Variable | Description |
| :--- | :--- | :--- | :--- |
| next | `26x,f30.15` | `REFL` | Reference body length |
| next | `26x,i16` | `ISOL` | Wave diffraction solution. Switch: 1 or 2. 1 if incident and scattering potentials are combined; 2 if only scattering.  |
| next | `23x,i16` | `IRSP` | Irregular frequencies: 1 for to remove them, 0 to keep them |
| next | `23x,i16` | `NTHREAD` | Number of OpenMP threads; setting it to 1 removes parallelization |
| (optional) | keyword, then free-format integer | `ISOLV` | Linear solver: 1 = direct LU (default), 2 = GMRES with H-matrix acceleration. Line is detected by the `Solver_type` keyword; the first integer after it is read and a trailing `# comment` is allowed. If the line is absent, ISOLV defaults to 1; if the value cannot be read, a warning is printed and ISOLV defaults to 1. See [iterative_solvers.md](iterative_solvers.md) for full details. |
| (optional) | keyword, then free-format integer | `LIDON` | Damping lid: 0 = none (default), 1 = read `DampingLid.pnl`. The line sits after `Solver_type`, or in its place when `Solver_type` is absent, and is read the same way: the `Damping_lid` keyword, the first integer after it, a trailing `# comment` allowed. When both lines are present, `Damping_lid` must follow `Solver_type`; in the other order the field-point section is misread and the run fails with a Fortran read error. A value other than 0 or 1, or no readable integer, stops the run (`Damping_lid must be 0 or 1`). See [DampingLid.pnl](#dampinglidpnl). |

**Pressure and Elevation**

As shown in the table below, the first two lines after the last value from the previous section are ignored.

| Line | Format | Variable | Description |
| :--- | :--- | :--- | :--- |
| next | free | - | Comment or blank line (ignored) |
| next | free | - | Comment or blank line (ignored) |
| next | `27x,i16` | `NFP` | Number of field points to output field pressure and elevation. |
| next `NFP` lines | free | `XFP(I,1:3)` | Coordinates of field points (1 set per line) |

**Output Files**

Every simulation outputs the hydrodynamic pressure and free-surface elevation on the bodies for the radiation and diffraction problems. These quantities are written in WAMIT format to the file `PressureElevation.6p` by default.

If `separate_wamit_diffraction_and_radiation_files = true`, the radiation and diffraction results are written to separate files instead:
- Radiation: `PressureElevationRadiation_<1_to_number_bodies>.6p` (one file per body)
- Diffraction: `PressureElevationDiffraction.6p`

| Line | Format | Variable | Description |
| :--- | :--- | :--- | :--- |
| next | free | - | Comment or blank line (ignored) |
| next | free | - | Comment or blank line (ignored) |
| next | `50x,i16` | `separate_wamit_diffraction_and_radiation_files` | Switch: 0 or 1. 0 for false, 1 for true. |


## Hydrostatic.in

Hydrostatic and damping properties are specified in the `Hydrostatic[_i].in` files. They are read by the subroutine `ReadHydroStaticMulti` in [HydroStaticMulti.f90](../src/HydroStaticMulti.f90). There is one file per body.


| Line | Format | Variable | Description |
| :--- | :--- | :--- | :--- |
| 1 | free | - | Comment or blank line (ignored) |
| 2 | `6(2x,E12.5)` | `XG_MULTI(BODY_N,1:3)` | Center of gravity coordinates (x, y, z) |
| 3 | free | - | Comment or blank line (ignored) |
| 4-9 | `6(2x,E12.5)` | `MATX_MULTI(BODY_N,1:6,1:6)` | 6×6 mass matrix, one row per line |
| 10 | free | - | Comment or blank line (ignored) |
| 11-16 | `6(2x,E12.5)` | `BLNR_MULTI(BODY_N,1:6,1:6)` | 6×6 external linear damping matrix, one row per line |
| 17 | free | - | Comment or blank line (ignored) |
| 18-23 | `6(2x,E12.5)` | `BQDR_MULTI(BODY_N,1:6,1:6)` | 6×6 external quadratic damping matrix, one row per line |
| 24 | free | - | Comment or blank line (ignored) |
| 25-30 | `6(2x,E12.5)` | `CRS_MULTI(NBODY,1:6,1:6)` | 6×6 hydrostatic restoring matrix, one row per line |
| 31 | free | - | Comment or blank line (ignored) |
| 32-37 | `6(2x,E12.5)` | `KSTF_MULTI(NBODY,1:6,1:6)` | 6×6 external restoring matrix, one row per line |


## HullMesh.pnl

Node and element data on the surface of each body are specified in the `HullMesh[_i].pnl` files. There is one file per body.

The files are read in [HAMS_Prog.f90](../src/HAMS_Prog.f90) and in the subroutine `ReadBodyMeshMulti` [ReadPanelMeshMulti.f90](../src/ReadPanelMeshMulti.f90). The expected file structure is described below.

A few notes:
- Numbers in each line can be separated by whitespaces but **not commas**.
- Any trailing text on each line must be preceded by `!`.
- The number of vertex indices must match `NTND` exactly.

| Line(s) | Format | Variable | Description |
| :------ | :------ | :-------- | :----------- |
| 1–3 | free | — | Comment or blank lines (ignored) |
| 4 | free | `NELEM, NTND, ISX, ISY` | Number of panels, number of nodes, X-symmetry flag, Y-symmetry flag |
| 5-6 | free | — | Comment or blank lines (ignored) |
| next `NTND` | free | `M` and `XYZ_LOCAL_MULTI(bodyID,nodeID,1:3)` | Node coordinates. Each line contains: `node_ID  x_local  y_local  z_local` |
| next 3 lines | free | — | Comment or blank (ignored) |
| next `NTND` lines | free | `M`, `NCN_MULTI(bodyID,nodeID)` and `NCON_MULTI(bodyID,nodeID,1:numvertices)`| Panel definitions. Each line contains: panel number, number of vertices and (vertex1, vertex2, ...)


## WaterPlaneMesh.pnl

Node and element data on the inner water plane for multiple bodies are specified in the `WaterPlaneMesh[_i].pnl` files. There is one file per body.

The files are read in [HAMS_Prog.f90](../src/HAMS_Prog.f90) and in the subroutine `ReadWTPLMeshMulti` in [ReadPanelMeshMulti.f90](../src/ReadPanelMeshMulti.f90). The file structure is exactly the same as the one for the HullMesh.pnl files, which is described above.

The node coordinates of the waterplane mesh are stored in `iXYZ_LOCAL_MULTI(bodyID,nodeID,1:3)`. The number of vertices in each panel is stored in `iNCN_MULTI(bodyID,nodeID)` and the vertices (vertex1, vertex2, ...) are stored in `iNCON_MULTI(bodyID,nodeID,1:numvertices)`.


## DampingLid.pnl

Read only when `ControlFile.in` has `Damping_lid 1`. GMRES (`Solver_type 2`) is not supported with a lid. `If_remove_irr_freq 1` with a lid is experimental for a single body (see below), and stops for several bodies.

The lid is a list of triangular or quadrilateral panels on the free surface z = 0, with the same file structure as `WaterPlaneMesh.pnl` (three header lines, then `NELEM NTND ISX ISY`, the nodes, the panels). There is one file for any number of bodies, because a lid belongs to no body. Its nodes are in global coordinates, not in the local frame of a body. Holes and separate patches are allowed, since the lid is only a list of panels. The vertex order of a panel does not matter, because the lid has no normal.

After the panels (after the `# End Definition of Node Relations` line) the file has a damping block, in the layout of the WAMIT `.dmp` file:

```
    # Start Definition of Damping   ! number of dampers, then: first panel  last panel  epsilon (m/s)
         2
         1        40    0.2000000000
        41        80    0.1000000000
    # End Definition of Damping
```

The lines are read free-format, and a trailing `! comment` is allowed. The second line is the number of dampers, NDAMPERS. Each damper line gives IP1 and IP2, the first and last panel number of `DampingLid.pnl` (1 to the number of lid panels), and EDAMPER, epsilon in m/s. The ranges must cover every lid panel exactly once. With `ISX` or `ISY` the images of a panel take the epsilon of the panel. epsilon is constant over frequency, epsilon = 0 is allowed, and there is no default. One log line per damper is printed at the start of the run: its range, number of panels and epsilon.

The lid may touch a body or lie away from every body. The file is read by `READ_LID_MESH` in [DampingLid.f90](../src/DampingLid.f90). Unlike the `WaterPlaneMesh.pnl` read in `HAMS_Prog.f90`, it does not overwrite `ISX` and `ISY`: it compares them with the hull's.

What the lid does: it is the dissipative lid of WAMIT (WAMIT User Manual v7.5, Section 12.8 and eqs. 15.91 to 15.93). On the lid the free surface carries a pressure p/rho = epsilon dphi/dz, proportional to the vertical velocity, which takes work from the fluid and so absorbs energy. With the time factor e^{-i omega t} the free-surface condition on the lid becomes dphi/dz = nu g phi/(g - i omega epsilon), nu = omega**2/g. epsilon = 0 gives the open free surface (no lid), and a large epsilon a rigid lid (dphi/dz = 0). Derived from the condition (spec, "Formulation"): the power absorbed is 1/2 rho epsilon abs(dphi/dz)**2 per unit area, never negative. The formulation is in [developer_documentation.md](developer_documentation.md#damping-lid).

This lid replaces the lid that held phi = 0 on the panels (design `docs/superpowers/specs/2026-11-05-damping-lid-design.md`, now superseded by `docs/superpowers/specs/2026-11-06-epsilon-damping-lid-design.md`). That lid reflected the wave and absorbed no energy. The `DampingLid.pnl` files written for it have no damping block and now stop the run (`no damping block`).

Choosing epsilon. The WAMIT manual gives no recommended value. The values below are starting points only: epsilon depends on the geometry and the flow that the lid stands for, and must be calibrated against model tests or CFD.

| Source | epsilon (m/s) | Use |
| :--- | ---: | :--- |
| WAMIT v7.5 test cases test02, test07a, test17c, test22c | 0.4; 1.0; 0.1 and 0.4; 1.0 | examples of the input, not calibrated values |
| Gourlay (2019) | 0.2 | gap between side-by-side vessels |

Some papers give the lid damping in the form of Chen, a dimensionless epsilon_Chen. For small damping, Gourlay (2019) gives epsilon ~ g epsilon_Chen / omega. This is for information only: it depends on the frequency, while the epsilon of `DampingLid.pnl` is constant over frequency.

The run stops with a message and exit code 1 when:

| Condition | Message contains |
| :--- | :--- |
| `Damping_lid 1` and no `DampingLid.pnl` | `DampingLid.pnl input file missing` |
| a node with abs(z) > 1e-10 | `lid node is not on z = 0` |
| `ISX` or `ISY` of the lid differ from the hull's | `lid symmetry differs from the hull` |
| `Solver_type 2` | `damping lid needs Solver_type 1` |
| a lid panel centroid lies inside a hull | `lid panel inside a hull` |
| more than one body and `If_remove_irr_freq 1` | `damping lid with several bodies needs If_remove_irr_freq 0` |
| `Damping_lid` other than 0 or 1, or no integer after it | `Damping_lid must be 0 or 1` |
| a line of `DampingLid.pnl` that cannot be read (header, counts, node, panel) | `cannot read DampingLid.pnl`, with the kind of line and its number |
| a lid panel with other than 3 or 4 vertices | `a lid panel in DampingLid.pnl must have 3 or 4 vertices` |
| a vertex number outside 1 to the number of nodes | `a vertex of a lid panel in DampingLid.pnl is not a node number` |
| a lid panel with a repeated vertex | `a lid panel in DampingLid.pnl has a repeated vertex` |
| a lid panel of zero area (area below 1e-12 times the square of its size) | `a lid panel in DampingLid.pnl has zero area` |
| no `# Start Definition of Damping` line after the panels | `no damping block` |
| number of dampers below 1 or above the number of lid panels, or unreadable | `number of dampers` (or `cannot read the damping`), with the line number |
| a damper range with first > last, or outside 1 to the number of lid panels | `damper panel range`, with the line number |
| a lid panel in two dampers | `lid panel in two dampers`, with the panel number |
| a lid panel in no damper | `lid panel in no damper`, with the panel number |
| a negative epsilon | `negative damping`, with the line number |
| a damper line that cannot be read | `cannot read the damping`, with the line number |
| a zero frequency in the frequency list | `zero frequency with a damping lid`, with its number in the list |

The zero frequency that `0_inf_frequency_limits 1` adds to the run is not in the list and does not stop it, since the infinite-frequency limit cannot be given without it (see Zero frequency below).

The inside-hull check sums the hull's double-layer Rankine integrals and their images in z = 0 at each lid centroid. This is a solid angle: 4 pi in magnitude inside a closed hull and 0 outside, and the run stops above 2 pi. A lid inside a hull is what the water plane mesh with `If_remove_irr_freq` is for. Only the centroids are tested, so a panel that straddles the waterline with its centroid outside the hull passes.

Finite depth: every integral over a lid panel takes the finite-depth part as a smooth remainder sampled at R/h of at least 0.03, because `FINGREEN3D` is unreliable at small R/h. Measured on z = 0 against a direct sum of the eigenfunction series ([README](../test/data/damping_lid/README.md), "Remainder"): `FINGREEN3D` is off by a near constant 7.9e-6 between 0.0005 h and 0.05 h, and below 0.003 h it scatters, to 1.1e-2 at 0.001 h and nu = 1. Lid-to-lid pairs are protected by the 0.03 h floor. Derived: the hull panels are not, since the lid rows and the field output evaluate them with `FINGREEN3D` at the true distance, so lid centroids or field points within about 0.01 h of a hull panel centroid get the scatter. See [developer_documentation.md](developer_documentation.md#damping-lid).

`If_remove_irr_freq 1` with a single-body lid is experimental. The system then has more rows than unknowns and is solved by least squares (`ZGELS`). Measured for `lid_ann_160` at epsilon 0.2 m/s: the worst change of A11, A33, B11, B33 and abs(X1) against the same run without `If_remove_irr_freq` is 0.32% (B11). With the phi = 0 lid it was 1.1%, and the cause was not found; that study was not repeated with the epsilon lid ([README](../test/data/damping_lid/README.md), "Symmetry and IRSP").

Zero frequency: a zero in the frequency list stops a run with the lid on, as in WAMIT (WAMIT User Manual v7.5, Section 12.8). The zero frequency that `0_inf_frequency_limits 1` adds is exempt: there, and at the infinite frequency, the lid is skipped and the hull is solved alone. Derived: at omega = 0 the lid term alpha = i omega epsilon nu/(g - i omega epsilon) is 0, so the lid has no effect there.

Output with a lid. No file format changes. A, B, the excitation and the RAOs are integrated over the hull only, as in WAMIT (WAMIT User Manual v7.5, Section 12.8). Field points:

- The field potential includes the lid term, a single layer of strength alpha_l (phi_l + phi_I,l) on each lid panel, with phi_I,l the incident potential in diffraction only.
- On z = 0 inside or on the edge of a lid panel of epsilon > 0, the elevation is eta = i omega phi/(g - i omega epsilon), with phi the total potential; elsewhere it is i omega phi/g, as without a lid. The pressure in the fluid has its usual form.
- Field points at z = 0 get radiation elevation rows in `PressureElevation.6p` only when a lid is on; they are not written without one.

What changes in the results, derived (spec, "Limits that follow"): B is the power radiated to the far field plus the power absorbed on the lid. B can still fall when the lid is added, because the lid can cut the far-field radiation by more than it absorbs (measured: B11 of `lid_ann_160` falls from 206.25 without the lid to 199.71 at 0.2 m/s, README "Single-body solve with a lid"). The Haskind relation over the hull does not hold for epsilon > 0, because the lid absorbs: do not compute the excitation from the radiation potentials when a lid is on (measured heave residual 4.6e-2 at 0.2 m/s, against 2.1e-3 without the lid).

The number of lid panels is printed next to the hull panel counts and written to `ErrorCheck.txt`, only when the lid is on, so the output of a run without a lid does not change.

## gen_mod (Optional)

Generalized mode files allow users to override one or more of the standard six rigid-body normals (surge, sway, heave, roll, pitch, yaw) with a custom shape function for any body. This enables the use of generalized modes within the BEM framework — for example, replacing a rigid-body DOF with a hinge rotation or a flexible deformation mode. There is at most one file per body per DOF. These files are optional; if no file is present for a given body and DOF, the standard rigid-body normal is retained. This approach is utilized since often times, not all 6 DOFs of a body are utilized.

The files are read in the subroutine `CalNormals` in [ReadPanelMesh.f90](../src/ReadPanelMesh.f90) for single-body cases and `CalNormalsMulti` in [ReadPanelMeshMulti.f90](../src/ReadPanelMeshMulti.f90) for multi-body cases. After the standard translational and rotational normals are computed, the program scans the input directory for generalized mode files. If a file is found, the corresponding column in `DXYZ_P` (single-body) or `DXYZ_MULTI_P` (multi-body) is overwritten with the values from the file.

The files must be placed in the same input directory as the other input files and follow the naming convention:

```
gen_mod_{body_number}_{dof}.txt
```

where `body_number` is the body index (always `1` for single-body cases; `1, 2, ...` for multi-body cases) and `dof` is the degree of freedom to override (`1` = surge, `2` = sway, `3` = heave, `4` = roll, `5` = pitch, `6` = yaw). For example, `gen_mod_1_3.txt` overrides heave normals for body 1 and `gen_mod_2_5.txt` overrides pitch normals for body 2.

### What the values mean

The values are the normal of that DOF, the term that DOF puts in the radiation boundary condition. They are used exactly as given: they are not rescaled, normalised or recomputed, so a value is not limited to the range of a unit normal.

They affect the places that use the DOF's normal:
- the radiation right side of that DOF (`RBC_RIGHT`, `RBC_IRR` and their multi-body versions),
- the excitation force of that DOF, and the added mass and damping in its row (`EFORCE`, `RFORCE` and their multi-body versions),
- that mode's own term in pressure and elevation at field points (`CalPotential`, `CalPotentialMulti`).

Through these, the mode's potential also sets the column of that DOF in added mass and damping. All other entries of added mass, damping and excitation do not depend on the mode, and equal those of a run without it.

The values do not change the geometry. Code that needs the panel's direction reads the geometric normal, `GEONRM_P` or `GEONRM_MULTI_COMB`, stored before the override. This covers dG/dn on the left side, the incident wave's normal derivative in the diffraction right side, and the kernel, incident term and on-body test at field points.

The mode replaces the rigid-body motion of that DOF, so the rigid-body results of that DOF are given up. The mass and restoring of the DOF are whatever `Hydrostatic.in` holds in that row and column, so motions are meaningful only if these are set for the mode. With a symmetry plane, the sign factors in the right side and the force integrals differ by DOF, so the mode must have that DOF's parity. No case in the repository tests a mode with a symmetry plane.

A few notes:
- The number of lines must match the number of panels (`NELEM`) for the corresponding body.
- Panel IDs must match the panel numbering in the corresponding `HullMesh` file.
- Numbers in each line are separated by whitespace and read with free format.

| Column | Format | Variable | Description |
| :----- | :----- | :------- | :---------- |
| 1 | free | `GEN_PID` | Panel ID |
| 2 | free | `GEN_NORMAL` | Generalized normal value for that panel |

### Generating Generalized Mode Files

A Python helper script (`generalized_modes.py`) is provided to generate these files. The script reads the hull mesh file (`HullMesh[_i].pnl`), computes panel centroid coordinates, reads the LCS from `ControlFile.in` for the corresponding body and transforms the centroids to global coordinates. The generalized normal at each panel centroid is then computed as `n_j = nx * X + ny * Y + nz * Z`, where `X`, `Y`, `Z` are the global centroid coordinates and `nx`, `ny`, `nz` are the user-defined shape function components. The body number is extracted automatically from the mesh filename (e.g., `HullMesh_2.pnl` → body 2).