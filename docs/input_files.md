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

## SecondOrder.in (Optional)

`SecondOrder.in` switches on the mean wave drift calculation. The file is optional: if it is absent, or `mean_drift = 0`, HAMS-MREL runs exactly as before. It is read from the input directory, right after `ControlFile.in`, and the run stops with a message if an option is unknown, has an invalid value, or conflicts with the first-order set-up.

Each line is `key = value`. Text after `!` is a comment, blank lines are ignored, keys are case-insensitive and may appear in any order. Keys that are left out take their default.

| Key | Values | Default | Description |
| :-- | :----- | :------ | :---------- |
| `mean_drift` | `0`, `1`, `2` | `0` | 0 is off, 1 is unidirectional, 2 is all heading pairs. |
| `drift_methods` | `far`, `near`, or both separated by a blank or comma | `far` | Far-field (momentum) and near-field (pressure integration) methods. |
| `velocity_method` | `differencing`, `source` | `differencing` | How the velocity on the body surface is obtained, for the near-field method. `differencing` takes the tangential part from the potential on the neighbouring panels, and the normal part from the known normal derivative. `source` solves for a source distribution with the same normal derivative, and takes the velocity from its integral: one more dense factorization and one more set of Green function evaluations per frequency. The potential itself stays that of the direct solve. The solve keeps the factorization and the gradient integrals of every symmetry block in memory, 64 nb^2 ns bytes with nb panels in the base mesh and ns symmetry blocks: 236 MB at 1920 panels and 1.2 GB at 4320 with no symmetry, about as much again as the first-order arrays. Irregular frequencies are not removed for this method, since `IRSP = 1` is rejected, and the run stops if a matrix is exactly singular. |
| `waterline_method` | `extrapolate`, `field` | `extrapolate` | How the potential at the waterline is obtained, for the near-field method. `extrapolate` takes it from the adjacent panel and its gradient. `field` evaluates the boundary integral formula at a point close to the waterline, out from the body and down. The point is in the fluid if the wall at the waterline widens downward by no more than 45 degrees from the vertical, and if no other body is closer to it than `waterline_offset` times the segment length, so keep the gap between bodies larger than that. On a wider wall the run stops with a message. It costs one set of Green function evaluations per waterline segment and frequency. It converges more slowly than `extrapolate`, and in finite depth it can jump with the mesh or the offset: see `test/data/drift_cyl/README.md`. |
| `waterline_offset` | positive real | `0.01` | With `waterline_method = field`: how far the point lies from the waterline, out from the body and down, each as a fraction of the waterline segment length. |
| `motion_state` | `free`, `fixed` | `free` | Whether the body responds to the waves or is held fixed. |
| `n_azimuth` | integer, at least 1 | `360` | Number of azimuth points in the far-field integration. |
| `far_tolerance` | real | `1d-4` | Tolerance of the far-field integration. |
| `write_terms` | `0`, `1`, `2` | `1` | 0 writes nothing extra, 1 writes the term file, 2 also writes panel diagnostics. |
| `wl_z_tolerance` | real | `1d-6` | Tolerance for panels at the waterline, relative to the body size. |

Example:

```
! Mean drift with both methods
mean_drift     = 1
drift_methods  = far near
n_azimuth      = 360
```

The options are checked against `ControlFile.in` and the mesh files. A run is rejected when:
- `ISOL = 2`, because the incident potential is needed.
- The water depth `H` is above 500 m. Set `H < 0` for infinite depth instead, because the incident-wave velocity routines treat depths above 500 m differently.
- `drift_methods` includes `near` and any generalized mode file (`gen_mod_{body}_{dof}.txt`) is present.
- `velocity_method = source` with `ISOLV = 2` (GMRES), or with `IRSP = 1` (irregular frequency removal).
