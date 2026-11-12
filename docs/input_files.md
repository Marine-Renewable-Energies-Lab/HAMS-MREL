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


## CSMesh.pnl (only with `cs_shape = file`)

With the control-surface mean drift and `cs_shape = file`, the control surface C of body `b` is read from `CSMesh_{b}.pnl` in the input directory. There is one file per body, numbered `1 … numberofbodies` without leading zeros, also for a single body (`CSMesh_1.pnl`). If a file is missing, the run stops with `cs_shape = file needs CSMesh_<b>.pnl in the input directory.`

The file has exactly the structure of `HullMesh.pnl` (see above), with these rules:
- `ISX = ISY = 0`. C is always given in full.
- Panels have 3 or 4 vertices.
- The coordinates are global, not body-local.
- **The panel normals must point towards the body.** HAMS takes the normal of a quadrilateral as (x4 - x2) x (x3 - x1) and of a triangle as (x1 - x2) x (x3 - x2), so the node order sets the direction. The run checks the sign of the sum over the panels of n . (x_c - x_b) A, with x_b the centre of the waterplane bounding box of the body; if it is not negative, the run stops with `CSMesh_<b>.pnl: the panel normals must point towards the body.`
- **Containment is the user's job.** C must enclose the body and no other body. The run does not test this for a file. It only checks that no point of C lies within one panel size of a hull panel centroid.
- C must cut the free surface z = 0 in a closed line. Gaps in that waterline are reported as warnings.

`write_control_surface` writes the generated C in this format, so a generated file can be edited and read back.

## FSMesh.pnl (optional, with the control-surface method)

The vertical force and the roll and pitch moments of the control-surface method also need the free surface F between the hull waterline and Γc, the line where C meets z = 0 (Chen 2007, §3.1). By default F of body `b` is generated: quadrilateral strips from each hull waterline segment out to Γc, along rays from the centre of the waterplane bounding box, with `cs_fs_rows` rows graded towards the hull (the inner row is 1/4 of the outer). The rays only cover F when the waterline is star-shaped about that centre. If it is not, the run goes on without F for that body: its F_z, M_x and M_y are not written, and the log says `Warning, body <b>: the waterline is not star-shaped about its centre; give FSMesh_<b>.pnl. No free surface F, so F_z, M_x, M_y of body <b> are not computed.` The other bodies keep their F and their rows I = 3, 4, 5.

`FSMesh_{b}.pnl` in the input directory replaces the generated F of body `b`. It has the structure of `CSMesh_{b}.pnl` (`ISX = ISY = 0`, 3 or 4 vertices, global coordinates), with these rules. A given file is a request for F, so a file that cannot be read or breaks the first rule stops the run with its message:
- All nodes lie in z = 0. The run stops with `FSMesh_<b>.pnl: F must lie in the free surface z = 0.` otherwise.
- The node order does not matter: the normal of F is always (0, 0, −1), into the fluid.
- F must fill the free surface between the hull waterline and Γc. The run does not test this.

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
| `drift_methods` | any of `far`, `near`, `control`, separated by blanks or commas | `far` | Far-field (momentum), near-field (pressure integration) and control-surface methods. For the vertical force on a hull with sharp edges (a right-angled bottom edge), prefer `control`: on a truncated cylinder (radius 10 m, draft 10 m, depth 100 m, 7168 panels, 0.70 to 1.10 rad/s), control F_z is within 1.44% of a semi-analytical reference, where near-field is 12% to 17% off on the fixed hull and 23% to 506% off on the free hull, with the wrong sign at 0.70 rad/s (see [the accuracy of the vertical loads](second_order_validation.md#vertical-force-and-roll-and-pitch-moments-f_z-m_x-and-m_y)). |
| `velocity_method` | `differencing`, `source` | `differencing` | How the velocity on the body surface is obtained, for the near-field method. `differencing` takes the tangential part from the potential on the neighbouring panels, and the normal part from the known normal derivative. `source` solves for a source distribution with the same normal derivative, and takes the velocity from its integral: one more dense factorization and one more set of Green function evaluations per frequency. The potential itself stays that of the direct solve. The solve keeps the factorization and the gradient integrals of every symmetry block in memory, 64 nb^2 ns bytes with nb panels in the base mesh and ns symmetry blocks: 236 MB at 1920 panels and 1.2 GB at 4320 with no symmetry, about as much again as the first-order arrays. Irregular frequencies are not removed for this method, since `IRSP = 1` is rejected, and the run stops if a matrix is exactly singular. |
| `waterline_method` | `extrapolate`, `field` | `extrapolate` | How the potential at the waterline is obtained, for the near-field method. `extrapolate` takes it from the adjacent panel and its gradient. `field` evaluates the boundary integral formula at a point close to the waterline, out from the body and down. The point is in the fluid if the wall at the waterline widens downward by no more than 45 degrees from the vertical, and if no other body is closer to it than `waterline_offset` times the segment length, so keep the gap between bodies larger than that. On a wider wall the run stops with a message. It costs one set of Green function evaluations per waterline segment and frequency. It converges more slowly than `extrapolate`, and in finite depth it can jump with the mesh or the offset: see `test/data/drift_cyl/README.md`. |
| `waterline_offset` | positive real | `0.01` | With `waterline_method = field`: how far the point lies from the waterline, out from the body and down, each as a fraction of the waterline segment length. |
| `motion_state` | `free`, `fixed` | `free` | Whether the body responds to the waves or is held fixed. |
| `n_azimuth` | integer, at least 1 | `360` | Number of azimuth points in the far-field integration. |
| `far_tolerance` | real | `1d-4` | Tolerance of the far-field integration. |
| `write_terms` | `0`, `1`, `2` | `1` | 0 writes nothing extra, 1 writes the term file, 2 also writes panel diagnostics. |
| `wl_z_tolerance` | real | `1d-6` | Tolerance for panels at the waterline, relative to the body size. |
| `cs_velocity` | `green`, `source` | `green` | With `control` or `qtf_b2 = control`: how the velocity on the control surface is obtained. `green` differentiates the potential by finite differences of the Green function. `source` solves for a source distribution and takes the velocity from its integral. In finite depth above ν·h = 50 (ν = ω²/g), the Green function of `green` still uses the series acceleration (Limes) that makes field values jump at small horizontal distances, so its values on C may scatter, and the log warns at each such frequency; there, check against `far`, or use deep water (`H < 0`) where kh is large. |
| `cs_shape` | `cylinder`, `box`, `file` | `cylinder` | With `control` or `qtf_b2 = control`: the shape of the control surface. `cylinder` and `box` are generated around each body. `file` reads `CSMesh_{b}.pnl` for each body. |
| `cs_margin` | real, at least 0 | `0` | With `control` or `qtf_b2 = control`: the horizontal gap in metres between the hull and the control surface. 0 means 4 times the median hull panel size. |
| `cs_bottom_margin` | real, at least 0 | `0` | With `control` or `qtf_b2 = control`: the gap in metres below the deepest hull point. 0 means 4 times the median hull panel size. |
| `cs_panel_size` | real, at least 0 | `0` | With `control` or `qtf_b2 = control`: the panel size of the control surface in metres. 0 means the median hull panel size. |
| `cs_fd_step` | real, above 0 and at most 0.5 | `0.02` | With `cs_velocity = green`, and for ∇ψ on C and F with `qtf_b2 = control`: the finite-difference step, as a fraction of the local control-surface panel size. Above about 0.3, the centroids of the top row of C (a quarter panel below z = 0) come within one step of z = 0, and their dφ/dz is then taken from the free-surface condition, which is off by order k z below z = 0 [derived]. Keep the default. |
| `cs_fs_rows` | integer, at least 0 | `0` | Number of rows of the free-surface patch F between the waterline and the control surface. 0 means the nearest integer to (distance of Γc − distance of the waterline)/`cs_panel_size`, at least 2; the distances are mean distances from the centre of the waterplane bounding box. Not used when `FSMesh_{b}.pnl` is given. |
| `qtf` | `0`, `1`, `2`, `3` | `0` | Quadratic part of the second-order QTFs: 0 off, 1 difference frequency (`QTF_Quadratic.10d`), 2 sum frequency (`QTF_Quadratic.10s`), 3 both. Needs `mean_drift = 1` or `2` and `near` in `drift_methods`. Heading pairs follow `mean_drift`: β1 = β2 for 1, every pair for 2. The first-order state of every frequency and heading is kept in memory, 128 bytes per panel each; the log gives the size. See [QTFs](#qtfs-quadratic-part). |
| `qtf_terms` | `quadratic`, `body`, `full` | `quadratic` | Which parts of the QTF. `quadratic`: the quadratic part only (products of first-order quantities). `body`: also the body forcing of the second-order potential, by the indirect (Haskind) method with the radiation potentials ψ at the difference frequencies, written to `QTF_Indirect.11d`. Difference QTF only; the free-surface integral of the second-order potential is left out. `full`: `body` plus the free-surface integral (Phase 6a, finite depth only, H > 0), so that `QTF_Indirect.11d` is the complete indirect-method QTF; it implies everything `body` does, and `qtf_b2`, `qtf_edge_weights` and `qtf_edge_angle` apply to it as to `body`. **Status:** the keys are read and checked, but the free-surface terms are not yet computed, so `full` gives the output of `body` for now. Needs `qtf = 1` or `3`: with `qtf = 2` the run stops, and with `qtf = 0` `body` and `full` do nothing. See [QTFs, potential part](#qtfs-potential-part). |
| `qtf_max_dw` | real, at least 0 | `0` | Largest \|ω_i − ω_j\| in rad/s written to `QTF_Quadratic.10d` and `QTF_Indirect.11d`. 0 writes every pair. With `qtf_terms = body` it also limits the difference frequencies at which ψ is solved. The sum-frequency file always has every pair. Choose it by a convergence study of the slow-drift response, widening the bandwidth until the result no longer changes, as Zhang et al. (2026) do for the full difference-frequency QTF of a floating wind turbine. The validation used 0.3 rad/s. |
| `qtf_b2` | `control`, `body` | `control` | With `qtf_terms = body`: where the surface piece of the second-derivative term B2 is integrated. `control` moves it from the hull onto the control surface C and the free surface F around each body by a Gauss identity (Cong et al. 2020; Chen 2007, App. C); the rest of B2 stays on the hull. `body` integrates it on the hull, which is inaccurate next to sharp edges ([validation](second_order_validation.md#qtf-potential-part)). `control` builds C and F with the `cs_*` keys, also when `drift_methods` has no `control`, and needs a C with a bottom. Ignored, with a note in the log, without `qtf_terms = body` or `full` and `qtf = 1` or `3`. |
| `qtf_edge_weights` | `off`, `on` | `off` | With `qtf_b2 = control`: product weights on the hull panels next to a sharp edge, for the hull term of B2 ([SharpEdge.f90](../src/SharpEdge.f90), velocity ~ ρ^p with p = π/α_f − 1, α_f the fluid wedge angle). Off by default: with HAMS's panel velocities it moved B2 further from the exact reference ([validation](second_order_validation.md#tried-and-not-kept)). Ignored, with a note, without `qtf_b2 = control`. |
| `qtf_edge_angle` | real, above 0 and below 180 | `45` | With `qtf_edge_weights = on`: a mesh edge is sharp when the normals of its two panels differ by more than this many degrees. |
| `qtf_fs_radius` | real, at least 0 | `0` | With `qtf_terms = full`: the radius R in metres of the circle around all bodies that splits the free surface into an inner mesh (waterline to R) and an outer region (r > R, partial waves). 0 sets the value automatically. Ignored, with a note in the log, without `qtf_terms = full` and `qtf = 1` or `3`. |
| `qtf_fs_rows` | integer, at least 0 | `0` | With `qtf_terms = full`: the number of radial rows of the generated inner mesh. 0 sets the value automatically. Ignored, with a note, as `qtf_fs_radius`. |
| `qtf_fs_modes` | integer, at least 0 | `0` | With `qtf_terms = full`: the number of angular modes M of the outer region. 0 sets the value automatically. Ignored, with a note, as `qtf_fs_radius`. |
| `qtf_fs_evanescent` | integer, at least 0 | `0` | With `qtf_terms = full`: the number N_ev of evanescent vertical modes of the outer region, besides the propagating mode. `0` sets N_ev = max(8, ceil(9 h / (pi (R - r_wl)))), r_wl the largest distance from the centre of the circle to a waterline node, so that the slowest mode kept decays by about e^-9 from the waterline to R [derived: k_n is about n pi / h]. The log gives the value used. Ignored, with a note, as `qtf_fs_radius`. |
| `qtf_field` | `0`, `1` | `0` | 1 writes `QTF_FieldPartial.txt`: the second-order pressure and free-surface elevation at the field points of `ControlFile.in`, for the pairs of `QTF_Quadratic.10d` and `.10s`. Needs `qtf` 1, 2 or 3 and at least one field point. **Partial:** the second-order diffracted and radiated waves are omitted, and near the bodies they are not small. See [Partial second-order field](#partial-second-order-field-at-the-field-points). |

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
- `drift_methods` includes `control` and any generalized mode file is present.
- `cs_velocity = source` with `ISOLV = 2` (GMRES), or with `IRSP = 1` (irregular frequency removal).
- `cs_margin`, `cs_bottom_margin`, `cs_panel_size` or `cs_fs_rows` is negative, or `cs_fd_step` is outside (0, 0.5].
- `qtf` is above 0 with `mean_drift = 0`, or without `near` in `drift_methods`; `qtf` is outside 0 to 3, or `qtf_max_dw` is negative.
- `qtf_terms = body` with `qtf = 2`: the body forcing is computed for the difference QTF only, since the sum frequency needs the free-surface integral. `qtf_terms = full` with `qtf = 2` stops as well: `qtf_terms = full needs the difference QTF (qtf = 1 or 3); the sum frequency is Phase 6b.`
- `qtf_terms = full` with `qtf = 1` or `3` and a deep-water run (`H < 0`): `qtf_terms = full needs a finite water depth (the outer free-surface integral uses the finite-depth modes); set a large finite depth instead.`
- `qtf_fs_radius`, `qtf_fs_rows`, `qtf_fs_modes` or `qtf_fs_evanescent` is negative.
- `qtf_b2`, `qtf_edge_weights` or `qtf_terms` has another value than those listed, or `qtf_edge_angle` is not between 0 and 180.
- With `qtf_terms = body`, `qtf = 1` or `3` and `qtf_b2 = control`: a body has no free surface F inside its C (`qtf_b2 = control: body <b> has no free surface F inside its control surface ...`), or its C has no bottom (`qtf_b2 = control: the control surface of body <b> has no bottom ...`). Unlike the control-surface drift, which then leaves out F_z, M_x and M_y of that body, B2 has no fallback. C and F are built and checked as for `drift_methods = control`, with the same stops. For a generated C, set `cs_bottom_margin` below h minus the draft to keep a bottom; or set `qtf_b2 = body`.
- `qtf_field` is not 0 or 1, or is 1 with `qtf = 0` or with no field points in `ControlFile.in`; or a field point is above the free surface, on a hull or inside a body (see [the point rules](#partial-second-order-field-at-the-field-points)).

A key that does not apply is not an error: `qtf_b2` without `qtf_terms = body` or `full` and the difference QTF, the `qtf_fs_*` keys without `qtf_terms = full` and the difference QTF, and `qtf_edge_weights` or `qtf_edge_angle` without `qtf_b2 = control`, are ignored, and the log says so (`QTF note: ... ignored`).

### Methods and outputs

From the first-order solution, HAMS-MREL computes the mean (time-averaged) second-order loads on bodies in regular waves, by three methods:

- **Far-field** (`far`, momentum conservation): F_x, F_y and M_z of the whole system, from the far-field amplitude of the waves the bodies send out. Output `MeanDriftMomentum.8`.
- **Near-field** (`near`, pressure integration over the mean wetted surface): all six loads of each body, as the sum of six terms (waterline, flare, velocity, displacement, rotation, hydrostatic). Output `MeanDriftPressure.9`.
- **Control surface** (`control`, middle-field, Chen 2007): all six loads of each body from the field on a closed surface C around it. F_z, M_x and M_y also need the free surface between the hull and C. Output `MeanDriftControlSurface.7`.

With `mean_drift = 2`, each gives the loads of every pair of headings at the same frequency, for bidirectional seas. The terms of `near` and `control` go to `MeanDriftTerms.txt` (with `write_terms` 1 or 2). The files are described in [Output Files](output_files.md#wamit-format). When `far` runs with another method, the log prints their relative difference for F_x, F_y and M_z (M_z moved to the origin). Use it as a mesh check: a large difference means the other result has not converged.

The defaults `velocity_method = differencing` and `waterline_method = extrapolate` were kept after the [convergence study](second_order_validation.md#mean-drift-convergence-study). `differencing` is the closer to far-field on hulls without sharp edges: within 0.24% from 0.70 to 1.80 rad/s on a free rounded truncated cylinder at 7680 panels, where `source` is off by 0.12% to 2.3%. There is no switch of method by frequency.

### Which method to use

- **F_x, F_y, M_z of the whole system:** far-field is the reference on any hull.
- **Rounded hulls** (no sharp edges at the bilge or corners): near-field with the defaults, for the loads per body and for F_z, M_x, M_y. Measured on a rounded truncated cylinder, a rounded box and a bottom-mounted cylinder.
- **Hulls with sharp edges**, such as a flat bottom meeting a vertical wall:
  - F_z, M_x, M_y: use `control`. On a sharp truncated cylinder its F_z is within 2.3% of a [semi-analytical reference](second_order_validation.md#semi-analytical-reference-sharp-truncated-cylinder) at 4032 panels and 1.4% at 7168, fixed and free. Near-field F_z is 12% to 17% off there on the fixed hull and 23% to 506% on the free hull, with the wrong sign at 0.70 rad/s.
  - Horizontal loads per body: near-field converges slowly. `source` is the closer up to about ka = 1.5 on the cylinder, and on the box at the three frequencies run; `differencing` is the closer from about ka = 2.0 on the cylinder, and for M_y. Check either against far-field. ka is inferred from the deep-water dispersion relation, a = 10 m.
  - Free sharp hull at low frequency: use far-field for F_x. At 0.70 rad/s near-field F_x has the wrong sign with `differencing` at 7168 panels, and `source` is 10% off.
  - Sharp box: near-field F_y and M_z are not reliable with either velocity method on the meshes tested (F_y +12.5% with `source` at 7040 panels).
- **`source`** costs about twice the run time (+98% at 7680 panels, ifx, 32 threads) and cannot remove irregular frequencies: `IRSP = 1` is rejected, and at the two irregular frequencies tested it is 32% to 198% off. Keep away from them, or use `differencing` with a lid there.
- **Waterline:** keep `extrapolate`. `field` is first order where `extrapolate` is second order; see [Limits](#limits).
- **Settings:** `GRN_KIND = 8` changed no printed digit of F_x. `ISOLV = 2` (GMRES) matches LU to within 0.003 in F_x / (ρg) on the cases tried.

The evidence is in [Second-order validation](second_order_validation.md).

### Mesh guidance

- Refine towards the waterline. Thin top rows converge now that the free-surface image term is integrated over the panel ([Developer Documentation](developer_documentation.md#free-surface-image-term)): on a rounded cylinder on springs at 1.10 rad/s, near-field is +1.36%, +0.40%, +0.21% from far-field at 480, 1920, 4320 panels.
- Round corners and the bilge where the real hull allows. On a box, that took near-field against far-field from 23.4% (sharp, 7040 panels) to 0.41% (rounded, 8000 panels) for the force vector at 1.00 rad/s. It makes a different hull: judge it against its own far-field.
- If near-field is needed on a sharp edge, grade the rows strongly towards the edge. With `differencing`, rows of size (i/n)^4 at the edge (and quadratic at the waterline) took fixed near-field F_z on the sharp cylinder from 15% to 3.1% at 4032 panels, converging at about order 1.5 against 0.7 with cosine rows. This does not help `source`, and the free hull near heave resonance stays far off: see [graded meshes](second_order_validation.md#near-field-on-meshes-graded-towards-a-sharp-edge).
- Keep every hull row deeper than `wl_z_tolerance` times the body size below z = 0, or the run takes it for a lid panel and stops.
- Check convergence with two meshes, panels per side doubled. The rate of fall per doubling is a better check than one difference.
- Irregular frequencies: with `differencing`, use `IRSP = 1` with a lid near them. On one cylinder at 0.95 and 1.10 rad/s, a lid with thin rim panels moved near-field F_x by at most 0.42% against `IRSP = 0`.
- Control surface: the defaults of `cs_*` were used in all the studies. For small hulls in a multi-body case, a finer C (`cs_panel_size` half the default) took the sum over two bodies from 2.28% to 0.58% of far-field; see [the fixtures](../test/data/drift_cyl/README.md#control-surface-method).

### Limits

- **Generalized modes:** `near` and `control` are rejected when a `gen_mod` file is present. Far-field runs, but no integration test covers it.
- **Narrow gaps:** not tested. In the test fixtures, bodies are at least 20 m apart centre to centre. With `waterline_method = field`, no other body may be closer to the field point than `waterline_offset` times the segment length.
- **Multi-body:** far-field gives only the total of all bodies. Near-field and control per body are checked against the far-field total in the two-body fixtures.
- **`waterline_method = field` in finite depth:** term I has been seen to jump with the mesh or the offset: on 32 x 16 rows at depth 10 and offset 0.01, -11.40% at 1.0 rad/s against -1.93% and -2.06% at offsets 0.005 and 0.05 (see [the fixtures](../test/data/drift_cyl/README.md)). That was measured before the finite-depth Green function used Linton's method for R/h < 0.02 (v2.3.3), which replaces the series behind such jumps for ν·h ≤ 50. It has not been measured again since. The run logs a warning when this method is used in finite depth.
- **Control surface, F_z, M_x, M_y:** validated with `cs_shape = cylinder` only, and on finer meshes at 0.70 to 1.10 rad/s only. Not computed when C reaches the sea bed or the free surface between the hull and C cannot be built (the log says why).
- **QTFs:** the quadratic part comes from the near-field method. With `qtf_terms = body` the difference QTF also has the body forcing of the second-order potential (as Nemoh's HASBO), but not the free-surface integral (Nemoh's HASFS, the free-surface part of WAMIT's indirect method; WAMIT v7.5 user manual, ch. 15), so `QTF_Indirect.11d` is still incomplete. The sum-frequency QTF has the quadratic part only. See [QTFs](#qtfs-quadratic-part) and [the potential part](#qtfs-potential-part).
- **QTF potential part:** validated against an exact reference on one hull, the sharp truncated cylinder, heading 0, in finite depth, single body, free with and without external damping ([validation](second_order_validation.md#qtf-potential-part)). `qtf_b2 = control` was run there at 448 and 1792 panels with the default C, and at 4032 with C fixed. On a second hull, the columns of the OC4 semisubmersible, it was compared with Nemoh at one mesh only, with no exact reference. Multi-body and half-mesh runs are covered by integration tests only.
- **Second-order field at the field points (`qtf_field`):** partial. The second-order diffracted wave, and for a free body the wave radiated by its second-order motion, are omitted; near the bodies they are not small. The complete field is a later phase. See [Partial second-order field](#partial-second-order-field-at-the-field-points).

### QTFs, quadratic part

With `qtf` set, the near-field terms are evaluated for every pair of frequencies (ω_i, ω_j) of `ControlFile.in`, not only for ω_i = ω_j. Every term of the mean drift is a product of two first-order quantities, so the same formulas give the difference-frequency part (with the conjugate of the second factor) and the sum-frequency part (without it); the diagonal of `QTF_Quadratic.10d` is the mean drift of `MeanDriftPressure.9`. The definition, the conjugation and the scaling are in [Output Files](output_files.md#wamit-format). The zero and infinite frequencies of `0_inf_frequency_limits = 1` are left out.

- Use a uniform frequency list fine enough for the slow-drift periods of interest: the QTF is only known on the grid.
- The same mesh guidance as for near-field mean drift applies, with `velocity_method = source` on hulls with sharp edges, as for the mean drift. On a sharp cylinder at 448 panels both QTFs are within about one point of Nemoh's quadratic part, and their median error against an exact reference is 0.1% to 9% ([validation](second_order_validation.md#quadratic-qtfs)).
- `pre_postprocessing_python/SecondOrder/qtf_tools.py` reads the files, gives Newman's approximation from the diagonal (Newman 1974) and the spectrum of the slowly varying force.

### QTFs, potential part

With `qtf_terms = body` (and `qtf = 1` or `3`), the difference QTF also gets the body forcing of the second-order potential, by the indirect method (Molin 1979): the second-order body condition on the mean surface (Ogilvie 1983) and the second-order incident wave (Sharma and Dean 1981) are integrated against the radiation potentials ψ_k at the difference frequency Ω = ω_i − ω_j. `QTF_Indirect.11d` holds the quadratic part plus these terms, with the rows of `QTF_Quadratic.10d` ([Output Files](output_files.md#wamit-format)). The free-surface integral of the second-order potential is left out (Phase 6), so the result is the analogue of Nemoh's DUOK + HASBO, not of a complete QTF. The formulation is in [Developer Documentation](developer_documentation.md#qtf-potential-part).

Example:

```
mean_drift     = 1
drift_methods  = far near
velocity_method = source
qtf            = 1
qtf_terms      = body
qtf_max_dw     = 0.3
```

**Cost.** ψ is needed at every distinct difference Ω = ω_i − ω_j > 0 of the frequency list within `qtf_max_dw`. An Ω on the list takes ψ from the radiation solve of that frequency. Every other Ω costs one radiation-only solve after the frequency loop (no diffraction, no output file). The log gives the count: `QTF body forcing: psi at N difference frequencies, M from the main grid and K radiation-only solves`. On a uniform list of step Δω the distinct Ω are multiples of Δω, so `qtf_max_dw` sets their number (30 for 0.01 rad/s steps and `qtf_max_dw = 0.3`); a non-uniform list can need up to N(N−1)/2 of them [derived: the distinct differences of N values].

**Memory.** The store also holds ψ at every Ω (6 × the number of bodies complex values per panel), the potential part of every pair, and with `qtf_b2 = control` the fields at the points of C and F. It is allocated at set-up, and the log gives its size: `QTF: with psi and the potential part the store takes <size> MB`, with the share of C and F for `qtf_b2 = control`. Measured on the sharp cylinder with one heading and 101 frequencies, with `qtf_b2 = body` (no C and F share): 8.1 MB at 448 panels and 29.4 MB at 1792. With the default `qtf_b2 = control` the store is larger by the fields on C and F.

**B2 and the control surface.** With the default `qtf_b2 = control`, each run also takes one field pass over the points of C and F per kept frequency (the total velocity, all headings) and one per difference frequency (∇ψ_k of every mode). The pass costs about N_hull × N_C Green function evaluations [derived: one influence row of the hull per point]. With the default `cs_*` keys, the margins and the panel size of C come from the median hull panel size, so C moves in and is refined with the hull, and this cost grows faster than the first-order solve. Measured on the free sharp cylinder (Snellius, 32 threads, 101 frequencies, 30 difference frequencies): the B2 pass took 6.5 s per frequency at 448 panels (1416 points on C and F), 69 s at 1792 (3212 points) and 288 s at 4032 (5648 points). The whole run took 15 min 34 s at 448 panels and 2 h 49 min at 1792; at 4032 it would have taken 11 to 12 h, and with C held fixed it took 3 h 58 min.

- **Keep the defaults** for accuracy. B2 converges with the hull mesh and C together ([validation](second_order_validation.md#qtf-potential-part)).
- **To limit the cost** on a fine hull, hold C fixed: set `cs_margin`, `cs_bottom_margin` and `cs_panel_size` explicitly, for instance to the values the default gives on a coarser mesh. On the free sharp cylinder at 1792 panels, C and F of the 448-panel default (7.84, 7.84 and 1.96 m) took the run from 2 h 49 min to 1 h 18 min, with B2 within 6.9 / 13.8 / 6.9% (F_x / F_z / M_y) of the exact reference against 7.3 / 14.2 / 7.1% with the default. With the same C at 4032 panels, B2 kept converging (4.1 / 9.9 / 5.3%, and 0.82 / 1.63 / 1.09% away from the resonance). This is an option for cost, not needed for accuracy.
- **Fixed bodies** (`motion_state = fixed`): B1, B2 and B3 are zero, since each is a product with the motion [derived: the formulation in [Developer Documentation](developer_documentation.md#qtf-potential-part)]. `qtf_b2 = body` then gives the same result without the passes over C and F [derived].
- `qtf_edge_weights` stays `off`: see [what was tried](second_order_validation.md#tried-and-not-kept).

**Lightly damped resonances.** A free body with radiation damping only can have a sharp resonance. On the free sharp cylinder the coupled surge-pitch mode at 0.77 rad/s has a peak of about 11 m/m in surge. There the computed surge and pitch are 3% to 5.5% off at 1792 panels, and the QTF, which is quadratic in the motions, inherits about twice that at the pairs with one frequency in the peak ([validation](second_order_validation.md#the-resonance-and-the-motions)). The second-order terms themselves are within 1% to 2% away from it. Where viscous damping matters, as at a sharp bilge, add it as external linear damping in `Hydrostatic.in` (lines 11 to 16): it lowers the peak and with it this error. Measured on the free sharp cylinder with external damping in heave and pitch of 1.4 and 2.1 times the radiation damping at the 0.77 rad/s peak: the potential part was within about 3% of the exact reference over all pairs at 1792 panels, against up to 13% without the damping ([validation](second_order_validation.md#damped-free-body)). Model tests of a floating wind turbine also found its low-frequency resonant response sensitive to viscous damping (Zhang et al. 2026).

**Threads.** The field pass over C and F runs its loop over the points on `NTHREAD` threads. The total velocity comes from `fields_green`, or from `fields_source` with `cs_velocity = source` (both in [FieldVelocity.f90](../src/FieldVelocity.f90)); ∇ψ_k always comes from `fields_green`. In `fields_source` the loop over the panels is serial. In `fields_green` the loop over the panels inside it (`field_influence_row`, [WaterlineField.f90](../src/WaterlineField.f90)) is also an OpenMP loop, which runs on one thread there as long as nested parallelism is off, the OpenMP default. If nested parallelism is switched on (for instance `OMP_MAX_ACTIVE_LEVELS` above 1), up to `NTHREAD` × `NTHREAD` threads are started, and nothing limits that. Leave nested parallelism off.

Sources: Chen, X.-B. (2007). Middle-field formulation for the computation of wave-drift loads. J. Eng. Math. 59, 61-82. Newman, J. N. (1974). Second-order, slowly-varying forces on vessels in irregular waves. Proc. Int. Symp. Dynamics of Marine Vehicles and Structures in Waves, London, 182-186. Molin, B. (1979). Second-order diffraction loads upon three-dimensional bodies. Applied Ocean Research 1(4), 197-202. Ogilvie, T. F. (1983). Second-order hydrodynamic effects on ocean platforms. International Workshop on Ship and Platform Motions, Berkeley. Sharma, J. N., Dean, R. G. (1981). Second-order directional seas and associated wave forces. SPE Journal 21(1), 129-140. Cong, P., Teng, B., Chen, L., Gou, Y. (2020). A novel solution to the second-order wave radiation force on an oscillating truncated cylinder based on the application of control surfaces. Ocean Engineering 204, 107278. WAMIT v7.5 user manual, ch. 15. Zhang, Y., Tosdevin, T., Tang, T., Stansby, P. (2026). Second-order difference-frequency effects on floating wind turbine responses: full-QTF modelling and experimental comparison. arXiv:2610.27940.

### Partial second-order field at the field points

With `qtf_field = 1`, HAMS-MREL writes `QTF_FieldPartial.txt`: the second-order pressure at the field points below z = 0 and the second-order free-surface elevation at the points on z = 0, for every frequency and heading pair of `QTF_Quadratic.10d` (difference frequency) and `.10s` (sum frequency). The field points are those of `ControlFile.in` (`Number_of_field_points`). The kernels follow `qtf`, the heading pairs `mean_drift`, and `qtf_max_dw` limits the difference pairs, as for the QTFs. The columns, the conjugation and the scaling are in [Output Files](output_files.md#wamit-format).

**The file is partial.** It holds two parts:
- the quadratic part: products of the total first-order field at the point (incident, diffracted, and radiated by the motions of a free body);
- the second-order incident wave of each pair (Sharma and Dean 1981); for one wave with itself at the sum frequency this is the second-order Stokes wave (Dean and Dalrymple 1991).

It omits the second-order diffracted wave and, for a free body, the wave radiated by its second-order motion. Near the bodies the omitted part is not small:
- On the hull, the second-order diffracted potential cancels the normal velocity of the second-order incident potential (Molin 1979), so next to the hull it is of the size of the incident part.
- On the free surface it is forced by products of the first-order field that contain the scattered and radiated first-order waves (Molin 1979). Near a body that is not small against the wave length these waves are of the size of the incident wave, so this forcing is of the size of the one that makes the second-order incident wave. Kim and Yue (1989, 1990) solve this complete problem for an axisymmetric body.

So the omitted part comes from the bodies and falls off with distance from them, as their first-order waves do. It is small only where those waves are small against the incident wave, for example around a hull small against the wave length: on a hull of radius 0.05 m the file gives the Stokes wave to a relative error of at most 6.2e-5 ([validation](second_order_validation.md#partial-second-order-field)). Use the file for checks, and for the field away from the bodies with this in mind. The complete field, with the second-order diffracted and radiated waves, is a later phase.

The equations are Bernoulli's, expanded about z = 0 (Dean and Dalrymple 1991):
- pressure at a fixed point below z = 0: p⁽²⁾ = −ρ[∂φ⁽²⁾/∂t + ½|∇φ⁽¹⁾|²], the dynamic second-order pressure only. A point between the wave troughs and z = 0 is not always wet, so keep pressure points below the troughs;
- elevation at a point on z = 0: η⁽²⁾ = −(1/g)[∂φ⁽²⁾/∂t + ½|∇φ⁽¹⁾|² + η⁽¹⁾ ∂²φ⁽¹⁾/∂z∂t],

with φ⁽²⁾ the second-order incident potential. The formulas per frequency pair are in [Developer Documentation](developer_documentation.md#partial-second-order-field-qtffieldf90).

**Points.** As for `PressureElevation.6p`, a point with |z| ≤ 1e-6 is an elevation point (its z is set to 0) and a point with z < −1e-6 a pressure point. The run stops with a message naming the point and the body when a point is:
- above the free surface (z > 1e-6);
- on a hull: within one panel size of the centroid of a hull panel, with the size of that panel;
- inside a body: the hull panels (symmetry images included) and their mirror image in z = 0 enclose the point, judged by the winding number (the solid angle they subtend at the point, over 4π, is more than 0.5). This also rejects a point on z = 0 inside the waterplane of a surface-piercing body. The test needs a watertight hull.

**Route and cost.** The first-order potential and its gradient at the points always come from the Green function, whatever `velocity_method` or `cs_velocity` is, with central differences of step δ = 0.02 times the size of the nearest hull panel (the default of `cs_fd_step`). At each frequency of the QTF grid, each stencil point takes one influence row (the Green function and its normal derivative from every hull panel), which serves all headings and is not stored. A pressure point takes 7 rows: the point and x ± δ, y ± δ, z ± δ. An elevation point takes 5: the point and x ± δ, y ± δ, with ∂φ/∂z from the linear free-surface condition. A pressure point less than δ below z = 0 also takes 5 rows. On a 480-panel hull, 6 points (3 of each type) and 2 headings took 1.5 to 1.8 s per frequency on a loaded laptop. The log gives this time at each frequency.

- The store keeps φ and ∇φ at each point, 64 bytes per point, frequency and heading. It is allocated at set-up, so a run that does not fit stops before the first-order solves, and the log gives its size. `qtf_field = 1` also keeps the first-order states of every heading of a frequency, as `mean_drift = 2` does.
- In finite depth above ν·h = 50 (ν = ω²/g), the log warns at each such frequency that the field values may scatter, for the reason given under `cs_velocity`.
- `pre_postprocessing_python/SecondOrder/qtf_field_plots.py` maps the file and writes grids of field points: see [Python Scripts](python_scripts.md#partial-second-order-field-maps).

Sources: Dean, R. G. and Dalrymple, R. A. (1991). Water Wave Mechanics for Engineers and Scientists. World Scientific. Kim, M.-H. and Yue, D. K. P. (1989). The complete second-order diffraction solution for an axisymmetric body. Part 1. Monochromatic incident waves. J. Fluid Mech. 200, 235-264. Kim, M.-H. and Yue, D. K. P. (1990). Part 2. Bichromatic incident waves and body motions. J. Fluid Mech. 211, 557-593. Molin, B. (1979). Second-order diffraction loads upon three-dimensional bodies. Applied Ocean Research 1, 197-202. Sharma, J. N. and Dean, R. G. (1981). Second-order directional seas and associated wave forces. Soc. Petroleum Eng. J. 21, 129-140.
