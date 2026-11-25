# Input Files

HAMS-MREL requires a set of input files to run. The number files depend on whether you are analyzing a single body or multiple bodies:
- `ControlFile.in`
- `Hydrostatic[_i].in`
- `HullMesh[_i].pnl` (for a single body, optionally with dipole panels for zero-thickness plates)
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
| next | `26x,i16` | `ISOL` | What the diffraction outputs contain: 1 = total (incident + scattered), 2 = scattered only. Both values solve for the same scattered potential; with 1 the incident potential is added afterwards. This affects the exciting force (`.3`), and through it the motions (`.4`), and the pressure and elevation (`.6p`); added mass and damping (`.1`) do not depend on it. **Not the same as in the original HAMS**, see the note below this table. |
| next | `23x,i16` | `IRSP` | Irregular frequencies: 1 for to remove them, 0 to keep them |
| next | `23x,i16` | `NTHREAD` | Number of OpenMP threads; setting it to 1 removes parallelization |
| (optional) | keyword, then free-format integer | `ISOLV` | Linear solver: 1 = direct LU (default), 2 = GMRES with H-matrix acceleration. Line is detected by the `Solver_type` keyword; the first integer after it is read and a trailing `# comment` is allowed. If the line is absent, ISOLV defaults to 1; if the value cannot be read, a warning is printed and ISOLV defaults to 1. See [iterative_solvers.md](iterative_solvers.md) for full details. |

**`ISOL`: difference from the original HAMS.** In the original HAMS, both values give the total exciting force, by two formulations: `ISOL = 1` solves for the scattered potential and adds the incident potential in the force, and `ISOL = 2` solves for the total potential directly, with the incident potential on the right-hand side (`AssbMatx.f90`, `ASSB_DBC`, and `PotentWavForce.f90`, `EFORCE`, of github.com/YingyiLiu/HAMS). HAMS-MREL always solves for the scattered potential (`AssbMatx.f90`, `ASSB_DBC`, here), and the force routine is unchanged, so here `ISOL = 2` gives the scattered (diffraction) part only, without the Froude-Krylov force. A control file from the original HAMS with `ISOL = 2` should therefore use `ISOL = 1` in HAMS-MREL to get the same total force; `ISOL = 1` means the same in both programs.

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
| 4 | free | `NELEM, NTND, ISX, ISY[, NDIPOLE]` | Number of hull panels, number of nodes, X-symmetry flag, Y-symmetry flag and, optionally, the number of dipole panels `NDIPOLE` (single body only, see [Dipole panels](#dipole-panels-optional-single-body)). Without the fifth field `NDIPOLE` is 0 and the file is read as before |
| 5-6 | free | — | Comment or blank lines (ignored) |
| next `NTND` | free | `M` and `XYZ_LOCAL_MULTI(bodyID,nodeID,1:3)` | Node coordinates. Each line contains: `node_ID  x_local  y_local  z_local` |
| next 3 lines | free | — | Comment or blank (ignored) |
| next `NELEM` lines | free | `M`, `NCN_MULTI(bodyID,nodeID)` and `NCON_MULTI(bodyID,nodeID,1:numvertices)`| Panel definitions. Each line contains: panel number, number of vertices and (vertex1, vertex2, ...) |
| next line (only if `NDIPOLE` > 0) | free | — | Comment line, `# Start Definition of Dipole Panels` (ignored) |
| next `NDIPOLE` lines (only if `NDIPOLE` > 0) | free | `M`, `DNCN(panelID)` and `DNCON(panelID,1:numvertices)` | Dipole panel definitions, in the same format as the hull panels. The vertex numbers refer to the same node list as the hull panels |

### Dipole panels (optional, single body)

Thin plates, such as heave plates and fins, can be given as sheets of zero thickness instead of closed bodies. Each dipole panel carries the jump of the potential across the sheet, μ = φ(+n side) − φ(−n side), instead of a potential. This is a normal-dipole sheet (Hess and Smith 1964, Sec. 4.2; the dipole panels of the WAMIT User Manual). A sheet has no interior, and there are no faces close together. Close faces matter: in the heave-plate lid study of 2026-11-24 (IRSP = 1, panels 0.75 m wide), the lid moved the heave added mass above 3.3 rad/s by 27 % for a closed plate 0.1 m thick, by 2 % for one 0.5 m thick, and by 0.0 % for the column alone.

The dipole panels go in `HullMesh.pnl`, after the hull panels: give their number as the fifth field `NDIPOLE` of line 4, then, after the last hull panel, one comment line and the `NDIPOLE` dipole panel lines (see the table above). Hull and dipole panels share one node list, so a node may be used by both. The first number of each panel line is not used; the dipole panels are numbered from 1 by convention. Lines 1 to 3 are skipped whatever they hold, so the comment on line 3 may mention the dipole panels. `NELEM` must be at least 1: a body made of dipole panels only (a free disc or plate, `NELEM = 0`) is not supported, and the run stops at input. Example of the panel part of a file with 192 hull panels and 96 dipole panels:

```
  # Start Definition of Node Relations   ! panel_number  number_of_vertices   Vertex1_ID   Vertex2_ID   Vertex3_ID   (Vertex4_ID)
     1     4       1       2       3       4
   ...
   192     4     741     742     743     744
    # Start Definition of Dipole Panels
     1     4     745     746     747     748
   ...
    96     4    1125    1126    1127    1128
    # End Definition of Node Relations
```

The block is read by `ReadBodyMesh` and set up by `DIPOLE_SETUP`, both in [ReadPanelMesh.f90](../src/ReadPanelMesh.f90). The panel normal follows the node order, as for hull panels. Either direction is allowed: flipping it changes the sign of μ and nothing else. Without the fifth field, or with `NDIPOLE = 0`, nothing changes.

The run stops at input if:

- any node of a dipole panel is less than one panel size below z = 0. Sheets must be fully submerged, and the wave part of the dipole equations is tested only for sheets deeper than their panel size;
- a dipole panel lies in the plane of a hull panel and overlaps it (the centroid of either inside the other). A sheet cannot close a body or lie on its surface. For a heave plate, the column's bottom disc is part of the hull, and the sheet is only the annulus outside it, in the same plane;
- a dipole centroid lies within 1e-6 panel sizes of a hull centroid (a sheet must not coincide with hull panels);
- the symmetry flags of the hull are not 0, `ISOLV` is not 1, there is more than one body (a fifth field above 0 in a `HullMesh_i.pnl` stops the run), or `SecondOrder.in` is active. Version 1 supports one body, no symmetry, the direct solver and first order only.

A sheet may meet the hull along an edge, for example a heave-plate annulus at the bottom of a column. Nothing is added at the junction. Example, the heave plate of Kurnia and Ducrozet (2023, Sec. 4.2.3), a column of radius 10 m and draft 9.9 m with a plate of radius 15 m: the hull panels of `HullMesh.pnl` are the column wall down to z = −9.9 m and the column's bottom disc (r < 10 m) at z = −9.9 m; its dipole panels are only the annulus 10 m ≤ r ≤ 15 m at z = −9.9 m. The bottom disc and the annulus lie in one plane and do not overlap. Grade the annulus toward its free outer rim. On this plate, with 0.75 m panels around, the added mass in heave fell by 3.5 % of the peak of the WAMIT curve from 7 to 28 radial rows (7, 14, 28 rows with radii r = 10 + 5 sin(π i / 2N) m), and by about 1 % from 14 to 28; outer rows 1:95 (7.9 mm by 0.75 m) were fine. The sheet replaces a closed plate 0.1 m thick; without the lid the two agree within 1.2 % in A11, A33, A55, B33 and B55 on the same panel size (see the [developer documentation](developer_documentation.md#dipole-zero-thickness-panels)). At a free edge μ goes to zero like the square root of the distance to the edge (Lamb 1932, Sec. 108, for a disc). Constant panels, collocated at their centroids, converge only as O(h) there. For a disc of radius a, A33 is 13 % above the exact (8/3) ρ a³ with 96 panels, 8.2 % above with 384, 4.4 % with 1536 and 2.7 % with 3072 (6, 12, 24 and 48 rings). Grading the rings toward the rim (r_i = a sin(π i / 2N)) changes this little. So use fine sheet meshes, check a refined mesh, and, where one number is needed, extrapolate: 2 A(fine) − A(coarse), with the fine mesh refined twice in each direction, was within 0.7 % for the disc. See the [developer documentation](developer_documentation.md#dipole-zero-thickness-panels).

The fixture `test/data/dipole/column_plate/nolid/input/HullMesh.pnl` is a column with a heave-plate annulus in this layout.

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

### Lid mesh requirements

With `IRSP = 1`, the waterplane mesh (the "lid") only supplies extra collocation points. HAMS-MREL places one point at each lid panel centroid and adds one equation there: Green's identity for a point inside the body, which has no free term. The hull potentials remain the only unknowns, and the hull and lid equations are solved together by least squares (`AssbMatx_irr.f90`, `AssbMatx_irrMulti.f90`). The lid panels are never sources. Their normals, areas and sizes are not used, so only the centroids matter. The rules are:

- **All nodes on z = 0 exactly.** The reader stops on any lid node with |z| > 1e-10 (`ReadWTPLMesh`, `ReadWTPLMeshMulti`).
- **Lid panels go only in `WaterPlaneMesh[_i].pnl`, never in `HullMesh[_i].pnl`.** A hull panel on z = 0 is a source panel in the free surface, and `CalNormals` stops on it.
- **Cover the whole interior waterplane** of each body that cuts the free surface. The lid edge should follow the hull waterline; using the hull waterline nodes as the outer lid nodes is the simplest way to do that.
- **Keep the outer centroids away from the hull wall.** The distance from each lid centroid to the wall should be at least about 0.3 to 0.5 times the height of the hull's waterline panel row. Near-field terms of the Green function are partly evaluated at one point per source panel, and that is poor when a lid point lies much closer to a hull panel than the panel is large. A ring of lid panels built on the hull waterline nodes, with a ring width about equal to the waterline row height, puts its centroids at about half that height from the wall.
- **Normals pointing up (+z)** are the convention used in the examples. The solver does not use the lid normals.

Why the distance matters: on a column with a heave plate (waterline rows 0.68 m high, 832 hull panels), a lid made by clipping a square grid to the waterline circle left thin triangles at the wall, with centroids 2 to 5 cm from it. With that lid, heave damping and heave excitation between 2 and 4 s were off from an independent reference solution by 68% and 30% on average. Without a lid they were off by 12% and 9%, because of the irregular frequency in that band. A lid of three rings on the 32 waterline nodes, with outer centroids 0.38 m from the wall, brought them to 5% and 0.7%. At the zero- and infinite-frequency limits both lids gave the no-lid added mass to within 0.04%, and between 5 and 30 s to within 0.3% on added mass and 1.1% on damping.


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
| `moment_origin` | `fixed`, `moving` | `fixed` | The point the second-order moments are taken about, for the mean drift of all three methods and both QTFs. `fixed`: `XR` of each body at its mean, space-fixed position, as WAMIT writes i = 4, 5, 6 (WAMIT v7.6 manual, sec. 5.2 and eqs. 15.53 to 15.57). `moving`: `XR` moving with the body, as WAMIT writes i = -4, -5, -6 and as Nemoh and HydroStar give them. The two differ by $P(\xi \times F^{(1)})$, with $F^{(1)}$ the first-order fluid force. For a free body with its centre of gravity at `XR` this is zero for the mean drift but not for the QTFs: up to 80 % of the largest pitch QTF on OC3 Hywind near the heave resonance. The file formats are the same, and the log says which origin was used. The `QTF_Motion` files always use `moving`, whatever this key says. The point is always `XR` of each body, the point the first-order moments are about too; there is no separate key for it. To compare with another solver, set `XR` to that solver's reference point. See [Moment reference point](output_files.md#moment-reference-point-note). |
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
| `qtf_terms` | `quadratic`, `body`, `full` | `quadratic` | Which parts of the QTF. `quadratic`: the quadratic part only (products of first-order quantities). `body`: also the body forcing of the second-order potential, by the indirect (Haskind) method with the radiation potentials ψ at the difference frequencies, written to `QTF_Indirect.11d`. The difference QTF (`qtf = 1` or `3`) by the indirect method, and the sum QTF (`qtf = 2` or `3`) by the indirect method too, written to `QTF_Indirect.11s` with or without `qtf_direct = 1` (Phase 6b; with `qtf_direct = 1` also by the direct method). With `body` the free-surface integral of the second-order potential is left out, and for the sum the log warns (see [Indirect sum frequency](#qtfs-potential-part)). `full`: `body` plus the free-surface integral (Phase 6a, finite depth only, H > 0), so that `QTF_Indirect.11d` is the complete indirect-method QTF; it implies everything `body` does, and `qtf_b2`, `qtf_edge_weights` and `qtf_edge_angle` apply to it as to `body`. With `full` the free-surface terms are in `.11d`, and `QTF_PotentialTerms.txt` gives them as two more terms, `FS_inner` and `FS_outer` ([Output files](output_files.md)). A mesh from `FSInnerMesh.pnl` must reach the circle R: the log gives the share of the free surface inside R that it covers, warns when it is off by more than 1% (short of R, or overlapping), and stops above 10%. **Accuracy:** on the sharp truncated cylinder, free body, the free-surface terms on the default generated inner mesh are within 12 / 4.3 / 2.4% (F_x / F_z / M_y) of the exact value at m1 (448 panels), 4.5 / 3.2 / 0.9% at m2 and 2.2 / 1.3 / 0.4% at m3, away from the surge–pitch resonance ([validation](second_order_validation.md#qtf-free-surface-integral)). They get more accurate with a finer hull mesh, because the inner mesh takes a column at every hull waterline node. On a coarse hull, check them on a finer hull mesh; a larger `qtf_fs_radius` or a smaller `qtf_fs_panel` is a second check (at m1 these gave 2–3.5%, but more radial rows alone, `qtf_fs_rows`, made it worse). The log notes this on every run with the generated mesh. `FSInnerMesh.pnl` gives a mesh of your own. Needs `qtf` 1, 2 or 3, with or without `qtf_direct = 1`; with `qtf = 0` `body` and `full` do nothing. With the sum (`qtf` 2 or 3) `full` needs finite depth, as for the difference. See [QTFs, potential part](#qtfs-potential-part). |
| `qtf_max_dw` | real, at least 0 | `0` | Largest \|ω_i − ω_j\| in rad/s written to `QTF_Quadratic.10d`, `QTF_Indirect.11d` and `QTF_Direct.12d`. 0 writes every pair. With `qtf_terms = body` it also limits the difference frequencies at which ψ is solved. It does not limit the sum-frequency files: `QTF_Quadratic.10s` has every pair, and `qtf_max_sw` limits `QTF_Direct.12s`. Choose it by a convergence study of the slow-drift response, widening the bandwidth until the result no longer changes, as Zhang et al. (2026) do for the full difference-frequency QTF of a floating wind turbine. The validation used 0.3 rad/s. |
| `qtf_max_sw` | real, at least 0 | `0` | Largest ω_i + ω_j in rad/s of the sum-frequency pairs (`qtf = 2` or `3` with `qtf_terms = body` or `full`, with or without `qtf_direct = 1`). 0 keeps every pair. The limit is inclusive (to 1e-9 relative) and applies to `QTF_Indirect.11s`, `QTF_PotentialTermsSum.txt` and `QTF_Direct.12s`, not to `QTF_Quadratic.10s` (`sum_pair_index`, [QTFDriver.f90](../src/QTFDriver.f90)); a sum row of `QTF_Field.txt` outside it has `has_direct` = 0. With `If_remove_irr_freq = 0` (IRSP = 0) the sums at or above the estimated lowest irregular frequency are left out as well, so the cutoff is the smaller of the two; the log names it and the number of pairs left out, and says that `If_remove_irr_freq = 1` keeps them. The estimate is the minimum over the closed waterline loops of √(g j₀₁ / a_eq), with j₀₁ = 2.404826 and a_eq = √(A/π), A the area of the loop: a lower bound for a wall-sided hull (Faber–Krahn: the disc of equal area has the lowest first Dirichlet eigenvalue; Pólya and Szegő 1951) and 1.536 rad/s for the sharp cylinder of radius 10 m, whose exact value is 1.548. Radiation-only solves are needed at every sum frequency kept, so a smaller `qtf_max_sw` also cuts their number. |
| `qtf_direct` | `0`, `1` | `0` | With `qtf_terms = body` or `full` and `qtf = 1`, `2` or `3`: `1` also solves the second-order boundary-value problem for the diffracted and radiated potential directly, at every difference frequency including those on the main grid, after the frequency loop (one more matrix assembly and factorisation per such frequency). With `body` the solve has the body forcing only; with `full` it has the free-surface forcing of φ_D as well, so that `.12d` is the complete direct-method QTF (finite depth only, as `full` always; slower, see [Direct method](#qtfs-potential-part)). It writes `QTF_Direct.12d` (and with `write_terms` 1 or 2 `QTF_DirectTerms.txt`), see [output files](output_files.md). With `qtf` 2 or 3 it also solves at every sum frequency ω_i + ω_j and writes `QTF_Direct.12s` (and with `write_terms` 1 or 2 `QTF_DirectTermsSum.txt`), see the sum-frequency paragraph of [Direct method](#qtfs-potential-part); use it with `qtf_terms = full`. With `If_remove_irr_freq = 1` the mesh must have lid points, else the run stops. **`qtf = 3` with `qtf_terms = full`:** the sum is on by default, with `qtf_direct = 0` as with 1, and the inner free-surface mesh and the mode count M of the outer region then follow the largest sum frequency, not the largest frequency of the list (`fs_setup(..., omega_sum_max)`, [SecondOrderDriver.f90](../src/SecondOrderDriver.f90) lines 379-381 and 483). The difference files (`.11d`, `.12d`, `QTF_DirectTerms.txt`) are then equal to those of the same input with `qtf = 1` only when `qtf_fs_panel` and `qtf_fs_modes` are both given (and the mesh is generated); otherwise they differ a little, because the inner mesh and M are others. This holds with `qtf_direct = 0` too: the difference files are then `.11d` and `QTF_PotentialTerms.txt`. The run also costs more: up to 2N − 1 extra frequencies with their own solves, the sum pair loop, and an inner mesh with many more points (2600 measured on the sharp cylinder at m1, against about 186 for the difference alone). Give both keys, or use `qtf = 1` for the difference alone, to keep the cost and the difference files as before. Other values stop the run. |
| `qtf_motion` | `0`, `1` | `0` | `1` computes the second-order motion of free bodies from the second-order force, and writes the second-order motion files: `QTF_Motion.11d` and `.11s` from `QTF_Indirect.11d` and `.11s`, and with `qtf_direct = 1` also `QTF_Motion.12d` and `.12s` from `QTF_Direct.12d` and `.12s` ([output files](output_files.md#wamit-format)). With `qtf_field = 1` and `qtf_direct = 1` it also adds the wave radiated by that motion to `QTF_Field.txt`. Needs `qtf_terms = body` or `full` and `qtf = 1`, `2` or `3`; else the run stops. A fixed body (`motion_state = fixed`) gets no motion and no motion file, and the log says so. The motion takes the linear damping `BLNR` only, not `BQDR`, and an Ω whose matrix is singular (Ω = 0 for a body without mooring) is left out of the files: see [Second-order motion](#second-order-motion). Other values stop the run. |
| `qtf_b2` | `control`, `body` | `control` | With `qtf_terms = body`: where the surface piece of the second-derivative term B2 is integrated. `control` moves it from the hull onto the control surface C and the free surface F around each body by a Gauss identity (Cong et al. 2020; Chen 2007, App. C); the rest of B2 stays on the hull. `body` integrates it on the hull, which is inaccurate next to sharp edges ([validation](second_order_validation.md#qtf-potential-part)). `control` builds C and F with the `cs_*` keys, also when `drift_methods` has no `control`, and needs a C with a bottom. Ignored, with a note in the log, without `qtf_terms = body` or `full` and `qtf` 1, 2 or 3. `control` moves B2 to the control surface for the difference QTF and for the sum QTF `QTF_Indirect.11s`, with any `qtf` (`uses_b2_control`, SecondOrderDriver.f90:891-900). `potential_part_sum` (QTFDriver.f90:1994-2005) then takes the control-surface form, as `potential_part` does for the difference. So the same input gives the same `.11s` for `qtf = 2` and `qtf = 3`. Since Phase 6b, `qtf = 2` with the default `qtf_b2 = control` also builds C and F (and writes `CSMesh_<b>.pnl` when `drift_methods` has no `control`) and stops on a hull without a C with a bottom (for example a bottom-mounted hull). Set `qtf_b2 = body` for such a hull. |
| `qtf_edge_weights` | `off`, `on` | `off` | With `qtf_b2 = control`: product weights on the hull panels next to a sharp edge, for the hull term of B2 ([SharpEdge.f90](../src/SharpEdge.f90), velocity ~ ρ^p with p = π/α_f − 1, α_f the fluid wedge angle). Off by default: with HAMS's panel velocities it moved B2 further from the exact reference ([validation](second_order_validation.md#tried-and-not-kept)). Ignored, with a note, without `qtf_b2 = control`. |
| `qtf_edge_angle` | real, above 0 and below 180 | `45` | With `qtf_edge_weights = on`, and for the B2 one-point mask of the direct method (`qtf_direct = 1`): a mesh edge is sharp when the normals of its two panels differ by more than this many degrees. |
| `qtf_fs_radius` | real, at least 0 | `0` | With `qtf_terms = full`: the radius R in metres of the circle around all bodies that splits the free surface into an inner mesh (waterline to R) and an outer region (r > R, partial waves). 0 sets the value automatically: r_x + max(0.5 r_x, 2 × the largest hull panel), r_x the larger of the waterline radius and the hull radius (the largest horizontal distance of a hull node from the centre). R must exceed both, else the run stops (a hull wider below the waterline would put the partial waves inside the body). Ignored, with a note in the log, without `qtf_terms = full` and `qtf` 1, 2 or 3 (`uses_fs`, SecondOrderDriver.f90:211-213 and 887). |
| `qtf_fs_rows` | integer, at least 0 | `0` | With `qtf_terms = full`: the number of radial rows of the generated inner mesh. 0 sets the value automatically. Ignored, with a note, as `qtf_fs_radius`. |
| `qtf_fs_panel` | real, at least 0 | `0` | With `qtf_terms = full`: the panel size in metres of the generated inner mesh. 0 sets λ/10 at the largest frequency of the run (λ from the finite depth). The size sets the number of chords of the polygon on the circle R, ceil(2πR/size) (at least 8), and the least number of radial rows (the outer row is not thicker than the size; a larger `qtf_fs_rows` is kept); the columns along the waterline follow the hull panels as before. It must be below R − r_wl (the width of the inner free surface, r_wl the largest distance from the centre to a waterline node), else the run stops with `qtf_fs_panel = <size> m must be below R - r_wl = <width> m`. Not used with `FSInnerMesh.pnl`. Ignored, with a note, as `qtf_fs_radius`. With the sum on (`qtf` 2 or 3, with or without `qtf_direct = 1`) the default follows the largest sum frequency (SecondOrderDriver.f90:379-381 and 483), and giving this key together with `qtf_fs_modes` keeps the difference files equal to those of `qtf = 1` (see `qtf_direct`). |
| `qtf_fs_modes` | integer, at least 0 | `0` | With `qtf_terms = full`: the number of angular modes M of the outer region. 0 sets the value automatically. Ignored, with a note, as `qtf_fs_radius`. With the sum on (`qtf` 2 or 3, with or without `qtf_direct = 1`) the automatic value follows the largest sum frequency (see `qtf_direct`). |
| `qtf_fs_evanescent` | integer, at least 0 | `0` | With `qtf_terms = full`: the number N_ev of evanescent vertical modes of the outer region, besides the propagating mode. `0` sets N_ev = max(8, ceil(9 h / (pi (R - r_wl)))), r_wl the largest distance from the centre of the circle to a waterline node, so that the slowest mode kept decays by about e^-9 from the waterline to R [derived: k_n is about n pi / h]. The log gives the value used. Ignored, with a note, as `qtf_fs_radius`. |
| `qtf_fs_circle` | `samples`, `panels` | `panels` | **On by default** (`panels`). With `qtf_terms = full`: how the partial waves of the outer region (r > R) are found. `panels` takes their amplitudes from the hull integral, by Graf's addition theorem (Abramowitz and Stegun 9.1.79, 9.6.3 and 9.6.4): no field is sampled on the circle, and the mode check uses every 8th of the circle samples. `samples` is the old route: the field is computed at every sample on the circle and decomposed in modes. On the panels route the check samples take the finite-depth field, so above ν h = 50 the mode check can warn on the samples' own scatter (the Limes series of the finite-depth function), while the ν h > 50 warning of the free-surface passes is not printed when those passes take deep-water rows; the coefficients do not depend on the check. Measured (sharp cylinder, m2, fixed and free body, sum and difference kernels, against the exact references; `validation-runs/fs-speedups/results/g1_compare.txt`): 126 rows, none worse than the samples route by more than max(0.5 percentage points, 10 %); the largest change is in the sum kernel P_FS F_z of the fixed body, 3.38 % (samples) to 2.36 % (panels). On a 448-panel cylinder the panels route is 3.6e-9 from the exact direct sum, against 1.7e-6 for the samples (unit test `pw_from_panels_matches_samples`). On OC4 (WAMIT model, 5022 panels) `FS_outer` moves by at most 2e-4 of its maximum. **The deep-water field rows (`field_deep_kh`) apply on the samples route only:** the panels route uses the finite-depth field for the partial waves. At the default `field_deep_kh = 10` this moves the results little (OC4: `FS_outer` within 2e-4; whole-run results moved by more, up to 6.6 % below 1.55 rad/s and 9.8 % above in mode 1 of `QTF_Indirect.11s`, in a run that also took deep-water rows on C and F, and the cause is open until a re-run, see the [validation document](second_order_validation.md#oc4-with-the-new-defaults)); with a small `field_deep_kh` it can (fixture `qtf_psi_reuse`, `field_deep_kh = 0.25`: 1.2e-2 of the largest value in `QTF_Indirect.11d`), and then `qtf_fs_circle = samples` keeps the deep-water rows. The run stops if a panel (with its images) is not inside the circle. Ignored, with a note in the log, without `qtf_terms = full` and `qtf` 1, 2 or 3. |
| `qtf_fs_inner` | `stencil`, `source` | `stencil` | With `qtf_terms = full`: how the fields on the inner free-surface mesh are found. `stencil` is the route of the earlier versions: the second derivative of the potential comes from a stencil around each point, with the near-waterline fit. `source` takes them from the source strengths and moves the second derivative onto the boundary by an integration by parts. **The source route is off by default because it failed the accuracy check against exact data (sharp m1 hull, 448 panels, fixed body): the inner integral at the sum pair (1.4, 1.0) rad/s is 3.58 % off against 2.05 % for the stencil, ∇φ at 0.1 m from the waterline 5.84 % against 4.92 %, and its free-surface terms move by 16 to 25 % between the 220- and 960-panel hulls, against 0.3 to 1 % for the stencil** (`validation-runs/fs-speedups/gate/G1_table.md`). It stops, with a message, when combined with `qtf_direct = 1`, `IRSP = 1` or `ISOLV = 2`, and when the generated or given inner mesh does not close on the waterline (a miss of more than 1e-6 segment lengths). Ignored, with a note in the log, without `qtf_terms = full` and `qtf` 1, 2 or 3. |
| `field_deep_kh` | real, at least 0 | `10` | **On by default** (`10`). With a finite depth, the free-surface field points of the second-order passes (z = 0: the free-surface fields and the waterline fit, and the free-surface ψ) use the deep-water Green function instead of the finite-depth one at a frequency with k0 h at or above this value (k0 the wave number, h the depth). This makes these passes cheaper. **Free surface only:** below z = 0 the deep-water function misses the seabed images of the finite-depth one, so at control-surface bottom points it is 2.3e-3 to 9.3e-3 of the row maximum off on an OC4-like case and 0.15 to 0.19 off with 1 m of bottom clearance, against 7e-5 for the finite-depth function (unit test `field_row_deep_seabed`); ∇ψ and the fields on the control surface C and the free surface F therefore keep the finite-depth function. Also not switched: the hull solve, the waterline rows of the near-field drift, the field points of `qtf_field` and the direct method, and the partial waves on the circle when `qtf_fs_circle = panels` (the default; they take the finite-depth field, and the log says so when k0 h of the largest frequency reaches this value). The deep-water rows on the circle apply with `qtf_fs_circle = samples` only. `0` switches it off; with `qtf_fs_circle = samples` as well it gives the old results exactly (the finite-depth function at every frequency). Measured with the key on against `field_deep_kh = 0` (h = 100 m, sum frequencies up to k0 h = 45; test `qtf_deep_kh_within_limits`): `QTF_Indirect.11s` 7.7e-4 and `QTF_PotentialTermsSum.txt` 1.07e-3 of the largest value of the file, `QTF_Indirect.11d` 1.2e-6 and `QTF_PotentialTerms.txt` 7.1e-5. **Accuracy (z = 0):** the deep-water rows are within 0.13 % of the row maximum from k0 h = 10 up to ν h = 50 (ν = ω²/g). Above ν h = 50 they are more accurate than the finite-depth field mode, measured against a direct series sum; at k0 h = 40 the finite-depth mode itself is 3.3e-3 off the direct sum, the deep-water row 8e-5 (`field_row_deep_within_tol`). At k0 h = 3.3 the difference is 1.6e-2. Whole-run results move by more than the rows: the OC4 run in the validation document (made with deep-water rows on C and F too) moved by up to 6.6 % below 1.55 rad/s and 9.8 % above, and its cause is open until a re-run. Sources: `validation-runs/fs-speedups/results/gbench/gbench_summary.txt` (benchmark of the two Green functions at the points of the OC4 WAMIT hull and circle) and `validation-runs/fs-speedups/results/gbench/outlier_check.txt`. A negative value stops the run. |
| `field_deep_draft` | real, at least 0 | `10` | **On by default** (`10`). The depth guard of `field_deep_kh`: the deep-water rows are used only where the depth is at least this many times the draft of the hull (its deepest node, all bodies together). In shallower water the free-surface passes keep the finite-depth Green function at every frequency, and the log says so once ("the finite-depth Green function is kept from k0 h = ..."). Reason: in water shallow relative to the draft the deep-water function misses the seabed by more than k0 h accounts for. On the sharp cylinder (radius 10 m, draft 10 m, fixed body; semi-analytical reference at each depth) the deep-water rows at k0 h >= 10 moved the free-surface part of the sum QTF, F_z, from 0.72 % to 7.29 % of the largest exact value at h = 25 m (2.5 drafts) and from 0.74 % to 2.96 % at h = 50 m (5 drafts); both fail the accuracy gate (at most 0.5 percentage points), and with finite-depth rows the error stays below 0.75 %. At h = 100 m (10 drafts) the same rows pass it (largest change 0.14 percentage points at k0 h = 5). So the default is 10 drafts, the shallowest depth measured to pass; between 5 and 10 drafts it was not measured. Sources: `validation-runs/fs-speedups/deep_inner/findings.md`. `0` switches the guard off. A negative value stops the run. |
| `qtf_field` | `0`, `1` | `0` | 1 writes `QTF_FieldPartial.txt`: the second-order pressure and free-surface elevation at the field points of `ControlFile.in`, for the pairs of `QTF_Quadratic.10d` and `.10s`. Needs `qtf` 1, 2 or 3 and at least one field point. **Partial:** the second-order diffracted and radiated waves are omitted, and near the bodies they are not small. With `qtf_direct = 1` it also writes `QTF_Field.txt`, which adds the direct part of the second-order potential (body forcing only) as separate columns. See [Partial second-order field](#partial-second-order-field-at-the-field-points). |

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
- `qtf = 2` with `qtf_direct = 0` is no longer an error: since Phase 6b it gives `QTF_Indirect.11s` (the check was removed from `SecondOrderInput.f90`).
- `qtf_terms = full` with `qtf` 1, 2 or 3 and a deep-water run (`H < 0`): `qtf_terms = full needs a finite water depth (the outer free-surface integral uses the finite-depth modes); set a large finite depth instead.`
- `field_deep_kh`, `field_deep_draft`, `qtf_fs_radius`, `qtf_fs_rows`, `qtf_fs_panel`, `qtf_fs_modes` or `qtf_fs_evanescent` is negative, or (with a generated inner mesh) `qtf_fs_panel` is not below R − r_wl.
- `qtf_fs_circle` or `qtf_fs_inner` has another value than those listed (`Invalid value "<word>" for key "<key>" in SecondOrder.in.`).
- `qtf_fs_inner = source` with `qtf_terms = full` and `qtf` 1, 2 or 3, and `ISOLV = 2` (`qtf_fs_inner = source needs the direct solver (ISOLV=1) in ControlFile.in.`), `IRSP = 1` (`qtf_fs_inner = source cannot be used with irregular frequency removal (IRSP=1) in ControlFile.in.`) or `qtf_direct = 1` (`qtf_fs_inner = source cannot be used with qtf_direct = 1 (the direct solve needs phi_zz at the inner points); use qtf_fs_inner = stencil.`), or with an inner mesh (`FSInnerMesh.pnl`) whose edge does not lie on the waterline nodes.
- `qtf_fs_circle = panels` with a panel (or its image) outside the circle of radius R, or with k_n R too large for the evanescent modes; set `qtf_fs_radius` larger or use `qtf_fs_circle = samples`.
- `qtf_direct` is not 0 or 1, or `qtf_direct = 1` without `qtf_terms = body` or `full` (`qtf_direct = 1 needs qtf_terms = body or full in SecondOrder.in.`), or with `qtf = 0` (`qtf_direct = 1 needs qtf = 1, 2 or 3 in SecondOrder.in.`); `qtf_max_sw` is negative (`qtf_max_sw must not be negative in SecondOrder.in.`).
- `qtf_direct = 1` and the direct solve did not run at a frequency of the list, so that `.12d` or `.12s` would be incomplete (`QTF direct: qtf_direct = 1, but the direct solve did not run at Omega = ... rad/s; QTF_Direct.12d or .12s would be incomplete.`). The check is in `check_direct_kept` ([QTFDriver.f90](../src/QTFDriver.f90)); it guards against a bug, not against a choice of input.
- `qtf_b2`, `qtf_edge_weights` or `qtf_terms` has another value than those listed, or `qtf_edge_angle` is not between 0 and 180.
- With `qtf_terms = body` or `full`, `qtf` 1, 2 or 3 and `qtf_b2 = control`: a body has no free surface F inside its C (`qtf_b2 = control: body <b> has no free surface F inside its control surface ...`), or its C has no bottom (`qtf_b2 = control: the control surface of body <b> has no bottom ...`). Unlike the control-surface drift, which then leaves out F_z, M_x and M_y of that body, B2 has no fallback. C and F are built and checked as for `drift_methods = control`, with the same stops. For a generated C, set `cs_bottom_margin` below h minus the draft to keep a bottom; or set `qtf_b2 = body`.
- `qtf_motion` is not 0 or 1 (`qtf_motion must be 0 or 1 in SecondOrder.in.`), or `qtf_motion = 1` without `qtf_terms = body` or `full` (`qtf_motion = 1 needs qtf_terms = body or full in SecondOrder.in.`) or with `qtf = 0` (`qtf_motion = 1 needs qtf = 1, 2 or 3 in SecondOrder.in.`). `qtf_motion = 1` with `motion_state = fixed` is accepted: the log says ` QTF note: qtf_motion = 1 with motion_state = fixed: no second-order motion is computed.` (`qtf_init`, [QTFDriver.f90](../src/QTFDriver.f90)).
- `qtf_field` is not 0 or 1, or is 1 with `qtf = 0` or with no field points in `ControlFile.in`; or a field point is above the free surface, on a hull or inside a body (see [the point rules](#partial-second-order-field-at-the-field-points)).

A key that does not apply is not an error: `qtf_b2` without `qtf_terms = body` or `full` and `qtf` 1, 2 or 3, the `qtf_fs_*` keys (also `qtf_fs_panel`) without `qtf_terms = full` and `qtf` 1, 2 or 3, and `qtf_edge_weights` or `qtf_edge_angle` without `qtf_b2 = control`, are ignored, and the log says so (`QTF note: ... ignored`).

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
- **Hulls with slender members through the waterline**, such as the OC4 semisubmersible with its inclined cross braces: check every hull velocity against far-field. On OC4 at 5022 panels (deep water, 0.8 to 1.85 rad/s), neither velocity method met the far-field check for F_x. Near-field F_x was off by up to 7.6% with `differencing` and up to 8.8% with `source`, with opposite signs at most frequencies. `control` with `cs_velocity = source` was +2% to +10% off, since it uses the same source solve. `control` with `cs_velocity = green` met the check: within 0.8% with C and F panels of about 1 m, and 0.3% at 0.5 m. With 3.3 m panels it was off by up to 8%. On such hulls use `control` with `cs_velocity = green` and C panels of 1 m or less, and check its F_x against far-field before reading F_z.
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
- **QTFs:** the quadratic part comes from the near-field method. With `qtf_terms = body` the difference QTF also has the body forcing of the second-order potential (as Nemoh's HASBO), but not the free-surface integral (Nemoh's HASFS, the free-surface part of WAMIT's indirect method; WAMIT v7.5 user manual, ch. 15), so `QTF_Indirect.11d` is incomplete. With `qtf_terms = full` (finite depth) it has the free-surface integral too and is the complete difference QTF of the indirect method; validated on one hull, the sharp cylinder, single body, heading 0 ([validation](second_order_validation.md#qtf-free-surface-integral)). The sum-frequency QTF has the quadratic part, and with `qtf_direct = 1` and `qtf_terms = full` also the second-order potential by the direct method (`QTF_Direct.12s`, the sum-frequency paragraph of [the potential part](#qtfs-potential-part)); an indirect sum-frequency QTF is Phase 6b. See [QTFs](#qtfs-quadratic-part) and [the potential part](#qtfs-potential-part).
- **QTF potential part:** validated against an exact reference on one hull, the sharp truncated cylinder, heading 0, in finite depth, single body, free with and without external damping ([validation](second_order_validation.md#qtf-potential-part)). `qtf_b2 = control` was run there at 448, 1792 and 4032 panels with the default C, and at 1792 and 4032 with C fixed; at 4032 panels the two give the same accuracy. On a second hull, the columns of the OC4 semisubmersible, it was compared with Nemoh at 1224 and 4608 panels and with a public WAMIT run at both, with no exact reference. Multi-body and half-mesh runs are covered by integration tests only.
- **Second-order field at the field points (`qtf_field`):** partial. The second-order diffracted wave, and for a free body the wave radiated by its second-order motion, are omitted; near the bodies they are not small. With `qtf_direct = 1`, `QTF_Field.txt` has the direct part of the diffracted wave at the difference frequencies and, with `qtf` 2 or 3, at the sum frequencies, without its free-surface forcing (with `qtf_terms = full` that forcing is in `.12d`, but solved apart and left out of the field), so it is still partial. With `qtf_motion = 1` it also has the radiated wave of the second-order motion. The complete field is a later phase. See [Partial second-order field](#partial-second-order-field-at-the-field-points).
- **Second-order motion (`qtf_motion`):** validated against an exact reference on one hull, the free sharp truncated cylinder with springs in surge, sway and yaw, heading 0, depth 100 m, at two meshes ([validation](second_order_validation.md#second-order-motion)). There the motion error is the force error: the part from A(Ω), B(Ω), M and C is at most 1.7% of the scale. Two bodies and the free-surface wave of the motion are covered by unit tests only. Quadratic damping (`BQDR`) is not used; see [Second-order motion](#second-order-motion).

### QTFs, quadratic part

With `qtf` set, the near-field terms are evaluated for every pair of frequencies (ω_i, ω_j) of `ControlFile.in`, not only for ω_i = ω_j. Every term of the mean drift is a product of two first-order quantities, so the same formulas give the difference-frequency part (with the conjugate of the second factor) and the sum-frequency part (without it); the diagonal of `QTF_Quadratic.10d` is the mean drift of `MeanDriftPressure.9`. The definition, the conjugation and the scaling are in [Output Files](output_files.md#wamit-format). The zero and infinite frequencies of `0_inf_frequency_limits = 1` are left out.

- Use a uniform frequency list fine enough for the slow-drift periods of interest: the QTF is only known on the grid.
- The same mesh guidance as for near-field mean drift applies, with `velocity_method = source` on hulls with sharp edges, as for the mean drift. On a sharp cylinder at 448 panels both QTFs are within about one point of Nemoh's quadratic part, and their median error against an exact reference is 0.1% to 9% ([validation](second_order_validation.md#quadratic-qtfs)).
- `pre_postprocessing_python/SecondOrder/qtf_tools.py` reads the files, gives Newman's approximation from the diagonal (Newman 1974) and the spectrum of the slowly varying force.

### QTFs, potential part

With `qtf_terms = body` (and `qtf = 1` or `3`), the difference QTF also gets the body forcing of the second-order potential, by the indirect method (Molin 1979): the second-order body condition on the mean surface (Ogilvie 1983) and the second-order incident wave (Sharma and Dean 1981) are integrated against the radiation potentials ψ_k at the difference frequency Ω = ω_i − ω_j. `QTF_Indirect.11d` holds the quadratic part plus these terms, with the rows of `QTF_Quadratic.10d` ([Output Files](output_files.md#wamit-format)). With `body` the free-surface integral of the second-order potential is left out, so the result is the analogue of Nemoh's DUOK + HASBO, not of a complete QTF. `qtf_terms = full` adds it (finite depth only): the integral of ψ_k times the second-order free-surface forcing over z = 0 outside the bodies, on an inner mesh out to a circle R and by partial waves beyond (Molin 1979; Kim and Yue 1989); see [Developer Documentation](developer_documentation.md#free-surface-integral-qtf_terms--full). The formulation is in [Developer Documentation](developer_documentation.md#potential-part-body-forcing).

**Indirect sum frequency (`qtf = 2` or `3`, Phase 6b).** With `qtf_terms = body` or `full` the sum-frequency QTF also has the indirect potential part T⁺_pot: the load of φ_I2 and the body-forcing terms (and with `full` the free-surface integral) against the radiation potentials ψ at the sum frequencies Ω = ω_i + ω_j, which equals the direct one by Green's second identity at Ω > 0 ([Developer Documentation](developer_documentation.md#sum-frequency-of-the-direct-method)). `QTF_Indirect.11s` is `QTF_Quadratic.10s` plus T⁺_pot ([Output files](output_files.md#wamit-format)). It needs no direct solve: the default `qtf_direct = 0` is enough, because `sum_on` in `qtf_init` ([QTFDriver.f90](../src/QTFDriver.f90)) has no `qtf_direct` term, and the check that stopped `qtf = 2` with `qtf_direct = 0` was removed from `SecondOrderInput.f90`.

- **Use `qtf_terms = full`.** With `body` the sum has no free-surface integral. The log warns (`QTF warning: qtf_terms = body with the sum frequency leaves out the free-surface integral`, `sum_warnings` in `QTFDriver.f90`), as WAMIT advises against a sum-frequency QTF without it.
- **Deep water.** `qtf_terms = full` with `H < 0` stops, for the sum as for the difference: `qtf_terms = full needs a finite water depth (the outer free-surface integral uses the finite-depth modes); set a large finite depth instead.` (`SecondOrderInput.f90`). The sum with `body` has no such stop.
- **Cutoff.** The sum stops at the first irregular frequency, as `.12s` does: with `If_remove_irr_freq = 0` the pairs with ω_i + ω_j at or above the estimate of the lowest irregular frequency, and with `qtf_max_sw` those above it, have no row in `.11s` (`sum_pair_index` = 0); `.10s` keeps them. The log names the irregular-frequency cutoff and the number of pairs it leaves out (nothing is logged for pairs cut by `qtf_max_sw`), and with `qtf_direct = 0` says they are left out of `QTF_Indirect.11s` (`sum_warnings`).
- **Frequencies and cost.** ψ and a radiation-only solve are needed at every distinct ω_i + ω_j kept. With `qtf = 2` and `qtf_direct = 0` these are the only extra solves, and no direct solve runs (`qtf_radiation_omegas`, `RadiationOnlySolves` in `HAMS_Prog.f90`). On a grid of N frequencies that is up to N(N + 1)/2 more Ω, which `qtf_max_sw` limits (counts in the sum paragraph of the direct method below). The pair loop of T⁺_pot is one more indirect pair loop (`potential_part_sum`).
- **Memory.** The store adds T⁺_pot, 16 · 6 · (bodies) · (frequencies)² · (heading pairs) bytes; with `qtf_direct = 0` that is the only sum store, and with 1 T⁺_dir comes too (`nsum_pairs`, `SecondOrderDriver.f90`). The log gives the size.
- **Free-surface mesh.** With `full` the generated inner mesh and M follow the largest sum frequency; see `qtf_fs_panel` and `qtf_fs_modes`.
- **Output.** `QTF_Indirect.11s` always, and `QTF_PotentialTermsSum.txt` with `write_terms` 1 or 2.

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

**Memory.** The store also holds ψ at every Ω (6 × the number of bodies complex values per panel), the potential part of every pair, and with `qtf_b2 = control` the fields at the points of C and F. It is allocated at set-up, and the log gives its size: `QTF: with psi and the potential part the store takes <size> MB`, with the share of C and F for `qtf_b2 = control`. Measured on the sharp cylinder with one heading and 101 frequencies, with `qtf_b2 = body` (no C and F share): 8.1 MB at 448 panels and 29.4 MB at 1792. With the default `qtf_b2 = control` the store is larger by the fields on C and F. With `qtf_direct = 1` it also holds T_dir and, with `write_terms` 1 or 2, the two parts of it (see the direct method below).

**B2 and the control surface.** With the default `qtf_b2 = control`, each run also takes one field pass over the points of C and F per kept frequency (the total velocity, all headings) and one per difference frequency (∇ψ_k of every mode). The pass costs about N_hull × N_C Green function evaluations [derived: one influence row of the hull per point]. With the default `cs_*` keys, the margins and the panel size of C come from the median hull panel size, so C moves in and is refined with the hull, and this cost grows faster than the first-order solve. Measured on the free sharp cylinder (Snellius, 32 threads, 101 frequencies, 30 difference frequencies): the B2 pass took 6.5 s per frequency at 448 panels (1416 points on C and F), 69 s at 1792 (3212 points) and 288 s at 4032 (5648 points). The whole run took 15 min 34 s at 448 panels and 2 h 49 min at 1792; at 4032 it took 11 h 52 min, and with C held fixed 3 h 58 min.

- **Keep the defaults** for accuracy. B2 converges with the hull mesh and C together ([validation](second_order_validation.md#qtf-potential-part)).
- **To limit the cost** on a fine hull, hold C fixed: set `cs_margin`, `cs_bottom_margin` and `cs_panel_size` explicitly, for instance to the values the default gives on a coarser mesh. On the free sharp cylinder at 1792 panels, C and F of the 448-panel default (7.84, 7.84 and 1.96 m) took the run from 2 h 49 min to 1 h 18 min, with B2 within 6.9 / 13.8 / 6.9% (F_x / F_z / M_y) of the exact reference against 7.3 / 14.2 / 7.1% with the default. With the same C at 4032 panels, B2 kept converging (4.1 / 9.9 / 5.3%, and 0.82 / 1.63 / 1.09% away from the resonance). The default C at 4032 panels gives B2 within 4.75 / 10.42 / 5.40% and took 11 h 52 min, so the fixed C takes about a third of the time for the same accuracy. This is an option for cost, not needed for accuracy.
- **Fixed bodies** (`motion_state = fixed`): B1, B2 and B3 are zero, since each is a product with the motion [derived: the formulation in [Developer Documentation](developer_documentation.md#potential-part-body-forcing)]. `qtf_b2 = body` then gives the same result without the passes over C and F [derived].
- `qtf_edge_weights` stays `off`: see [what was tried](second_order_validation.md#tried-and-not-kept).

**Lightly damped resonances.** A free body with radiation damping only can have a sharp resonance. On the free sharp cylinder the coupled surge-pitch mode at 0.77 rad/s has a peak of about 11 m/m in surge. There the computed surge and pitch are 3% to 5.5% off at 1792 panels, and the QTF, which is quadratic in the motions, inherits about twice that at the pairs with one frequency in the peak ([validation](second_order_validation.md#the-resonance-and-the-motions)). The second-order terms themselves are within 1% to 2% away from it. Where viscous damping matters, as at a sharp bilge, add it as external linear damping in `Hydrostatic.in` (lines 11 to 16): it lowers the peak and with it this error. Measured on the free sharp cylinder with external damping in heave and pitch of 1.4 and 2.1 times the radiation damping at the 0.77 rad/s peak: the potential part was within about 3% of the exact reference over all pairs at 1792 panels, against up to 13% without the damping ([validation](second_order_validation.md#damped-free-body)). Model tests of a floating wind turbine also found its low-frequency resonant response sensitive to viscous damping (Zhang et al. 2026).

**Threads.** The field pass over C and F runs its loop over the points on `NTHREAD` threads. The total velocity comes from `fields_green`, or from `fields_source` with `cs_velocity = source` (both in [FieldVelocity.f90](../src/FieldVelocity.f90)); ∇ψ_k always comes from `fields_green`. In `fields_source` the loop over the panels is serial. In `fields_green` the loop over the panels inside it (`field_influence_row`, [WaterlineField.f90](../src/WaterlineField.f90)) is also an OpenMP loop, which runs on one thread there as long as nested parallelism is off, the OpenMP default. If nested parallelism is switched on (for instance `OMP_MAX_ACTIVE_LEVELS` above 1), up to `NTHREAD` × `NTHREAD` threads are started, and nothing limits that. Leave nested parallelism off.

**Direct method (`qtf_direct = 1`).** The potential part can also be had by solving the second-order boundary-value problem for the diffracted and radiated potential φ_D at each difference frequency, with the body forcing (`qtf_terms = body`) or with the body forcing and the free-surface forcing (`qtf_terms = full`, below), and integrating φ_I2 + φ_D over the hull. It is written to `QTF_Direct.12d` ([Output Files](output_files.md#wamit-format)) beside `QTF_Indirect.11d`, which stays. It is a check of the indirect result. It is not a more accurate alternative on a sharp hull: both methods integrate B2 on the hull there, with different discretisations. The two loads are equal by Green's theorem to the mesh error ([Developer Documentation](developer_documentation.md#direct-method)).

```
mean_drift     = 1
drift_methods  = far near
velocity_method = source
qtf            = 1
qtf_terms      = body
qtf_b2         = body
qtf_max_dw     = 0.3
qtf_direct     = 1
```

- **Compare with `qtf_b2 = body`.** The direct method forms B2 by the Stokes form with the Green function and does not read `qtf_b2`. `.12d` equals `.11d` only when `.11d` uses the hull form, so set `qtf_b2 = body` for the comparison. With the default `qtf_b2 = control`, `.11d` has the control-surface B2 and differs from `.12d` by the difference of the two forms of B2.
- **Cost.** The solve runs after the frequency loop, at every difference frequency within `qtf_max_dw`, including those on the main grid (`RadiationOnlySolves` then covers every Ω). A difference frequency that is not on the main grid has its matrix from the radiation-only solve that ψ needs anyway, so the direct method adds only the right sides and one block solve for all its pairs. One that is on the main grid costs one more matrix assembly and factorisation. The cost per pair grows with the square of the number of panels. With `ISOLV = 2` every Ω also gets an O(N³) LU factorisation, on top of the iterative solve. Measured: 7.2 ms per pair on a laptop with 8 threads on the 448-panel sharp cylinder, 5% of that run (5130 pairs); the 4032-panel sharp run on Snellius (24 CPUs, 101 frequencies, 30 difference frequencies) took 2 h 26 min in all ([validation](second_order_validation.md#direct-method-against-indirect-method)).
- **Memory.** The direct method stores T_dir, 16 · 6 · (bodies) · (frequencies)² · (heading pairs) bytes, the size of the potential part of the pairs. With `write_terms` 1 or 2 it stores two more arrays of that size (`tdir_i2` and `tdir_d`, for `QTF_DirectTerms.txt`). The logged store size counts the potential part and T_dir, not these two. Each array is 1.0 MB for one body, 101 frequencies and one heading, so it only matters for many frequencies and heading pairs. The work arrays of the solve are about 25 complex values per panel and pair, in blocks of at most 256 MB.
- **Irregular frequencies.** `If_remove_irr_freq = 1` works with the direct method, if the mesh has lid points; without them the run stops. On the sharp cylinder the lid removes the spike of the direct-indirect difference at the irregular frequency ([validation](second_order_validation.md#irregular-frequency-removal)).
- **Accuracy.** On the rounded hull the direct and indirect loads agree to 0.23% (forces) and 0.35% (moments) at 4320 panels. On the sharp cylinder with the default meshes the direct load is closer to the exact potential part than the indirect one with `qtf_b2 = body` (F_x / F_z / M_y 13 / 13 / 28% against 21 / 21 / 49% at 4032 panels, largest over the pairs, resonance included). The default indirect method, `qtf_b2 = control`, is closer still: 6.9 / 9.1 / 5.8% on the same pairs ([validation](second_order_validation.md#free-body-the-hull-and-control-surface-forms-of-b2)). For accuracy on a sharp hull use the default `.11d`, and `.12d` as an independent check.
- **Sum frequency.** `qtf` 2 or 3 with `qtf_direct = 1` also gives the sum-frequency QTF: see the paragraph below.
- **Sum frequency (`qtf = 2` or `3` with `qtf_direct = 1`, Phase 7c-2).** The direct solve also runs at the sum frequencies Ω = ω_i + ω_j, with the same matrix and the sum form of the second-order forcing, and `QTF_Direct.12s` is `QTF_Quadratic.10s` plus the load of φ_I2 + φ_D at Ω. By Green's second identity at Ω > 0 this load equals the indirect sum-frequency potential part [derived, [Developer Documentation](developer_documentation.md#sum-frequency-of-the-direct-method)]. The run also computes the indirect sum load of the same pairs from the same ψ at Ω, writes it as `QTF_Indirect.11s` (see the paragraph on the indirect sum frequency above), and with `write_terms` 1 or 2 also beside the direct parts in `QTF_DirectTermsSum.txt`: that is the pair-by-pair check of `.12s`.
  - **Use `qtf_terms = full`.** With `body` the sum has no free-surface integral. WAMIT advises against using a sum-frequency QTF without it, and the log warns. `full` needs finite depth, as for the difference.
  - **Frequencies.** ψ and a radiation-only solve are needed at every distinct ω_i + ω_j kept, besides the difference frequencies; on a grid of N frequencies that is up to N(N + 1)/2 more Ω, which `qtf_max_sw` limits. A frequency that is both a sum and a difference is solved once. On a uniform grid of 101 frequencies (the grid of the Phase 7b cases) the sums are 201 Ω and 5151 pairs per heading pair [counted].
  - **Pairs.** T⁺_ji = T⁺_ij (see `.10s`), so one of each pair is solved and both rows are written. The diagonal ω_i = ω_j is a regular pair at Ω = 2ω_i; its difference counterpart is the mean drift and has no potential part.
  - **Irregular frequencies.** The sums reach 2ω_max, which can pass the lowest irregular frequency of the hull. With `If_remove_irr_freq = 1` and lid points every sum is kept. With 0, the sums at or above the estimate of `qtf_max_sw` are left out of `.12s`, and the log names the cutoff and the number of pairs.
  - **Free-surface mesh.** With a generated inner mesh the panel size follows the shortest wave of the sum, λ(Ω_max)/10, so the inner mesh and the setup per Ω are much larger than for the difference alone; the log gives the numbers (`QTF direct free surface: ...`). With `qtf_fs_panel` or a `FSInnerMesh.pnl` the mesh is the one of the difference run and the sum is less accurate (about 1.5% of max |T⁺_pot| at Ω = 1.4 rad/s on the 448-panel cylinder, against 0.6% with the mesh built for the sum; test `direct_fs_sum_equals_indirect`) [measured]. Against the exact sum QTF of the sharp truncated cylinder (sharp m2, fixed and free), the free-surface part of `.12s` is within 2.0% (max) and 0.5% (median) of the exact one at every sum frequency up to 2.8 rad/s. This needs the near-waterline fit of the inner fields (commit 47d5911), which is always on; before it the error grew with Ω⁺ to a median of 17% at 2.8 rad/s. The log line `QTF free surface: N inner points within one waterline segment of a waterline take the fit ...` gives the points fitted and any left on the stencil because another waterline is close ([validation](second_order_validation.md#snellius-validation-of-the-sum-frequency)).
  - **Memory.** The store adds T⁺_dir and T⁺_pot, 16 · 6 · (bodies) · (frequencies)² · (heading pairs) bytes each, and with `write_terms` 1 or 2 two more (three with `full`); the logged size counts the sum pairs.
  - **Field.** With `qtf_field = 1` the sum rows of `QTF_Field.txt` get the direct part, i ρ Ω φ_D and (i Ω/g) φ_D with Ω = ω_i + ω_j, body forcing only as for the difference. The memory of the sum rows is given in the field section.
  - **Accuracy, as far as measured.** On the laptop the direct sum agrees with the indirect sum of the same run to the mesh error of the hull: on the sharp cylinder the gap, as a share of max |P_BODY|, is 8.8 / 5.9 / 4.5% at 112 / 448 / 1008 panels (the same rate as the difference, so a mesh error at the sharp edge); on the rounded hull with a coarse bottom it is 12.7% in the moments (1.6% in the forces) on the laptop, and 9.8 / 2.7 / 1.2% (max, all pairs) at 480 / 1920 / 4320 panels on Snellius, so it is mesh error. Against the exact sum-frequency reference of the sharp cylinder, `.12s` is within 5.4 / 5.0 / 15.3% (F_x / F_z / M_y, max, about 1% median) at 4032 panels, away from the surge-pitch resonance; with the near-waterline fit the free-surface part is within 2% (max) of the exact one and the fixed-body `.12s` within 4.4% (sharp m2; Snellius, 2026-11-18). See the [validation](second_order_validation.md#direct-method-at-the-sum-frequency).
- **Free-surface forcing (`qtf_terms = full`).** The free surface z = 0 outside the waterlines forces φ_D through Q_D, the quadratic terms of the second-order free-surface condition (the same Q_D as in the indirect method). With `full` the solve adds the term (1/g) ∬ G Q_D dS over the free surface to its right side, at every hull panel and lid point, so `.12d` is `.10d` plus the whole load of φ_I2 + φ_D, and equals `.11d` of `qtf_terms = full` (with `qtf_b2 = body`) to the mesh error ([Developer Documentation](developer_documentation.md#direct-method)). The free surface is the one of the indirect method: the inner mesh to the circle R (`qtf_fs_*` keys, or `FSInnerMesh.pnl`) and, beyond R, partial waves. It needs a finite depth, `qtf = 1` or `3`, and the memory and runtime of both `full` and the direct solve. Use it to check `.11d` of `full` with an independent route, or when the free-surface part of φ_D is wanted on its own: with `write_terms` 1 or 2 `QTF_DirectTerms.txt` gives it as the third part D_FS. The setup per Ω builds the Green function from every hull panel and lid point to every inner panel, and its Graf amplitudes for the outer part; it stores 16 · (rows) · (inner panels + (2M + 1)(N_ev + 1)) bytes for them, which the log gives (`QTF direct free surface: ... MB`), and the near-panel integrals cost up to 64 Green function calls per near pair, so on a production hull the setup is under 20 s per Ω at 4032 rows (1.3 s at 448 rows, 24 threads, Snellius) and takes 0.8 GB at 4032 rows and 360 inner panels. Per pair the outer moments of Q_D take 0.44 s on average (50 to 220 s for the 142 to 200 pairs of one Ω on the sharp cylinder, whatever the hull mesh; 1.9 to 2.5 s for two pairs on the 480-panel hull in 100 m of water, 8 threads, laptop) and the contraction under 0.2 s. The number of evanescent modes grows like h/(R − r_x), so a small hull in deep water, or a `qtf_fs_radius` close to the hull, needs far more memory: the log gives `QTF direct free surface: orders M, N_ev up to ..., about ... MB per Omega` before the frequency loop, and the run stops with a message if the arrays do not fit. The test cases run in 35 to 71 s on a quiet laptop with 8 threads (up to 134 s under load; 480 panels, depth 20 m). The cost per Ω is measured on Snellius in the [validation](second_order_validation.md#direct-method-with-the-free-surface-forcing). With field points (`qtf_field = 1`) the free-surface part of φ_D is solved apart from the body part, so the direct columns of `QTF_Field.txt` stay the body part (φ_I2 + the body-forced φ_D, as its header says) and `.12d` still has the whole load.

Sources: Chen, X.-B. (2007). Middle-field formulation for the computation of wave-drift loads. J. Eng. Math. 59, 61-82. Newman, J. N. (1974). Second-order, slowly-varying forces on vessels in irregular waves. Proc. Int. Symp. Dynamics of Marine Vehicles and Structures in Waves, London, 182-186. Molin, B. (1979). Second-order diffraction loads upon three-dimensional bodies. Applied Ocean Research 1(4), 197-202. Ogilvie, T. F. (1983). Second-order hydrodynamic effects on ocean platforms. International Workshop on Ship and Platform Motions, Berkeley. Sharma, J. N., Dean, R. G. (1981). Second-order directional seas and associated wave forces. SPE Journal 21(1), 129-140. Cong, P., Teng, B., Chen, L., Gou, Y. (2020). A novel solution to the second-order wave radiation force on an oscillating truncated cylinder based on the application of control surfaces. Ocean Engineering 204, 107278. WAMIT v7.5 user manual, ch. 15. Zhang, Y., Tosdevin, T., Tang, T., Stansby, P. (2026). Second-order difference-frequency effects on floating wind turbine responses: full-QTF modelling and experimental comparison. arXiv:2610.27940.

### Second-order motion

With `qtf_motion = 1` each free body moves under the second-order force of every row of the QTF files. Per second-order frequency Ω the run solves

[−Ω²(M + A(Ω)) − iΩ(B(Ω) + B_ext) + C + K_ext] ξ⁽²⁾ = F⁽²⁾,

with M = `MATX`, B_ext = `BLNR`, C = `CRS` and K_ext = `KSTF` of `Hydrostatic.in`, and F⁽²⁾ the force of the row. Ω is ω_i − ω_j for the difference and ω_i + ω_j for the sum. A(Ω) and B(Ω) are the added mass and radiation damping at |Ω|, computed from the radiation potentials ψ that the potential part already solves at every Ω, so the motion needs no extra solve. The force of `QTF_Indirect.11d` and `.11s` gives `QTF_Motion.11d` and `.11s`; with `qtf_direct = 1` the force of `QTF_Direct.12d` and `.12s` also gives `QTF_Motion.12d` and `.12s` ([output files](output_files.md#wamit-format)). The derivation and the time convention are in [Developer Documentation](developer_documentation.md#second-order-motion).

```
mean_drift      = 1
drift_methods   = far near
velocity_method = source
qtf             = 3
qtf_terms       = full
qtf_direct      = 1
qtf_motion      = 1
```

- **The static diagonal and mooring.** On the difference diagonal (ω_i = ω_j) Ω = 0 and the equation is static: (C + K_ext) ξ⁽²⁾ = F̄, the offset under the mean drift force. For a body without mooring, C + K_ext has no restoring in surge, sway and yaw and is singular. The run checks the matrix of each Ω by its reciprocal condition number rcond (LAPACK `ZGECON`, 1-norm). If rcond < 1e-12, the pairs of that Ω are left out of the motion files, the file header names them (`# Left out: ...`), and the log has one line, for an unmoored body:
  ```
   QTF motion: Omega =   0.00000 rad/s, rcond 0.000E+00 < 1e-12 (modes 1 2 6 without restoring); its pairs are left out of the QTF_Motion files.
  ```
  Every other Ω gets one line ` QTF motion: Omega = ... rad/s, rcond ...` (`qtf_motion_part`, [QTFDriver.f90](../src/QTFDriver.f90)). To get the static diagonal, give K_ext in surge, sway and yaw; heave, roll and pitch have the hydrostatic C. A spring in surge alone leaves sway and yaw free, so Ω = 0 is still left out, surge included (integration test `motion_integration_rows`: modes 2 6 named). An Ω > 0 with a small but regular matrix, such as a lightly damped slow-drift resonance, is written as computed.
- **Example: the validation case.** The free sharp truncated cylinder (radius 10 m, draft 10 m, depth 100 m) gets springs for a natural period of 100 s and dampers of 5% of critical in surge, sway and yaw, with ω_n = 2π/100 = 0.0628 rad/s: K_11 = K_22 = 1.96063e4 N/m and B_11 = B_22 = 3.12044e4 N s/m (from M + A_11), K_66 = I_66 ω_n² = 6.20126e5 N m and B_66 = 0.1 I_66 ω_n = 9.86960e5 N m s (A_66 = 0 for an axisymmetric hull). The derivation is beside `SHARP_MOTION_K11` in `pre_postprocessing_python/SecondOrder/drift_convergence.py`. In `Hydrostatic.in` they are the (1,1), (2,2) and (6,6) entries of the external restoring and external linear damping matrices. The slow-drift resonance is then at Ω ≈ 0.063 rad/s; the validation grid adds 1.013 rad/s to 0.40–1.40 rad/s so that the pair (0.95, 1.013) has Ω = 0.063 rad/s. With these springs every Ω, Ω = 0 included, is in the files, and the smallest rcond was 2.6e-5 ([validation](second_order_validation.md#second-order-motion)).
- **Damping: give low-frequency damping as `BLNR`.** The motion takes B(Ω) + `BLNR` only. `BQDR` is a linearisation on the first-order motion amplitude, iterated in `SolveMotion`, and has no second-order form, so the motion does not use it. If any entry of `BQDR` is not zero, the log warns once (`body_terms_setup`, [SecondOrderDriver.f90](../src/SecondOrderDriver.f90)):
  ```
   QTF warning: the quadratic damping BQDR is not used by the second-order motion (qtf_motion = 1),
   which has no second-order form of it; the motion takes B(Omega) + BLNR. Give the low-frequency
   damping of the mooring or the hull as BLNR if it is needed there.
  ```
  The damping sets the height of a slow-drift peak: in the validation case the peak at Ω = 0.063 rad/s is about 10 times the static surge, as 1/(2ζ) with ζ = 5% from `BLNR` predicts ([validation](second_order_validation.md#second-order-motion)). Without `BLNR` the peak is set by B(Ω) alone.
- **Sum frequency.** `QTF_Motion.11s` has the rows of `.11s`. With `qtf_terms = body` the sum warning adds ` The same holds for QTF_Motion.11s, the motion from that force.`, and with `If_remove_irr_freq = 0` the cutoff warning adds ` QTF_Motion.11s leaves out the same pairs.` (`sum_warnings`, [QTFDriver.f90](../src/QTFDriver.f90)).
- **Several bodies.** M, C, K_ext and B_ext are block-diagonal, from the `Hydrostatic_i.in` file of each body; A(Ω) and B(Ω) couple the bodies (unit test `motion_two_bodies_coupled`).
- **The radiated wave in the field.** With `qtf_field = 1` and `qtf_direct = 1`, `QTF_Field.txt` adds the wave radiated by the motion of `QTF_Motion.12d` and `.12s`, φ_R = Σ_j (−iΩ ξ_j) ψ_j, to `Re_dir` and `Im_dir`, with p = iρΩ φ_R and η = (iΩ/g) φ_R as for φ_D. The column `has_rad` marks the rows that have it: 0 on rows without a direct part and on rows of an Ω left out of the motion. With `qtf_direct = 0` there is no second-order potential in the field, so no radiated wave is added (`qfield_radiated_active`, [QTFField.f90](../src/QTFField.f90)): adding a complete term to the partial field of `QTF_FieldPartial.txt` would mix the two.
- **Cost and memory.** Per Ω the run forms A and B from ψ and factors one 6N × 6N matrix (N bodies); this is small next to the radiation solves. The store adds ξ, 16 · 6N · (frequencies)² · (heading pairs) · (routes) bytes for the difference and as much again for the sum, with one route without `qtf_direct` and two with it; the log gives it (` QTF motion: of which ... MB for the second-order motions (qtf_motion = 1)`). With the radiated wave in the field, each Ω also takes one field pass of the 6N ψ columns over the points; its time is in the field time of the `QTF direct` log line, and ψ at the points (16 · (points) · 6N bytes per Ω) is not in the logged size.
- **Accuracy.** On the validation case the motion error is the force error: the part of the error from A(Ω), B(Ω) and the polygon's M and C is at most 1.7% of the scale at 448 panels and 0.03% at 1792 ([validation](second_order_validation.md#second-order-motion)). Converge the force first.

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

with φ⁽²⁾ the second-order incident potential. The formulas per frequency pair are in [Developer Documentation](developer_documentation.md#partial-second-order-field).

**Points.** As for `PressureElevation.6p`, a point with |z| ≤ 1e-6 is an elevation point (its z is set to 0) and a point with z < −1e-6 a pressure point. The run stops with a message naming the point and the body when a point is:
- above the free surface (z > 1e-6);
- on a hull: within one panel size of the centroid of a hull panel, with the size of that panel;
- inside a body: the hull panels (symmetry images included) and their mirror image in z = 0 enclose the point, judged by the winding number (the solid angle they subtend at the point, over 4π, is more than 0.5). This also rejects a point on z = 0 inside the waterplane of a surface-piercing body. The test needs a watertight hull.

**Route and cost.** The first-order potential and its gradient at the points always come from the Green function, whatever `velocity_method` or `cs_velocity` is, with central differences of step δ = 0.02 times the size of the nearest hull panel (the default of `cs_fd_step`). At each frequency of the QTF grid, each stencil point takes one influence row (the Green function and its normal derivative from every hull panel), which serves all headings and is not stored. A pressure point takes 7 rows: the point and x ± δ, y ± δ, z ± δ. An elevation point takes 5: the point and x ± δ, y ± δ, with ∂φ/∂z from the linear free-surface condition. A pressure point less than δ below z = 0 also takes 5 rows. On a 480-panel hull, 6 points (3 of each type) and 2 headings took 1.5 to 1.8 s per frequency on a loaded laptop. The log gives this time at each frequency.

- The store keeps φ and ∇φ at each point, 64 bytes per point, frequency and heading. It is allocated at set-up, so a run that does not fit stops before the first-order solves, and the log gives its size. `qtf_field = 1` also keeps the first-order states of every heading of a frequency, as `mean_drift = 2` does.
- In finite depth above ν·h = 50 (ν = ω²/g), the log warns at each such frequency that the field values may scatter, for the reason given under `cs_velocity`.
- `pre_postprocessing_python/SecondOrder/qtf_field_plots.py` maps the file and writes grids of field points: see [Python Scripts](python_scripts.md#partial-second-order-field-maps).

**Direct part (`qtf_field = 1` with `qtf_direct = 1`).** With both keys (`qtf_direct` needs `qtf` 1, 2 or 3), HAMS-MREL also writes `QTF_Field.txt`. Its 15 first columns are those of `QTF_FieldPartial.txt`, which is written as without `qtf_direct`, byte for byte. Three columns follow: `Re_dir`, `Im_dir` and `has_direct` (columns 16 to 18), and with `qtf_motion = 1` a fourth, `has_rad` (column 19; see [the radiated wave](#second-order-motion) below). The columns are in [Output Files](output_files.md#wamit-format).
- **What the direct part is.** The potential φ_D of the direct method ([above](#qtfs-potential-part)), solved at each difference frequency Ω = ω_i − ω_j (rows K = 1) or sum frequency Ω = ω_i + ω_j (rows K = 2), evaluated at the field points by Green's second identity over the hull: 4π φ_D(x) = ∫_S [G ∂φ_D/∂n − φ_D ∂G/∂n] dS, with the body forcing as the right side of ∂φ_D/∂n. It enters as the second-order incident potential does (Bernoulli's equation above): p = iρΩ φ_D and η = (iΩ/g) φ_D, per unit amplitude pair, with φ_D in place of φ⁽²⁾ [derived; same forms as the incident part in `QTFField.f90`].
- **It is not added to Re and Im.** The columns Re, Im, Re_quad, Im_quad, Re_inc, Im_inc are the partial value of `QTF_FieldPartial.txt`. The value with the direct part is Re + Re_dir and Im + Im_dir.
- **Rows with a direct part.** `has_direct` = 1 on the difference rows (K = 1) with ω_i ≠ ω_j for the pairs the direct solve covers (those within `qtf_max_dw`), and on the sum rows (K = 2) of the pairs of `QTF_Direct.12s` (the diagonal included, within `qtf_max_sw` and the irregular-frequency cutoff). The mean rows (ω_i = ω_j, K = 1) and the rows outside those limits have `has_direct` = 0 and `Re_dir` = `Im_dir` = 0: the direct solve has no such pair.
- **What is not included.** The free-surface forcing of φ_D is not in the direct part. With `qtf_terms = full` it is in `QTF_Direct.12d`, but its part of φ_D is solved apart and not evaluated at the field points (Green's formula over the hull alone holds only for the body-forced φ_D). The wave radiated by the second-order motion of a free body is in it only with `qtf_motion = 1`, on the rows with `has_rad` = 1 ([Second-order motion](#second-order-motion)). The sum of the columns is therefore not yet the complete second-order field, and it does not yet satisfy the second-order free-surface condition on z = 0. The check of that condition needs the free-surface forcing and is for a later phase.
- **Points.** The point rules above apply. A point on z = 0 off the hull keeps the factor 4π: the image term of G coincides with 1/r there, so the singularity is 2/r over the half space below z = 0 [derived; `QTFDirectField.f90`]. Approaching the hull, the field tends to the value on the hull panel: on the truncated cylinder it is within 1.5e-3 at 0.25 panel sizes from a panel centroid (measured in the unit tests, [validation](second_order_validation.md#direct-part-of-the-field)).
- **Cost and memory.** Per difference frequency, each point takes one influence row of the Green function over the hull and one B2 integral for each block of pairs (the pairs of a frequency are solved in blocks of columns that fit a memory budget; normally one block, so one row per point and frequency). The store keeps φ_D at each point, frequency pair and heading pair, 16 bytes each: 16 · (points) · (frequencies)² · (heading pairs) bytes, in addition to the store of the partial field. With the sum (`qtf` 2 or 3) there is a second store of the same size for the sum rows. It is allocated at set-up and the run stops with a message if it does not fit; use fewer frequencies, headings or field points, or `qtf_direct = 0`.
- **Finite depth.** The direct part at the points uses the Green function without the Limes correction, as the first-order field at the same points does, so isolated outliers of the series do not land on single points. At depth 10 on a truncated cylinder this differs from the route of the hull solve by 1.1e-5 relative in the influence rows (measured in the unit tests).

Sources: Dean, R. G. and Dalrymple, R. A. (1991). Water Wave Mechanics for Engineers and Scientists. World Scientific. Kim, M.-H. and Yue, D. K. P. (1989). The complete second-order diffraction solution for an axisymmetric body. Part 1. Monochromatic incident waves. J. Fluid Mech. 200, 235-264. Kim, M.-H. and Yue, D. K. P. (1990). Part 2. Bichromatic incident waves and body motions. J. Fluid Mech. 211, 557-593. Molin, B. (1979). Second-order diffraction loads upon three-dimensional bodies. Applied Ocean Research 1, 197-202. Sharma, J. N. and Dean, R. G. (1981). Second-order directional seas and associated wave forces. Soc. Petroleum Eng. J. 21, 129-140.
