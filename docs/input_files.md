# Input Files

HAMS-MREL requires a set of input files to run. The number files depend on whether you are analyzing a single body or multiple bodies:
- `ControlFile.in`
- `Hydrostatic[_i].in`
- `HullMesh[_i].pnl` (one per body, optionally with dipole panels for zero-thickness plates)
- `WaterPlaneMesh[_i].pnl` (only if analyzing irregular frequencies)
- `gen_mod_[i]_[dof].txt` (optional, for generalized modes)
- `SecondOrder.in` (optional, for second-order wave loads)
- `CSMesh_[i].pnl` (optional, the control surface, with `cs_shape = file`)
- `FSMesh_[i].pnl` (optional, the free surface inside the control surface)
- `FSInnerMesh.pnl` (optional, the inner free surface of the QTF free-surface integral)

For multi-body cases, `_i` denotes the body number (`i = 1 … numberofbodies`). For a single body, no index is used (`HullMesh.pnl`, `Hydrostatic.in` and `WaterPlaneMesh.pnl`). The control-surface files `CSMesh_[i].pnl` and `FSMesh_[i].pnl` always carry the body number, also for a single body.

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
| 4 | free | `NELEM, NTND, ISX, ISY[, NDIPOLE]` | Number of hull panels, number of nodes, X-symmetry flag, Y-symmetry flag and, optionally, the number of dipole panels `NDIPOLE` (see [Dipole panels](#dipole-panels-optional)). Without the fifth field `NDIPOLE` is 0 |
| 5-6 | free | — | Comment or blank lines (ignored) |
| next `NTND` | free | `M` and `XYZ_LOCAL_MULTI(bodyID,nodeID,1:3)` | Node coordinates. Each line contains: `node_ID  x_local  y_local  z_local` |
| next 3 lines | free | — | Comment or blank (ignored) |
| next `NELEM` lines | free | `M`, `NCN_MULTI(bodyID,nodeID)` and `NCON_MULTI(bodyID,nodeID,1:numvertices)`| Panel definitions. Each line contains: panel number, number of vertices and (vertex1, vertex2, ...) |
| next line (only if `NDIPOLE` > 0) | free | — | Comment line, `# Start Definition of Dipole Panels` (ignored) |
| next `NDIPOLE` lines (only if `NDIPOLE` > 0) | free | `M`, `DNCN(panelID)` and `DNCON(panelID,1:numvertices)` | Dipole panel definitions, in the same format as the hull panels. The vertex numbers refer to the same node list as the hull panels |

### Dipole panels (optional)

Thin plates, such as heave plates and fins, can be given as sheets of zero thickness instead of closed bodies. Each dipole panel carries the jump of the potential across the sheet, μ = φ(+n side) − φ(−n side), instead of a potential. This is a normal-dipole sheet (Hess and Smith 1964, Sec. 4.2; the dipole panels of the WAMIT User Manual). A sheet has no interior and no faces close together. A thin closed plate has both: its two faces lie a thickness apart, which needs very small panels to resolve, and its interior adds irregular frequencies and makes the irregular-frequency lid (`IRSP = 1`) hard to use.

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

With several bodies, each `HullMesh_i.pnl` may carry its own dipole block, in the same layout. The sheet belongs to body i: its nodes are in body i's node list and are moved by `LCS_i` with the hull nodes, its rotational normals are taken about `Reference_body_center_i`, and its pressure jump acts on body i only. In body i's modes it moves with body i, and with no other body. The added mass, damping and exciting forces then include the sheets in every block, the coupling between bodies included. The blocks are read by `ReadBodyMeshMulti` ([ReadPanelMeshMulti.f90](../src/ReadPanelMeshMulti.f90)) and set up by `DIPOLE_SETUP_MULTI` ([DipoleMulti.f90](../src/DipoleMulti.f90)). The checks below run on the sheets of all bodies against the hulls of all bodies, so a sheet of one body may not coincide with or overlap a hull panel of another.

The run stops at input if:

- any node of a dipole panel is less than one panel size below z = 0. Sheets must be fully submerged, and the wave part of the dipole equations is tested only for sheets deeper than their panel size;
- a dipole panel lies in the plane of a hull panel and overlaps it (the centroid of either inside the other). A sheet cannot close a body or lie on its surface. For a heave plate, the column's bottom disc is part of the hull, and the sheet is only the annulus outside it, in the same plane;
- a dipole centroid lies within 1e-6 panel sizes of a hull centroid (a sheet must not coincide with hull panels);
- the symmetry flags of the hull are not 0, `ISOLV` is not 1, field points are requested (`Number_of_field_points` above 0), or `SecondOrder.in` is active. Dipole panels support one or several bodies, no symmetry, the direct solver and first order only.

A sheet may meet the hull along an edge, for example a heave-plate annulus at the bottom of a column. Nothing is added at the junction. For a column of radius R with a heave plate of radius R_p at its bottom, z = −d: the hull panels of `HullMesh.pnl` are the column wall down to z = −d and the column's bottom disc (r < R) at z = −d; its dipole panels are only the annulus R ≤ r ≤ R_p at z = −d. The bottom disc and the annulus lie in one plane and do not overlap.

At a free edge μ goes to zero like the square root of the distance to the edge (Lamb 1932, Sec. 108, for a disc). Constant panels, collocated at their centroids, converge only as O(h) there, so the sheet mesh matters more than a hull mesh of the same size. Grade the panels toward the free rim, for example rings at r_i = R + (R_p − R) sin(π i / 2N), use fine sheet meshes, and check the result on a refined mesh. Where one number is needed, extrapolate: 2 A(fine) − A(coarse), with the fine sheet refined twice in each direction. See the [developer documentation](developer_documentation.md#dipole-zero-thickness-panels) for the method.

[Example 9](../examples/example9_single_body_dipole_heave_plate/Input/HullMesh.pnl) is a column with a heave-plate annulus in this layout.

## CSMesh.pnl (only with `cs_shape = file`)

With the control-surface mean drift and `cs_shape = file`, the control surface C of body `b` is read from `CSMesh_{b}.pnl` in the input directory. There is one file per body, numbered `1 … numberofbodies` without leading zeros, also for a single body (`CSMesh_1.pnl`). If a file is missing, the run stops with `cs_shape = file needs CSMesh_<b>.pnl in the input directory.`

The file has exactly the structure of `HullMesh.pnl` (see above), with these rules:
- `ISX = ISY = 0`. C is always given in full.
- Panels have 3 or 4 vertices.
- The coordinates are global, not body-local.
- **The panel normals must point towards the body.** HAMS takes the normal of a quadrilateral as (x4 - x2) x (x3 - x1) and of a triangle as (x1 - x2) x (x3 - x2), so the node order sets the direction. The run checks the sign of the sum over the panels of n . (x_c - x_b) A, with x_b the centre of the waterplane bounding box of the body; if it is not negative, the run stops with `CSMesh_<b>.pnl: the panel normals must point towards the body.`
- **Containment is the user's job.** C must enclose the body and no other body. The run does not test this for a file. It only checks that no point of C lies within one panel size of a hull panel centroid.
- C must cut the free surface z = 0 in a closed line. Gaps in that waterline are reported as warnings.

A generated C is written to the output directory in this format, as `CSMesh_<b>.pnl`, so it can be edited and read back.

## FSMesh.pnl (optional, with the control-surface method)

The vertical force and the roll and pitch moments of the control-surface method also need the free surface F between the hull waterline and Γc, the line where C meets z = 0 (Chen 2007, §3.1). By default F of body `b` is generated: quadrilateral strips from each hull waterline segment out to Γc, along rays from the centre of the waterplane bounding box, with `cs_fs_rows` rows graded towards the hull (the inner row is 1/4 of the outer). The rays only cover F when the waterline is star-shaped about that centre. If it is not, the run goes on without F for that body: its F_z, M_x and M_y are not written, and the log says `Warning, body <b>: the waterline is not star-shaped about its centre; give FSMesh_<b>.pnl. No free surface F, so F_z, M_x, M_y of body <b> are not computed.` The other bodies keep their F and their rows I = 3, 4, 5.

`FSMesh_{b}.pnl` in the input directory replaces the generated F of body `b`. It has the structure of `CSMesh_{b}.pnl` (`ISX = ISY = 0`, 3 or 4 vertices, global coordinates), with these rules. A given file is a request for F, so a file that cannot be read or breaks the first rule stops the run with its message:
- All nodes lie in z = 0. The run stops with `FSMesh_<b>.pnl: F must lie in the free surface z = 0.` otherwise.
- The node order does not matter: the normal of F is always (0, 0, −1), into the fluid.
- F must fill the free surface between the hull waterline and Γc. The run does not test this.

## FSInnerMesh.pnl (optional, with `qtf_terms = full`)

With `qtf_terms = full`, the inner free surface of the QTF free-surface integral, from the waterlines of all bodies out to the circle of radius R, is generated for one body whose waterline is star-shaped about its centre. `FSInnerMesh.pnl` in the input directory replaces the generated mesh. It is needed for two or more bodies, and for a waterline that is not star-shaped; without it such a run stops. There is one file for all bodies.

The file has the structure of `CSMesh_{b}.pnl` (`ISX = ISY = 0`, 3 or 4 vertices, global coordinates), with these rules:
- All nodes lie in z = 0. The run stops with `FSInnerMesh.pnl: the inner free surface must lie in z = 0.` otherwise.
- No panel centroid lies inside a waterplane.
- With `qtf_fs_radius = 0`, R is the largest distance of a node from the centre of the waterline bounding box. With `qtf_fs_radius` given, the outer nodes must lie within 0.1% of R. R must exceed the waterline and hull radius, as for the generated mesh.
- The mesh must cover the free surface between the waterlines and the circle. The log gives the share it covers, warns when it is off by more than 1% (short of R, or overlapping), and stops above 10%.

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

Why the distance matters: a lid made by clipping a square grid to the waterline circle leaves thin triangles at the wall, with centroids very close to it. Such a lid can make the damping and excitation near the irregular frequencies worse than no lid at all, while a lid of rings on the waterline nodes, with its outer centroids about half a row height from the wall, removes the irregular frequencies cleanly.


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

`SecondOrder.in` switches on the second-order wave loads: the mean drift, the quadratic transfer functions (QTFs), the second-order motion and the second-order field at the field points. The file is optional: if it is absent, or `mean_drift = 0`, HAMS-MREL computes first-order results only. It is read from the input directory, right after `ControlFile.in`, and the run stops with a message if an option is unknown, has an invalid value, or conflicts with the first-order set-up.

Each line is `key = value`. Text after `!` is a comment, blank lines are ignored, keys are case-insensitive and may appear in any order. Keys that are left out take their default.

**Mean drift**

| Key | Values | Default | Description |
| :-- | :----- | :------ | :---------- |
| `mean_drift` | `0`, `1`, `2` | `0` | 0 is off, 1 is unidirectional, 2 is all heading pairs. |
| `drift_methods` | any of `far`, `near`, `control`, separated by blanks or commas | `far` | Far-field (momentum), near-field (pressure integration) and control-surface methods. See [Which method to use](#which-method-to-use). |
| `motion_state` | `free`, `fixed` | `free` | Whether the body responds to the waves or is held fixed. |
| `velocity_method` | `differencing`, `source` | `differencing` | How the velocity on the hull is obtained for the near-field method: from the potential on the neighbouring panels, or from a source distribution solved for at each frequency. `source` cannot be used with `ISOLV = 2` or `IRSP = 1`. See [Near-field settings](#near-field-settings). |
| `waterline_method` | `extrapolate`, `field` | `extrapolate` | How the potential at the waterline is obtained for the near-field method: from the adjacent panel and its gradient, or from the boundary integral at a point close to the waterline. See [Near-field settings](#near-field-settings). |
| `waterline_offset` | positive real | `0.01` | With `waterline_method = field`: how far the point lies from the waterline, out from the body and down, each as a fraction of the waterline segment length. |
| `moment_origin` | `fixed`, `moving` | `fixed` | The point the second-order moments are taken about: `XR` of each body at its mean position (`fixed`), or `XR` moving with the body (`moving`). Applies to the mean drift of all three methods and to both QTFs. See [Moment reference point](#moment-reference-point). |
| `n_azimuth` | integer, at least 1 | `360` | Number of azimuth points in the far-field integration. |
| `far_tolerance` | real | `1d-4` | Tolerance of the far-field integration. |
| `write_terms` | `0`, `1`, `2` | `1` | 0 writes nothing extra, 1 writes the term files, 2 also writes panel diagnostics. |
| `wl_z_tolerance` | real | `1d-6` | Tolerance for panels at the waterline, relative to the body size. |

**Control surface**

These keys apply with `control` in `drift_methods`, and with `qtf_b2 = control` (the default when the QTF potential part is computed).

| Key | Values | Default | Description |
| :-- | :----- | :------ | :---------- |
| `cs_velocity` | `green`, `source` | `green` | How the velocity on the control surface is obtained: by finite differences of the Green function, or from a source distribution. `source` cannot be used with `ISOLV = 2` or `IRSP = 1`. See [Control-surface settings](#control-surface-settings). |
| `cs_shape` | `cylinder`, `box`, `file` | `cylinder` | The shape of the control surface C. `cylinder` and `box` are generated around each body; `file` reads [`CSMesh_{b}.pnl`](#csmeshpnl-only-with-cs_shape--file) for each body. |
| `cs_margin` | real, at least 0 | `0` | The horizontal gap in metres between the hull and C. 0 means 4 times the median hull panel size. |
| `cs_bottom_margin` | real, at least 0 | `0` | The gap in metres below the deepest hull point. 0 means 4 times the median hull panel size. |
| `cs_panel_size` | real, at least 0 | `0` | The panel size of C in metres. 0 means the median hull panel size. |
| `cs_fd_step` | real, above 0 and at most 0.5 | `0.02` | The finite-difference step of `cs_velocity = green` and of ∇ψ on C and F, as a fraction of the local panel size of C. Keep the default. |
| `cs_fs_rows` | integer, at least 0 | `0` | Number of rows of the free surface F between the waterline and C. 0 means (distance of Γc − distance of the waterline)/`cs_panel_size`, rounded, at least 2. Not used when [`FSMesh_{b}.pnl`](#fsmeshpnl-optional-with-the-control-surface-method) is given. |

**QTFs**

| Key | Values | Default | Description |
| :-- | :----- | :------ | :---------- |
| `qtf` | `0`, `1`, `2`, `3` | `0` | 0 off, 1 difference frequency, 2 sum frequency, 3 both. Needs `mean_drift = 1` or `2` and `near` in `drift_methods`. Heading pairs follow `mean_drift`: β1 = β2 for 1, every pair for 2. See [QTFs, quadratic part](#qtfs-quadratic-part). |
| `qtf_terms` | `quadratic`, `body`, `full` | `quadratic` | Which parts of the QTF. `quadratic`: the quadratic part only. `body`: also the body forcing of the second-order potential, by the indirect method (`QTF_Indirect.11d`, `.11s`). `full`: `body` plus the free-surface integral, the complete QTF; needs a finite depth. `body` and `full` need `qtf` 1, 2 or 3. See [QTFs, potential part](#qtfs-potential-part). |
| `qtf_max_dw` | real, at least 0 | `0` | Largest \|ω_i − ω_j\| in rad/s of the difference-frequency files. 0 writes every pair. It also limits the difference frequencies at which ψ is solved. |
| `qtf_max_sw` | real, at least 0 | `0` | Largest ω_i + ω_j in rad/s of the sum-frequency potential part (`.11s`, `.12s` and their term files). 0 keeps every pair. `QTF_Quadratic.10s` keeps every pair. See [Sum frequency](#sum-frequency). |
| `qtf_direct` | `0`, `1` | `0` | 1 also solves the second-order potential directly, at every difference and (with `qtf` 2 or 3) sum frequency, and writes `QTF_Direct.12d` and `.12s`. Needs `qtf_terms = body` or `full` and `qtf` 1, 2 or 3. See [Direct method](#direct-method). |
| `qtf_motion` | `0`, `1` | `0` | 1 computes the second-order motion of the free bodies from the second-order force and writes the `QTF_Motion` files. Needs `qtf_terms = body` or `full` and `qtf` 1, 2 or 3. See [Second-order motion](#second-order-motion). |
| `qtf_b2` | `control`, `body` | `control` | Where the surface piece of the term B2 of the potential part is integrated: on the control surface C and free surface F around each body (`control`), or on the hull (`body`). `control` needs a C with a bottom. Used with `qtf_terms = body` or `full`. See [B2 and the control surface](#b2-and-the-control-surface). |
| `qtf_edge_weights` | `off`, `on` | `off` | With `qtf_b2 = control`: product weights on the hull panels next to a sharp edge, for the hull term of B2. Keep it off. |
| `qtf_edge_angle` | real, above 0 and below 180 | `45` | With `qtf_edge_weights = on`, and for the direct method: a mesh edge is sharp when the normals of its two panels differ by more than this many degrees. |
| `qtf_field` | `0`, `1` | `0` | 1 writes the second-order pressure and elevation at the field points of `ControlFile.in` to `QTF_FieldPartial.txt` (and, with `qtf_direct = 1`, `QTF_Field.txt`). Needs `qtf` 1, 2 or 3 and at least one field point. See [Partial second-order field](#partial-second-order-field-at-the-field-points). |

**Free-surface integral**

These keys apply with `qtf_terms = full` and `qtf` 1, 2 or 3. See [Free-surface integral](#free-surface-integral).

| Key | Values | Default | Description |
| :-- | :----- | :------ | :---------- |
| `qtf_fs_radius` | real, at least 0 | `0` | The radius R in metres of the circle that splits the free surface into an inner mesh and an outer region of partial waves. 0 sets it from the size of the bodies. R must exceed the waterline and hull radius. |
| `qtf_fs_rows` | integer, at least 0 | `0` | The number of radial rows of the generated inner mesh. 0 sets it automatically. |
| `qtf_fs_panel` | real, at least 0 | `0` | The panel size in metres of the generated inner mesh. 0 sets λ/10 at the largest frequency. It must be below R − r_wl (r_wl the waterline radius). Not used with `FSInnerMesh.pnl`. |
| `qtf_fs_modes` | integer, at least 0 | `0` | The number of angular modes M of the outer region. 0 sets it automatically. |
| `qtf_fs_evanescent` | integer, at least 0 | `0` | The number N_ev of evanescent vertical modes of the outer region. 0 sets N_ev = max(8, ceil(9 h / (π (R − r_wl)))); the log gives the value used. |
| `qtf_fs_circle` | `panels`, `samples` | `panels` | How the partial-wave amplitudes outside R are found: from the hull panels by Graf's addition theorem (`panels`), or from the field sampled on the circle (`samples`). `panels` needs every panel and its images inside the circle. |
| `qtf_fs_inner` | `stencil`, `source` | `stencil` | How the fields on the inner mesh are found: by a finite-difference stencil with a near-waterline fit (`stencil`), or from source strengths (`source`). Keep `stencil`; `source` is less accurate and cannot be used with `qtf_direct = 1`, `IRSP = 1` or `ISOLV = 2`. |
| `field_deep_kh` | real, at least 0 | `10` | In finite depth, the free-surface points of the second-order passes use the deep-water Green function at frequencies with k0 h at or above this value, which is cheaper. 0 switches it off. See [Deep-water rows](#deep-water-rows). |
| `field_deep_draft` | real, at least 0 | `10` | The deep-water rows of `field_deep_kh` are used only where the depth is at least this many times the draft. 0 switches the guard off. |

Example:

```
! Mean drift with both methods
mean_drift     = 1
drift_methods  = far near
n_azimuth      = 360
```

A key that does not apply is not an error: `qtf_b2` without `qtf_terms = body` or `full` and `qtf` 1, 2 or 3, the `qtf_fs_*` keys without `qtf_terms = full` and `qtf` 1, 2 or 3, and `qtf_edge_weights` or `qtf_edge_angle` without `qtf_b2 = control` are ignored, and the log says so (`QTF note: ... ignored`).

### Checks at input

The options are checked against `ControlFile.in` and the mesh files. A run is rejected when:
- `ISOL = 2`, because the incident potential is needed.
- The water depth `H` is above 500 m. Set `H < 0` for infinite depth instead, because the incident-wave velocity routines treat depths above 500 m differently.
- `drift_methods` includes `near` or `control` and any generalized mode file (`gen_mod_{body}_{dof}.txt`) is present.
- `velocity_method = source` or `cs_velocity = source` with `ISOLV = 2` (GMRES), or with `IRSP = 1` (irregular frequency removal).
- `cs_margin`, `cs_bottom_margin`, `cs_panel_size` or `cs_fs_rows` is negative, or `cs_fd_step` is outside (0, 0.5].
- `qtf` is above 0 with `mean_drift = 0`, or without `near` in `drift_methods`; `qtf` is outside 0 to 3, or `qtf_max_dw` or `qtf_max_sw` is negative.
- `qtf_terms`, `qtf_b2`, `qtf_edge_weights`, `qtf_fs_circle` or `qtf_fs_inner` has another value than those listed (`Invalid value "<word>" for key "<key>" in SecondOrder.in.`), or `qtf_edge_angle` is not between 0 and 180.
- `qtf_terms = full` with `qtf` 1, 2 or 3 and a deep-water run (`H < 0`): `qtf_terms = full needs a finite water depth (the outer free-surface integral uses the finite-depth modes); set a large finite depth instead.`
- `field_deep_kh`, `field_deep_draft`, `qtf_fs_radius`, `qtf_fs_rows`, `qtf_fs_panel`, `qtf_fs_modes` or `qtf_fs_evanescent` is negative; `qtf_fs_radius` does not exceed the waterline or hull radius; or, with a generated inner mesh, `qtf_fs_panel` is not below R − r_wl (`qtf_fs_panel = <size> m must be below R - r_wl = <width> m`).
- `qtf_fs_inner = source` with `ISOLV = 2`, `IRSP = 1` or `qtf_direct = 1`, or with an inner mesh whose edge does not lie on the waterline nodes.
- `qtf_fs_circle = panels` with a panel (or its image) outside the circle of radius R, or with k_n R too large for the evanescent modes; set `qtf_fs_radius` larger or use `qtf_fs_circle = samples`.
- `qtf_direct` is not 0 or 1, or `qtf_direct = 1` without `qtf_terms = body` or `full` (`qtf_direct = 1 needs qtf_terms = body or full in SecondOrder.in.`) or with `qtf = 0` (`qtf_direct = 1 needs qtf = 1, 2 or 3 in SecondOrder.in.`).
- `qtf_direct = 1` with `IRSP = 1` and a mesh without lid points.
- `qtf_motion` is not 0 or 1 (`qtf_motion must be 0 or 1 in SecondOrder.in.`), or `qtf_motion = 1` without `qtf_terms = body` or `full` (`qtf_motion = 1 needs qtf_terms = body or full in SecondOrder.in.`) or with `qtf = 0` (`qtf_motion = 1 needs qtf = 1, 2 or 3 in SecondOrder.in.`). `qtf_motion = 1` with `motion_state = fixed` is accepted, and the log says that no second-order motion is computed.
- With `qtf_terms = body` or `full`, `qtf` 1, 2 or 3 and `qtf_b2 = control`: a body has no free surface F inside its C (`qtf_b2 = control: body <b> has no free surface F inside its control surface ...`), or its C has no bottom (`qtf_b2 = control: the control surface of body <b> has no bottom ...`), as for a bottom-mounted hull. Unlike the control-surface drift, which then leaves out F_z, M_x and M_y of that body, B2 has no fallback. For a generated C, set `cs_bottom_margin` below h minus the draft to keep a bottom; or set `qtf_b2 = body`.
- `qtf_field` is not 0 or 1, or is 1 with `qtf = 0` or with no field points in `ControlFile.in`; or a field point is above the free surface, on a hull or inside a body (see [the point rules](#partial-second-order-field-at-the-field-points)).
- The arrays of the run do not fit in memory. The stores are allocated at set-up, so such a run stops before the first-order solves, and the log gives their size.

### Methods and outputs

From the first-order solution, HAMS-MREL computes the mean (time-averaged) second-order loads on bodies in regular waves, by three methods:

- **Far-field** (`far`, momentum conservation): F_x, F_y and M_z of the whole system, from the far-field amplitude of the waves the bodies send out. Output `MeanDriftMomentum.8`.
- **Near-field** (`near`, pressure integration over the mean wetted surface): all six loads of each body, as the sum of six terms (waterline, flare, velocity, displacement, rotation, hydrostatic). Output `MeanDriftPressure.9`.
- **Control surface** (`control`, middle-field, Chen 2007): all six loads of each body from the field on a closed surface C around it. F_z, M_x and M_y also need the free surface between the hull and C. Output `MeanDriftControlSurface.7`.

With `mean_drift = 2`, each gives the loads of every pair of headings at the same frequency, for bidirectional seas. The terms of `near` and `control` go to `MeanDriftTerms.txt` (with `write_terms` 1 or 2). The files are described in [Output Files](output_files.md#wamit-format). When `far` runs with another method, the log prints their relative difference for F_x, F_y and M_z (M_z moved to the origin). Use it as a mesh check: a large difference means the other result has not converged.

The defaults are `velocity_method = differencing` and `waterline_method = extrapolate`. `differencing` is the closer to far-field on hulls without sharp edges. There is no switch of method by frequency.

### Which method to use

- **F_x, F_y, M_z of the whole system:** far-field is the reference on any hull.
- **Rounded hulls** (no sharp edges at the bilge or corners): near-field with the defaults, for the loads per body and for F_z, M_x, M_y.
- **Hulls with sharp edges**, such as a flat bottom meeting a vertical wall:
  - F_z, M_x, M_y: use `control`. Near-field F_z converges slowly there and can have the wrong sign on a free hull.
  - Horizontal loads per body: near-field converges slowly. `source` is the closer up to about ka = 1.5 on a cylinder, `differencing` from about ka = 2.0, and for M_y. Check either against far-field.
  - Free sharp hull at low frequency: use far-field for F_x; near-field F_x is not reliable there with either velocity method.
  - Sharp box: near-field F_y and M_z are not reliable with either velocity method.
- **Hulls with slender members through the waterline**, such as the OC4 semisubmersible with its inclined cross braces: check every hull velocity against far-field. Near-field with either velocity method, and `control` with `cs_velocity = source` (which uses the same source solve), may not meet the far-field check for F_x. Use `control` with `cs_velocity = green` and C and F panels of 1 m or less, and check its F_x against far-field before reading F_z.
- **`source`** costs about twice the run time and cannot remove irregular frequencies: `IRSP = 1` is rejected, and it is far off at an irregular frequency. Keep away from them, or use `differencing` with a lid there.
- **Waterline:** keep `extrapolate`. `field` is first order where `extrapolate` is second order; see [Limits](#limits).
- **Settings:** `GRN_KIND = 8` and `ISOLV = 2` (GMRES) change the mean drift very little against the defaults.

### Mesh guidance

- Refine towards the waterline. Thin top rows converge, because the free-surface image term is integrated over the panel ([Developer Documentation](developer_documentation.md#free-surface-image-term)).
- Round corners and the bilge where the real hull allows: near-field then converges much faster. It makes a different hull: judge it against its own far-field.
- If near-field is needed on a sharp edge, grade the rows strongly towards the edge, for example rows of size (i/n)^4 at the edge (and quadratic at the waterline) with `differencing`. This does not help `source`, and the free hull near heave resonance stays far off.
- Keep every hull row deeper than `wl_z_tolerance` times the body size below z = 0, or the run takes it for a lid panel and stops.
- Check convergence with two meshes, panels per side doubled. The rate of fall per doubling is a better check than one difference.
- Irregular frequencies: with `differencing`, use `IRSP = 1` with a lid near them.
- Control surface: start from the defaults of `cs_*`. For small hulls in a multi-body case, use a finer C (`cs_panel_size` half the default).

### Near-field settings

- **`velocity_method = source`** solves for a source distribution with the same normal derivative as the potential, and takes the velocity from its integral. The potential itself stays that of the direct solve. It costs one more dense factorization and one more set of Green function evaluations per frequency, and keeps the factorization and the gradient integrals of every symmetry block in memory: 64 nb² ns bytes, with nb panels in the base mesh and ns symmetry blocks, about as much again as the first-order arrays. The run stops if a matrix is exactly singular.
- **`waterline_method = field`** evaluates the boundary integral formula at a point `waterline_offset` segment lengths out from the body and down. The point is in the fluid if the wall at the waterline widens downward by no more than 45 degrees from the vertical; on a wider wall the run stops with a message. No other body may be closer to the point than `waterline_offset` times the segment length, so keep the gap between bodies larger than that. It costs one set of Green function evaluations per waterline segment and frequency. It converges more slowly than `extrapolate`, and in finite depth it can jump with the mesh or the offset; the run logs a warning when it is used in finite depth.

### Control-surface settings

- **`cs_velocity = green`** differentiates the potential on C by finite differences of the Green function, with the step `cs_fd_step`. In finite depth above ν·h = 50 (ν = ω²/g), the Green function still uses the series acceleration (Limes) that makes field values jump at small horizontal distances, so its values on C may scatter, and the log warns at each such frequency. There, check against `far`, or use deep water (`H < 0`) where kh is large.
- **`cs_fd_step`** above about 0.3 brings the centroids of the top row of C (a quarter panel below z = 0) within one step of z = 0. Their dφ/dz is then taken from the free-surface condition, which is off below z = 0. Keep the default.
- **`cs_fs_rows`** counts the rows from the waterline to Γc, the line where C meets z = 0, by their mean distances from the centre of the waterplane bounding box.
- The generated C of each body is written to the output directory as `CSMesh_<b>.pnl`, and can be edited and read back with `cs_shape = file`. With `qtf_b2 = control`, C and F are built and checked as for `drift_methods = control`, also when `drift_methods` has no `control`.

### Moment reference point

`moment_origin = fixed` takes the second-order moments about `XR` of each body at its mean, space-fixed position, as WAMIT writes i = 4, 5, 6 (WAMIT v7.6 manual, sec. 5.2 and eqs. 15.53 to 15.57). `moving` takes them about `XR` moving with the body, as WAMIT writes i = −4, −5, −6 and as Nemoh and HydroStar give them. The two differ by $P(\xi \times F^{(1)})$, with $F^{(1)}$ the first-order fluid force. For a free body with its centre of gravity at `XR` this is zero for the mean drift but not for the QTFs, where it can be large near a resonance. The file formats are the same, and the log says which origin was used. The `QTF_Motion` files always use `moving`, whatever this key says.

The point is always `XR` of each body, the point the first-order moments are about too; there is no separate key for it. To compare with another solver, set `XR` to that solver's reference point. See [Moment reference point](output_files.md#moment-reference-point-note) in the output files.

### Limits

- **Generalized modes:** `near` and `control` are rejected when a `gen_mod` file is present. Far-field runs with them.
- **Narrow gaps:** not tested. With `waterline_method = field`, no other body may be closer to the field point than `waterline_offset` times the segment length.
- **Multi-body:** far-field gives only the total of all bodies. Near-field and control give the loads per body.
- **`waterline_method = field` in finite depth:** term I can jump with the mesh or the offset. The finite-depth Green function uses Linton's method for R/h < 0.02, which removes the series behind such jumps for ν·h ≤ 50. The run logs a warning when this method is used in finite depth.
- **Control surface, F_z, M_x, M_y:** checked with `cs_shape = cylinder` only. Not computed when C reaches the sea bed or the free surface between the hull and C cannot be built (the log says why).
- **QTFs:** the quadratic part comes from the near-field method. With `qtf_terms = body` the QTF also has the body forcing of the second-order potential (as Nemoh's HASBO), but not the free-surface integral (Nemoh's HASFS, the free-surface part of WAMIT's indirect method; WAMIT v7.5 user manual, ch. 15), so it is incomplete. With `qtf_terms = full` (finite depth) it has the free-surface integral too and is complete; this is checked for a single body only. Multi-body and half-mesh runs of the potential part are covered by integration tests only.
- **Second-order field at the field points (`qtf_field`):** partial. The second-order diffracted wave, without its free-surface forcing, is added only with `qtf_direct = 1`, and the radiated wave of the second-order motion only with `qtf_motion = 1`. See [Partial second-order field](#partial-second-order-field-at-the-field-points).
- **Second-order motion (`qtf_motion`):** the motion is as accurate as the QTF it is solved from. Two bodies and the field of the radiated wave are covered by unit tests only. Quadratic damping (`BQDR`) is not used; see [Second-order motion](#second-order-motion).

### QTFs, quadratic part

With `qtf` set, the near-field terms are evaluated for every pair of frequencies (ω_i, ω_j) of `ControlFile.in`, not only for ω_i = ω_j. Every term of the mean drift is a product of two first-order quantities, so the same formulas give the difference-frequency part (with the conjugate of the second factor) and the sum-frequency part (without it); the diagonal of `QTF_Quadratic.10d` is the mean drift of `MeanDriftPressure.9`. The difference part is written to `QTF_Quadratic.10d` and the sum part to `QTF_Quadratic.10s`. The definition, the conjugation and the scaling are in [Output Files](output_files.md#wamit-format). The zero and infinite frequencies of `0_inf_frequency_limits = 1` are left out.

- Use a uniform frequency list fine enough for the slow-drift periods of interest: the QTF is only known on the grid.
- Choose `qtf_max_dw` by a convergence study of the slow-drift response, widening the bandwidth until the result no longer changes, as Zhang et al. (2026) do for the full difference-frequency QTF of a floating wind turbine.
- The first-order state of every frequency and heading is kept in memory, 128 bytes per panel each; the log gives the size.
- The same mesh guidance as for near-field mean drift applies, with `velocity_method = source` on hulls with sharp edges.
- `pre_postprocessing_python/SecondOrder/qtf_tools.py` reads the files, gives Newman's approximation from the diagonal (Newman 1974) and the spectrum of the slowly varying force.

### QTFs, potential part

With `qtf_terms = body`, the QTF also gets the body forcing of the second-order potential, by the indirect method (Molin 1979): the second-order body condition on the mean surface (Ogilvie 1983) and the second-order incident wave (Sharma and Dean 1981) are integrated against the radiation potentials ψ_k at the second-order frequency, Ω = ω_i − ω_j for the difference and Ω = ω_i + ω_j for the sum. `QTF_Indirect.11d` and `QTF_Indirect.11s` hold the quadratic part plus these terms, with the rows of `QTF_Quadratic.10d` and `.10s` ([Output Files](output_files.md#wamit-format)). With `body` the free-surface integral of the second-order potential is left out, so the result is the analogue of Nemoh's DUOK + HASBO, not a complete QTF. `qtf_terms = full` adds it (finite depth only): see [Free-surface integral](#free-surface-integral). The formulation is in [Developer Documentation](developer_documentation.md#potential-part-body-forcing).

Example:

```
mean_drift     = 1
drift_methods  = far near
velocity_method = source
qtf            = 1
qtf_terms      = body
qtf_max_dw     = 0.3
```

**Cost.** ψ is needed at every distinct Ω of the frequency list within `qtf_max_dw` (and, for the sum, `qtf_max_sw`). An Ω on the list takes ψ from the radiation solve of that frequency. Every other Ω costs one radiation-only solve after the frequency loop (no diffraction, no output file). The log gives the count: `QTF body forcing: psi at N difference frequencies, M from the main grid and K radiation-only solves`. On a uniform list of step Δω the distinct differences are multiples of Δω, so `qtf_max_dw` sets their number (30 for 0.01 rad/s steps and `qtf_max_dw = 0.3`); a non-uniform list can need up to N(N−1)/2 of them.

**Memory.** The store holds ψ at every Ω (6 × the number of bodies complex values per panel), the potential part of every pair, 16 · 6 · (bodies) · (frequencies)² · (heading pairs) bytes per kernel, and with `qtf_b2 = control` the fields at the points of C and F. It is allocated at set-up, and the log gives its size: `QTF: with psi and the potential part the store takes <size> MB`, with the share of C and F for `qtf_b2 = control`. For example, one cylinder with one heading and 101 frequencies, with `qtf_b2 = body`, takes about 8 MB at 448 panels and 30 MB at 1792.

**Lightly damped resonances.** A free body with radiation damping only can have a sharp resonance. Near the peak the error of the computed motions is largest, and the QTF, which is quadratic in the motions, inherits about twice that error at the pairs with one frequency in the peak. Where viscous damping matters, as at a sharp bilge, add it as external linear damping in `Hydrostatic.in` (lines 11 to 16): it lowers the peak and with it this error. Model tests of a floating wind turbine also found its low-frequency resonant response sensitive to viscous damping (Zhang et al. 2026).

**Fixed bodies** (`motion_state = fixed`): the terms B1, B2 and B3 are zero, since each is a product with the motion. `qtf_b2 = body` then gives the same result without the passes over C and F.

#### B2 and the control surface

With the default `qtf_b2 = control`, the surface piece of the second-derivative term B2 is moved from the hull onto the control surface C and the free surface F around each body by a Gauss identity (Cong et al. 2020; Chen 2007, App. C); the rest of B2 stays on the hull. On the hull B2 is inaccurate next to sharp edges. `control` is used for the difference and the sum QTF alike, so the same input gives the same `.11s` for `qtf = 2` and `qtf = 3`. C and F are built with the `cs_*` keys, and C needs a bottom: for a bottom-mounted hull set `qtf_b2 = body`.

Each run takes one field pass over the points of C and F per kept frequency (the total velocity, all headings) and one per second-order frequency (∇ψ_k of every mode). The pass costs about N_hull × N_C Green function evaluations. With the default `cs_*` keys, the margins and the panel size of C come from the median hull panel size, so C moves in and is refined with the hull, and this cost grows faster than the first-order solve. On a hull of a few thousand panels the B2 passes can take most of the run time.

- **Keep the defaults** for accuracy. B2 converges with the hull mesh and C together.
- **To limit the cost** on a fine hull, hold C fixed: set `cs_margin`, `cs_bottom_margin` and `cs_panel_size` explicitly, for instance to the values the default gives on a coarser mesh. B2 still converges with the hull mesh when C is held fixed, at a fraction of the time.
- **Threads.** The field pass runs its loop over the points on `NTHREAD` threads. With nested OpenMP parallelism switched on (for instance `OMP_MAX_ACTIVE_LEVELS` above 1), up to `NTHREAD` × `NTHREAD` threads are started, and nothing limits that. Leave nested parallelism off, the OpenMP default.

#### Free-surface integral

With `qtf_terms = full` (finite depth only, `H > 0`), the potential part also has the integral of ψ_k times the second-order free-surface forcing over z = 0 outside the bodies (Molin 1979; Kim and Yue 1989). The free surface is split by a circle of radius R (`qtf_fs_radius`) around all bodies: an inner mesh from the waterlines to R, and an outer region beyond R, where the field is written as partial waves with M angular modes (`qtf_fs_modes`) and N_ev evanescent vertical modes (`qtf_fs_evanescent`). `QTF_PotentialTerms.txt` (and `QTF_PotentialTermsSum.txt`) give the two parts as the terms `FS_inner` and `FS_outer` ([Output files](output_files.md)). The method is in [Developer Documentation](developer_documentation.md#free-surface-integral-qtf_terms--full).

- **Inner mesh.** For one body whose waterline is star-shaped about its centre, the inner mesh is generated: strips from the waterline out to a polygon on the circle, graded towards the waterline, with a column at every hull waterline node. Its panel size at the circle is `qtf_fs_panel`, or λ/10 at the largest frequency; it sets the number of chords on the circle, ceil(2πR/size) (at least 8), and the least number of radial rows. For two or more bodies, or a waterline that is not star-shaped, give the mesh as [`FSInnerMesh.pnl`](#fsinnermeshpnl-optional-with-qtf_terms--full).
- **Radius.** With `qtf_fs_radius = 0`, R = r_x + max(0.5 r_x, 2 × the largest hull panel), with r_x the larger of the waterline radius and the hull radius (the largest horizontal distance of a hull node from the centre). R must exceed both, else the run stops: a hull wider below the waterline would put the partial waves inside the body.
- **Accuracy.** The free-surface terms on the generated inner mesh get more accurate with a finer hull mesh, because the inner mesh follows the hull waterline nodes. On a coarse hull, check them on a finer hull mesh; a larger `qtf_fs_radius` or a smaller `qtf_fs_panel` is a second check (more radial rows alone, `qtf_fs_rows`, does not help). The log notes this on every run with the generated mesh.
- **Partial waves.** With `qtf_fs_circle = panels` (the default) the amplitudes of the partial waves come from the hull integral by Graf's addition theorem (Abramowitz and Stegun 9.1.79, 9.6.3 and 9.6.4), without sampling the field on the circle; a mode check uses every 8th of the circle samples. With `samples` the field is computed at every sample on the circle and decomposed in modes. Above ν h = 50 the mode check can warn on the scatter of the finite-depth field at the samples; the amplitudes do not depend on the check.
- **Inner fields.** With `qtf_fs_inner = stencil` (the default) the second derivative of the potential at the inner points comes from a finite-difference stencil, and next to the waterline from a fit. The log line `QTF free surface: N inner points within one waterline segment of a waterline take the fit ...` gives the points fitted and any left on the stencil because another waterline is close.

#### Deep-water rows

With a finite depth, `field_deep_kh` (default 10) makes the free-surface points of the second-order passes (z = 0: the free-surface fields, the waterline fit and the free-surface ψ) use the deep-water Green function instead of the finite-depth one at frequencies with k0 h at or above its value (k0 the wave number, h the depth). This makes these passes cheaper. From k0 h = 10 the deep-water rows are close to the finite-depth ones, and above ν h = 50 (ν = ω²/g) they are more accurate than the finite-depth series; at small k0 h they differ much more, so do not set the key far below the default.

- **Free surface only.** Below z = 0 the deep-water function misses the seabed images, so ∇ψ and the fields on C and F keep the finite-depth function. The hull solve, the waterline rows of the near-field drift, the field points of `qtf_field`, the direct method and, with `qtf_fs_circle = panels`, the partial waves on the circle keep it too (the log says so when k0 h of the largest frequency reaches the value). With `qtf_fs_circle = samples` the samples on the circle take the deep-water rows.
- **Depth guard.** The deep-water rows are used only where the depth is at least `field_deep_draft` (default 10) times the draft of the hull (its deepest node, all bodies together). In shallower water the free-surface passes keep the finite-depth Green function at every frequency, and the log says so once (`the finite-depth Green function is kept from k0 h = ...`).
- **Off.** `field_deep_kh = 0` uses the finite-depth function at every frequency.

#### Direct method

With `qtf_direct = 1` the potential part is also computed by solving the second-order boundary-value problem for the diffracted and radiated potential φ_D at each second-order frequency, with the body forcing (`qtf_terms = body`) or with the body and free-surface forcing (`qtf_terms = full`), and integrating φ_I2 + φ_D over the hull. It is written to `QTF_Direct.12d` and, with `qtf` 2 or 3, `QTF_Direct.12s` ([Output Files](output_files.md#wamit-format)), beside `QTF_Indirect.11d` and `.11s`, which stay. The two loads are equal by Green's theorem to the mesh error ([Developer Documentation](developer_documentation.md#direct-method)), so the direct method is a check of the indirect result. It is not a more accurate alternative on a sharp hull: for accuracy there use the default `.11d`, and `.12d` as an independent check.

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

- **Compare with `qtf_b2 = body`.** The direct method forms B2 by the Stokes form with the Green function and does not read `qtf_b2`. `.12d` equals `.11d` only when `.11d` uses the hull form, so set `qtf_b2 = body` for the comparison. With the default `qtf_b2 = control`, the two differ by the difference of the two forms of B2.
- **Cost.** The solve runs after the frequency loop, at every second-order frequency kept, including those on the main grid. An Ω that is not on the main grid has its matrix from the radiation-only solve that ψ needs anyway, so the direct method adds only the right sides and one block solve for all its pairs. An Ω on the main grid costs one more matrix assembly and factorisation. The cost per pair grows with the square of the number of panels. With `ISOLV = 2` every Ω also gets an LU factorisation, on top of the iterative solve.
- **Memory.** The direct method stores T_dir, 16 · 6 · (bodies) · (frequencies)² · (heading pairs) bytes per kernel. With `write_terms` 1 or 2 it stores two more arrays of that size (three with `full`), for `QTF_DirectTerms.txt` and `QTF_DirectTermsSum.txt`; the logged size does not count these. The work arrays of the solve are about 25 complex values per panel and pair, in blocks of at most 256 MB.
- **Irregular frequencies.** `If_remove_irr_freq = 1` works with the direct method if the mesh has lid points; without them the run stops. The lid removes the spike of the direct-indirect difference at an irregular frequency.
- **Free-surface forcing (`qtf_terms = full`).** The free surface outside the waterlines forces φ_D through Q_D, the quadratic terms of the second-order free-surface condition (the same Q_D as in the indirect method). The solve adds (1/g) ∬ G Q_D dS over the free surface to its right side, so `.12d` equals `.11d` of `qtf_terms = full` (with `qtf_b2 = body`) to the mesh error. The free surface is the one of the indirect method (the inner mesh to R and partial waves beyond), and with `write_terms` 1 or 2 `QTF_DirectTerms.txt` gives its part as D_FS. The set-up per Ω stores 16 · (rows) · (inner panels + (2M + 1)(N_ev + 1)) bytes, which the log gives (`QTF direct free surface: ... MB`); on a hull of a few thousand rows it takes seconds per Ω and under a GB. The number of evanescent modes grows like h/(R − r_x), so a small hull in deep water, or a `qtf_fs_radius` close to the hull, needs far more memory: the log gives `QTF direct free surface: orders M, N_ev up to ..., about ... MB per Omega` before the frequency loop, and the run stops with a message if the arrays do not fit.

#### Sum frequency

With `qtf` 2 or 3 and `qtf_terms = body` or `full`, the potential part is also computed at the sum frequencies Ω = ω_i + ω_j and written to `QTF_Indirect.11s`, and with `qtf_direct = 1` to `QTF_Direct.12s`. The indirect sum needs no direct solve. With `qtf_direct = 1` and `write_terms` 1 or 2, `QTF_DirectTermsSum.txt` gives the direct and indirect parts side by side, the pair-by-pair check of `.12s`.

- **Use `qtf_terms = full`.** With `body` the sum has no free-surface integral. WAMIT advises against using a sum-frequency QTF without it, and the log warns (`QTF warning: qtf_terms = body with the sum frequency leaves out the free-surface integral`). `full` needs a finite depth, as for the difference.
- **Frequencies and cost.** ψ and a radiation-only solve are needed at every distinct ω_i + ω_j kept, besides the difference frequencies. On a uniform grid of N frequencies the sums are 2N − 1 values and N(N + 1)/2 pairs per heading pair; `qtf_max_sw` limits them. A frequency that is both a sum and a difference is solved once.
- **Irregular frequencies.** The sums reach 2ω_max, which can pass the lowest irregular frequency of the hull. With `If_remove_irr_freq = 1` and lid points every sum is kept. With `If_remove_irr_freq = 0` the sums at or above an estimate of the lowest irregular frequency are left out of `.11s`, `.12s` and their term files, as are those above `qtf_max_sw`; `.10s` keeps them. The log names the cutoff and the number of pairs it leaves out. The estimate is the minimum over the closed waterline loops of √(g j₀₁ / a_eq), with j₀₁ = 2.404826 and a_eq = √(A/π), A the area of the loop: a lower bound for a wall-sided hull (Faber–Krahn; Pólya and Szegő 1951).
- **Pairs.** T⁺_ji = T⁺_ij, so one of each pair is solved and both rows are written. The diagonal ω_i = ω_j is a regular pair at Ω = 2ω_i.
- **Free-surface mesh with `qtf = 3`.** With `qtf_terms = full` and the sum on, the generated inner mesh and the mode count M follow the largest sum frequency, so the inner mesh has many more points than for the difference alone (more than ten times as many on a small hull). The difference files are then equal to those of `qtf = 1` only when `qtf_fs_panel` and `qtf_fs_modes` are both given; otherwise they differ a little. With `qtf_fs_panel` or `FSInnerMesh.pnl` the mesh is the one of the difference and the sum is less accurate. Give both keys, or use `qtf = 1` for the difference alone, to keep the cost down.
- **Memory.** The store adds the sum potential part, and with `qtf_direct = 1` T⁺_dir, of the same size as the difference stores; the logged size counts them.
- **Accuracy.** The direct sum agrees with the indirect sum of the same run to the mesh error of the hull, which falls as the hull mesh is refined; a coarse hull bottom makes the gap larger in the moments than in the forces.

Sources: Chen, X.-B. (2007). Middle-field formulation for the computation of wave-drift loads. J. Eng. Math. 59, 61-82. Newman, J. N. (1974). Second-order, slowly-varying forces on vessels in irregular waves. Proc. Int. Symp. Dynamics of Marine Vehicles and Structures in Waves, London, 182-186. Molin, B. (1979). Second-order diffraction loads upon three-dimensional bodies. Applied Ocean Research 1(4), 197-202. Ogilvie, T. F. (1983). Second-order hydrodynamic effects on ocean platforms. International Workshop on Ship and Platform Motions, Berkeley. Sharma, J. N., Dean, R. G. (1981). Second-order directional seas and associated wave forces. SPE Journal 21(1), 129-140. Kim, M.-H. and Yue, D. K. P. (1989). The complete second-order diffraction solution for an axisymmetric body. Part 1. Monochromatic incident waves. J. Fluid Mech. 200, 235-264. Cong, P., Teng, B., Chen, L., Gou, Y. (2020). A novel solution to the second-order wave radiation force on an oscillating truncated cylinder based on the application of control surfaces. Ocean Engineering 204, 107278. Pólya, G. and Szegő, G. (1951). Isoperimetric Inequalities in Mathematical Physics. Princeton University Press. WAMIT v7.5 user manual, ch. 15. Zhang, Y., Tosdevin, T., Tang, T., Stansby, P. (2026). Second-order difference-frequency effects on floating wind turbine responses: full-QTF modelling and experimental comparison. arXiv:2610.27940.

### Second-order motion

With `qtf_motion = 1` each free body moves under the second-order force of every row of the QTF files. Per second-order frequency Ω the run solves

[−Ω²(M + A(Ω)) − iΩ(B(Ω) + B_ext) + C + K_ext] ξ⁽²⁾ = F⁽²⁾,

with M = `MATX`, B_ext = `BLNR`, C = `CRS` and K_ext = `KSTF` of `Hydrostatic.in`, and F⁽²⁾ the force of the row. Ω is ω_i − ω_j for the difference and ω_i + ω_j for the sum. A(Ω) and B(Ω) are the added mass and radiation damping at |Ω|, computed from the radiation potentials ψ that the potential part already solves at every Ω, so the motion needs no extra solve. The force of `QTF_Indirect.11d` and `.11s` gives `QTF_Motion.11d` and `.11s`; with `qtf_direct = 1` the force of `QTF_Direct.12d` and `.12s` also gives `QTF_Motion.12d` and `.12s` ([output files](output_files.md#wamit-format)). A fixed body (`motion_state = fixed`) gets no motion and no motion file, and the log says so. The derivation and the time convention are in [Developer Documentation](developer_documentation.md#second-order-motion).

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
  Every other Ω gets one line ` QTF motion: Omega = ... rad/s, rcond ...`. To get the static diagonal, give K_ext in surge, sway and yaw; heave, roll and pitch have the hydrostatic C. A spring in surge alone leaves sway and yaw free, so Ω = 0 is still left out, surge included. An Ω > 0 with a small but regular matrix, such as a lightly damped slow-drift resonance, is written as computed.
- **Example.** A free truncated cylinder (radius 10 m, draft 10 m, depth 100 m) gets springs for a natural period of 100 s and dampers of 5% of critical in surge, sway and yaw, with ω_n = 2π/100 = 0.0628 rad/s: K_11 = K_22 = 1.96063e4 N/m and B_11 = B_22 = 3.12044e4 N s/m (from M + A_11), K_66 = I_66 ω_n² = 6.20126e5 N m and B_66 = 0.1 I_66 ω_n = 9.86960e5 N m s (A_66 = 0 for an axisymmetric hull). In `Hydrostatic.in` they are the (1,1), (2,2) and (6,6) entries of the external restoring and external linear damping matrices. The slow-drift resonance is then at Ω ≈ 0.063 rad/s; to put a pair on it, add a frequency 0.063 rad/s above one of the grid, for example 1.013 rad/s to a grid with 0.95 rad/s. With these springs every Ω, Ω = 0 included, is in the files.
- **Damping: give low-frequency damping as `BLNR`.** The motion takes B(Ω) + `BLNR` only. `BQDR` is a linearisation on the first-order motion amplitude and has no second-order form, so the motion does not use it. If any entry of `BQDR` is not zero, the log warns once:
  ```
   QTF warning: the quadratic damping BQDR is not used by the second-order motion (qtf_motion = 1),
   which has no second-order form of it; the motion takes B(Omega) + BLNR. Give the low-frequency
   damping of the mooring or the hull as BLNR if it is needed there.
  ```
  The damping sets the height of a slow-drift peak: about 1/(2ζ) times the static response, with ζ the damping ratio from `BLNR` (10 times for ζ = 5% in the example above). Without `BLNR` the peak is set by B(Ω) alone.
- **Sum frequency.** `QTF_Motion.11s` has the rows of `.11s`. With `qtf_terms = body` the sum warning adds ` The same holds for QTF_Motion.11s, the motion from that force.`, and with `If_remove_irr_freq = 0` the cutoff warning adds ` QTF_Motion.11s leaves out the same pairs.`
- **Several bodies.** M, C, K_ext and B_ext are block-diagonal, from the `Hydrostatic_i.in` file of each body; A(Ω) and B(Ω) couple the bodies.
- **The radiated wave in the field.** With `qtf_field = 1` and `qtf_direct = 1`, `QTF_Field.txt` adds the wave radiated by the motion of `QTF_Motion.12d` and `.12s`, φ_R = Σ_j (−iΩ ξ_j) ψ_j, to `Re_dir` and `Im_dir`, with p = iρΩ φ_R and η = (iΩ/g) φ_R as for φ_D. The column `has_rad` marks the rows that have it: 0 on rows without a direct part and on rows of an Ω left out of the motion. With `qtf_direct = 0` no radiated wave is added: adding a complete term to the partial field of `QTF_FieldPartial.txt` would mix the two.
- **Cost and memory.** Per Ω the run forms A and B from ψ and factors one 6N × 6N matrix (N bodies); this is small next to the radiation solves. The store adds ξ, 16 · 6N · (frequencies)² · (heading pairs) · (routes) bytes for the difference and as much again for the sum, with one route without `qtf_direct` and two with it; the log gives it (` QTF motion: of which ... MB for the second-order motions (qtf_motion = 1)`). With the radiated wave in the field, each Ω also takes one field pass of the 6N ψ columns over the points; its time is in the field time of the `QTF direct` log line, and ψ at the points (16 · (points) · 6N bytes per Ω) is not in the logged size.
- **Accuracy.** The motion error is mostly the force error: the part from A(Ω), B(Ω) and the polygon's M and C is small and falls with the mesh. Converge the force first.

### Partial second-order field at the field points

With `qtf_field = 1`, HAMS-MREL writes `QTF_FieldPartial.txt`: the second-order pressure at the field points below z = 0 and the second-order free-surface elevation at the points on z = 0, for every frequency and heading pair of `QTF_Quadratic.10d` (difference frequency) and `.10s` (sum frequency). The field points are those of `ControlFile.in` (`Number_of_field_points`). The kernels follow `qtf`, the heading pairs `mean_drift`, and `qtf_max_dw` limits the difference pairs, as for the QTFs. The columns, the conjugation and the scaling are in [Output Files](output_files.md#wamit-format).

**The file is partial.** It holds two parts:
- the quadratic part: products of the total first-order field at the point (incident, diffracted, and radiated by the motions of a free body);
- the second-order incident wave of each pair (Sharma and Dean 1981); for one wave with itself at the sum frequency this is the second-order Stokes wave (Dean and Dalrymple 1991).

It omits the second-order diffracted wave and, for a free body, the wave radiated by its second-order motion. Near the bodies the omitted part is not small:
- On the hull, the second-order diffracted potential cancels the normal velocity of the second-order incident potential (Molin 1979), so next to the hull it is of the size of the incident part.
- On the free surface it is forced by products of the first-order field that contain the scattered and radiated first-order waves (Molin 1979). Near a body that is not small against the wave length these waves are of the size of the incident wave, so this forcing is of the size of the one that makes the second-order incident wave. Kim and Yue (1989, 1990) solve this complete problem for an axisymmetric body.

So the omitted part comes from the bodies and falls off with distance from them, as their first-order waves do. It is small only where those waves are small against the incident wave, for example around a hull small against the wave length, where the file gives the second-order Stokes wave. Use the file for checks, and for the field away from the bodies with this in mind.

The equations are Bernoulli's, expanded about z = 0 (Dean and Dalrymple 1991):
- pressure at a fixed point below z = 0: p⁽²⁾ = −ρ[∂φ⁽²⁾/∂t + ½|∇φ⁽¹⁾|²], the dynamic second-order pressure only. A point between the wave troughs and z = 0 is not always wet, so keep pressure points below the troughs;
- elevation at a point on z = 0: η⁽²⁾ = −(1/g)[∂φ⁽²⁾/∂t + ½|∇φ⁽¹⁾|² + η⁽¹⁾ ∂²φ⁽¹⁾/∂z∂t],

with φ⁽²⁾ the second-order incident potential. The formulas per frequency pair are in [Developer Documentation](developer_documentation.md#partial-second-order-field).

**Points.** As for `PressureElevation.6p`, a point with |z| ≤ 1e-6 is an elevation point (its z is set to 0) and a point with z < −1e-6 a pressure point. The run stops with a message naming the point and the body when a point is:
- above the free surface (z > 1e-6);
- on a hull: within one panel size of the centroid of a hull panel, with the size of that panel;
- inside a body: the hull panels (symmetry images included) and their mirror image in z = 0 enclose the point, judged by the winding number (the solid angle they subtend at the point, over 4π, is more than 0.5). This also rejects a point on z = 0 inside the waterplane of a surface-piercing body. The test needs a watertight hull.

**Route and cost.** The first-order potential and its gradient at the points always come from the Green function, whatever `velocity_method` or `cs_velocity` is, with central differences of step δ = 0.02 times the size of the nearest hull panel. At each frequency of the QTF grid, each stencil point takes one influence row (the Green function and its normal derivative from every hull panel), which serves all headings and is not stored. A pressure point takes 7 rows: the point and x ± δ, y ± δ, z ± δ. An elevation point takes 5: the point and x ± δ, y ± δ, with ∂φ/∂z from the linear free-surface condition. A pressure point less than δ below z = 0 also takes 5 rows. The log gives the time at each frequency.

- The store keeps φ and ∇φ at each point, 64 bytes per point, frequency and heading. It is allocated at set-up, so a run that does not fit stops before the first-order solves, and the log gives its size. `qtf_field = 1` also keeps the first-order states of every heading of a frequency, as `mean_drift = 2` does.
- In finite depth above ν·h = 50 (ν = ω²/g), the log warns at each such frequency that the field values may scatter, for the reason given under [Control-surface settings](#control-surface-settings).
- `pre_postprocessing_python/SecondOrder/qtf_field_plots.py` maps the file and writes grids of field points: see [Python Scripts](python_scripts.md#partial-second-order-field-maps).

**Direct part (`qtf_field = 1` with `qtf_direct = 1`).** With both keys, HAMS-MREL also writes `QTF_Field.txt`. Its 15 first columns are those of `QTF_FieldPartial.txt`, which is written as without `qtf_direct`. Three columns follow: `Re_dir`, `Im_dir` and `has_direct` (columns 16 to 18), and with `qtf_motion = 1` a fourth, `has_rad` (column 19; see [Second-order motion](#second-order-motion)). The columns are in [Output Files](output_files.md#wamit-format).
- **What the direct part is.** The potential φ_D of the [direct method](#direct-method), solved at each difference frequency Ω = ω_i − ω_j (rows K = 1) or sum frequency Ω = ω_i + ω_j (rows K = 2), evaluated at the field points by Green's second identity over the hull: 4π φ_D(x) = ∫_S [G ∂φ_D/∂n − φ_D ∂G/∂n] dS, with the body forcing as the right side of ∂φ_D/∂n. It enters as the second-order incident potential does: p = iρΩ φ_D and η = (iΩ/g) φ_D, per unit amplitude pair.
- **It is not added to Re and Im.** The columns Re, Im, Re_quad, Im_quad, Re_inc, Im_inc are the partial value of `QTF_FieldPartial.txt`. The value with the direct part is Re + Re_dir and Im + Im_dir.
- **Rows with a direct part.** `has_direct` = 1 on the difference rows (K = 1) with ω_i ≠ ω_j for the pairs the direct solve covers (those within `qtf_max_dw`), and on the sum rows (K = 2) of the pairs of `QTF_Direct.12s` (the diagonal included, within `qtf_max_sw` and the irregular-frequency cutoff). The mean rows (ω_i = ω_j, K = 1) and the rows outside those limits have `has_direct` = 0 and `Re_dir` = `Im_dir` = 0.
- **What is not included.** The free-surface forcing of φ_D is not in the direct part. With `qtf_terms = full` it is in `QTF_Direct.12d`, but its part of φ_D is solved apart and not evaluated at the field points (Green's formula over the hull alone holds only for the body-forced φ_D). The wave radiated by the second-order motion of a free body is in it only with `qtf_motion = 1`, on the rows with `has_rad` = 1. The sum of the columns is therefore not the complete second-order field, and it does not satisfy the second-order free-surface condition on z = 0.
- **Points.** The point rules above apply. A point on z = 0 off the hull keeps the factor 4π: the image term of G coincides with 1/r there, so the singularity is 2/r over the half space below z = 0. Approaching the hull, the field tends to the value on the hull panel.
- **Cost and memory.** Per second-order frequency, each point takes one influence row of the Green function over the hull and one B2 integral for each block of pairs (normally one block). The store keeps φ_D at each point, frequency pair and heading pair, 16 bytes each: 16 · (points) · (frequencies)² · (heading pairs) bytes, in addition to the store of the partial field, and as much again for the sum rows with `qtf` 2 or 3. It is allocated at set-up and the run stops with a message if it does not fit; use fewer frequencies, headings or field points, or `qtf_direct = 0`.
- **Finite depth.** The direct part at the points uses the Green function without the Limes correction, as the first-order field at the same points does, so isolated outliers of the series do not land on single points.

Sources: Dean, R. G. and Dalrymple, R. A. (1991). Water Wave Mechanics for Engineers and Scientists. World Scientific. Kim, M.-H. and Yue, D. K. P. (1989). The complete second-order diffraction solution for an axisymmetric body. Part 1. Monochromatic incident waves. J. Fluid Mech. 200, 235-264. Kim, M.-H. and Yue, D. K. P. (1990). Part 2. Bichromatic incident waves and body motions. J. Fluid Mech. 211, 557-593. Molin, B. (1979). Second-order diffraction loads upon three-dimensional bodies. Applied Ocean Research 1, 197-202. Sharma, J. N. and Dean, R. G. (1981). Second-order directional seas and associated wave forces. Soc. Petroleum Eng. J. 21, 129-140.
