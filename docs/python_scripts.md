# Python Scripts for Pre and Post Processing

HAMS-MREL includes a set of Python scripts for generating input files, postprocessing output, and computing generalized mode shapes. The scripts are organized as follows:

```
HAMS_InputOutput/
├── Input/
│   ├── Single_body/
│   │   ├── HAMSInput.py
│   │   └── HAMS_single_bodies.py
│   ├── Multi-body/
│   │   ├── HAMSMultiInput.py
│   │   └── HAMS_multi_bodies.py
│   └── Mesh/
│       ├── hams_mesh.py
│       └── tests/test_hams_mesh.py
└── Output/
    ├── Single_body/
    │   ├── Example_single_body.py
    │   └── postprocessing_functions_hams_single_bodies.py
    └── Multi-body/
        ├── Multibodies_HAMS.py
        ├── Obtaining_RAOs_from_HAMS_PA_array.py
        └── post_processing_HAMS_multibodies.py
HAMS_to_WECSim/
├── HAMS_to_WECSim.py
└── postprocessing_functions_Mat.py
generalized_modes.py
SecondOrder/
├── make_mesh.py
├── make_drift_fixtures.py
├── drift_convergence.py
├── truncated_cylinder_reference.py
├── linton_evans_reference.py
├── qtf_field_plots.py
└── qtf_compare.py
```

Dependencies: `numpy`, `matplotlib`, `shapely`. The WEC-Sim conversion scripts additionally use `cmath`, `os`, and `datetime` (all standard library).

Where a script takes the gravitational acceleration `g`, its default is 9.80665 m/s², the standard value that the Fortran code uses everywhere (`G` in `Const_mod`, [src/WavDynMods.f90](../src/WavDynMods.f90)). Keep the two the same.

## Input Generation

### Single Body

**`HAMS_single_bodies.py`** contains two classes and two wrapper functions for generating the `ControlFile.in` and `Hydrostatic.in` files for a single-body simulation.

- `HAMSControlInput` — class that constructs the control file line by line. Methods write each section: water depth, wave frequencies, wave headings, solver options, pressure/elevation field points, and closing notes.
- `HAMSHydrostaticFile` — class that constructs the hydrostatic file. Writes the center of gravity, 6×6 mass matrix, external linear damping matrix, external quadratic damping matrix, hydrostatic restoring matrix, and external restoring matrix.
- `create_control_file_hams(...)` — wrapper function that instantiates `HAMSControlInput` and calls all write methods in the correct order.
- `create_hydrostatic_file(...)` — wrapper function that instantiates `HAMSHydrostaticFile` and writes the hydrostatic data.

**`HAMSInput.py`** is an example script that demonstrates how to use `HAMS_single_bodies.py`. The user sets simulation parameters as variables at the top of the script (water depth, frequency range, headings, rotation center, solver options, mass and stiffness matrices, etc.), then calls `create_control_file_hams` and `create_hydrostatic_file`. The script also creates the `Input/` and `Output/` directory structure and moves the generated files into `Input/`.

### Multi-Body

**`HAMS_multi_bodies.py`** contains the `HAMSInput` class and a `create_control_file_hams` wrapper function for generating the `ControlFile.in` for multi-body simulations. The class extends the single-body version with:

- `add_body_info()` — writes the number of bodies and the LCS (Local Coordinate System) coordinates and orientation for each body.
- Per-body rotation centers (`Reference_body_center_i`) in the solver options section.

**`HAMSMultiInput.py`** is an example script that demonstrates how to use `HAMS_multi_bodies.py`. It sets up a 10-body simulation with LCS coordinates and rotation centers per body. Note that the `Hydrostatic[_i].in` files must be created separately for multi-body cases. `Mesh/hams_mesh.py` (below) writes them together with the meshes.

### Meshes

**`Mesh/hams_mesh.py`** builds simple geometric meshes and writes the mesh and hydrostatic input files for one body, an array of identical bodies, or a body with dipole (zero-thickness) sheets. It needs `numpy`, and `matplotlib` for `--plot`. Each shape is a function that returns a `Body` and a subcommand:

| Subcommand | Function | Shape |
| :--- | :--- | :--- |
| `cylinder` | `cylinder(radius, draft, ...)` | vertical surface-piercing cylinder, wall and flat bottom |
| `box` | `box(length, width, draft, ...)` | rectangular surface-piercing box, length along x |
| `hemisphere` | `hemisphere(radius, ...)` | lower half of a sphere centred on the waterline |
| `sphere` | `sphere(radius, centre_depth, ...)` | sphere below the free surface |
| `spar` | `spar(profile, ...)` | body of revolution of a polyline profile `r z` from the waterline down, the bottom disc added: a stepped or tapered spar |
| `cone` | `cone(top_radius, bottom_radius, draft, ...)` | truncated cone |
| `heave-plate` | `heave_plate(column_radius, draft, plate_radius, ...)` | cylinder column with a dipole annulus at its bottom, sharing the nodes of the column's bottom ring |
| `--disc R DEPTH [X Y]` on any subcommand | `disc(radius, depth, ..., centre=(x, y))` | circular dipole sheet, added to the hull |
| `--plate L W DEPTH [X Y]` on any subcommand | `plate(length, width, depth, ..., centre=(x, y))` | rectangular dipole sheet, added to the hull |

A sheet is always added to a hull body, with `add_sheets(body, sheet, ...)` or the repeatable `--disc` and `--plate` options. It is never written alone, because HAMS stops at input when `HullMesh.pnl` has no hull panels (`NELEM = 0`). `write_input` refuses a sheet without a hull and says so. The sheets take the panel size, and `--edge-grading` grades them towards their rims. They add no volume, so `Hydrostatic.in` comes from the hull, or from `--mass`, `--cog` and `--gyration`.

```
python hams_mesh.py cylinder --radius 5 --draft 10 --panel-size 0.8 --waterline-grading 0.5 --edge-grading 0.5 --lid
python hams_mesh.py spar --profile 3.25 0 3.25 -4 4.7 -12 4.7 -120 --panel-size 1.5 --edge-grading 0.5
python hams_mesh.py heave-plate --column-radius 10 --draft 9.9 --plate-radius 15 --panel-size 0.75 --edge-grading 1
python hams_mesh.py cylinder --radius 3 --draft 6 --panel-size 0.5 --disc 2 12 --plate 4 2 10 8 0 --edge-grading 1
python hams_mesh.py cylinder --radius 2 --draft 4 --panel-size 0.4 --at 0 0 --at 20 0 --at 10 15 30 --lid --out-dir Input
```

- **Fineness.** `--panel-size` sets a target panel edge, from which the panel counts follow (around a circle, rounded up to a multiple of 4). The counts of each shape (`--ntheta`, `--nz`, `--nr`, `--nx`, `--ny`, `--nphi`, `--nplate`, and `--counts` per profile segment for `spar`) override it.
- **Grading.** `--waterline-grading G` makes the rows finer towards z = 0, and `--edge-grading G` towards sharp edges (the bilge corner of a cylinder or box, the corners of a spar profile) and towards the free rim of a dipole sheet. G runs from 0 (uniform) to 1 (sine spacing). A segment fine at one end has its nodes at (1 − G) u + G f(u), with f(u) = 1 − cos(πu/2) or sin(πu/2); fine at both ends, f(u) = (1 − cos πu)/2 with the larger G. The counts stay those of uniform spacing, so grading makes the end panels smaller than the panel size. Grade a dipole annulus or disc towards its rim, with G = 1 as recommended in [Input Files](input_files.md#dipole-panels-optional-single-body).
- **Outputs**, in `--out-dir` (default the current folder): `HullMesh.pnl`; `WaterPlaneMesh.pnl` with `--lid` (for `IRSP = 1`); and `Hydrostatic.in`. The script prints the volume, waterplane area, mass and restoring it wrote, and the `ControlFile.in` lines to paste: the rotation centre and `If_remove_irr_freq`.
- **Dipole panels** are written in `HullMesh.pnl`, with one node list shared with the hull. Line 4 is `NELEM NTND ISX ISY NDIPOLE`, and the fifth number is written only when `NDIPOLE > 0`. The hull panels are followed by the line `# Start Definition of Dipole Panels` and the `NDIPOLE` dipole panels, in the same panel format and numbered from 1. No `DipoleMesh.pnl` is written. Dipole sheets support one body only, `ISOLV = 1`, no symmetry and first order only, and the script reminds you of this.
- **Arrays.** With more than one `--at X Y [ANGLE]` (angle about z in degrees), the same body is written at each position as `HullMesh_i.pnl`, `WaterPlaneMesh_i.pnl` and `Hydrostatic_i.in`, and the `Number_of_bodies`, `LCS_i` and `Reference_body_center_i` lines are printed.
  - The meshes are in the body's own frame, placed by the LCS lines, as `ReadBodyMeshMulti` does: rotation about z, then the shift.
  - The centre of gravity, rotation centre and matrices of `Hydrostatic_i.in` are global, as in `examples/example1_multi_bodies`.
  - The LCS and rotation-centre lines use the fixed formats of `ReadControlFile` (`26x,4f12.3` and `28x,3f12.3`). The LCS angle is in radians with 3 decimals, so positions and angles are rounded to what HAMS reads before the hydrostatics are computed.
- **Hydrostatics.** The volume, the centre of buoyancy and the waterplane moments come from the panels themselves, by the divergence theorem on the hull closed at z = 0. This is exact for the flat panels, with a 7-point degree-5 rule on each triangle. By default the mass is a uniform solid filling the submerged volume: mass ρV, centre of gravity at the centre of buoyancy, and its inertia from the same panels. `--mass`, `--cog` (in the body frame) and `--gyration` (radii about the centre of gravity) replace it, with a warning if the mass is not ρV. The restoring is that of WAMIT, about the rotation centre (x, y, `--centre-z`), C33, C34, C35, C44, C45, C46, C55 and C56. External damping and stiffness are zero. ρ = 1000 kg/m³ and g = 9.80665 m/s², as `Const_mod`. For the uniform cylinder the matrices equal those of `SecondOrder/make_mesh.py` to 1e-9.
- **Checks before writing**:
  - normals as HAMS computes them (`CalTransNormals`) point out of the fluid;
  - the hull is closed apart from its waterline on z = 0, and every edge is shared by two panels in opposite order;
  - no node is above z = 0, and no hull panel centroid is on z = 0;
  - the dipole input checks of `DIPOLE_SETUP` ([ReadPanelMesh.f90](../src/ReadPanelMesh.f90)): every dipole node at least one panel size below z = 0, no dipole centroid on a hull centroid, and no overlap with a coplanar hull panel;
  - every lid node lies on z = 0, with a warning if a lid centroid is closer to the hull wall than 0.3 waterline rows ([Lid mesh requirements](input_files.md#lid-mesh-requirements)). The lid's outer ring is one waterline row wide for that reason.
- **`--plot PNG`** draws the mesh: hull grey, dipole sheets orange, lid blue.
- **Tests.** `python -m pytest pre_postprocessing_python/HAMS_InputOutput/Input/Mesh` runs 17 tests (about 10 s), one per shape. Two of them cover a disc and a plate added to a hull (sheet area exact, hydrostatics equal to the hull's) and the refusal of a sheet alone. They check the volume and waterplane area against the exact values within the panel error: equal to the inscribed polygon or frustum where that is exact, and within the sagitta bound for the spheres. They also read the files back as HAMS reads them, and test the checks, the lid, the hydrostatics about an offset centre and with a given mass, the multi-body set and its control lines.


## Output Postprocessing

### Single Body

**`postprocessing_functions_hams_single_bodies.py`** provides two functions for reading HAMS output files in Hams_format:

- `get_mass_and_damping_matrices_hams(added_mass_file_location, added_damping_file_location, dof)` — reads `OAMASS{i}.txt` and `ODAMPING{i}.txt` and returns the frequency list, added mass coefficients, and radiation damping coefficients for a given DOF.
- `get_external_force_hams(external_force_file_location, column)` — reads `OEXFOR{i}.txt` and returns the excitation force values for a given column (real or imaginary part).

**`Example_single_body.py`** is an example script that demonstrates postprocessing for a single point absorber in heave. It reads the HAMS output, plots added mass, radiation damping, and excitation force coefficients, solves the frequency-domain equation of motion with and without a PTO system, computes the natural frequency, and plots the RAOs and average absorbed power. The dynamic equation solved is:

```
[-ω²(m + a₃₃) + iω(b₃₃ + B_PTO) + c₃₃] · η₃₃ = Fₑ
```

### Multi-Body

**`post_processing_HAMS_multibodies.py`** provides functions for reading and reorganizing HAMS output for multi-body simulations:

- `get_mass_and_damping_matrices_hams_multibodies(...)` — reads added mass and radiation damping from HAMS output files for a given DOF. Same interface as the single-body version.
- `get_external_force_hams_multibodies(...)` — reads excitation forces from HAMS output files.
- `compile_mass_damping_matrices_multiple_bodies(n_bodies, length_lists, added_mass_final, radiation_damping_final)` — reorganizes the raw coefficient lists into a compiled format indexed by body pair. Returns lists of size `n_bodies²`, each containing the coefficient values across all frequencies.
- `compile_excitation_force_matrix_multiple_bodies(n_bodies, length_lists, exciting_force_final, n_directions)` — reorganizes excitation force lists indexed by body and wave direction. Returns a list of size `n_bodies × n_directions`.
- `hydrodynamic_coefficients_exciting_forces_all_bodies(...)` — main function that reads all HAMS output files for a DOF, optionally strips zero/infinite frequency entries, and returns the compiled added mass, radiation damping, and excitation force matrices along with the number of bodies and wave directions.
- `hydrodynamic_coefficients_exciting_forces_all_bodies_zero_inf(...)` — same as above but returns only the zero and infinite frequency coefficients.
- `obtain_specific_coefficient_matrix(added_mass_compiled, radiation_damping_compiled, n_bodies, x, y)` — retrieves added mass and radiation damping for a specific body pair (`x`, `y`).
- `obtain_specific_exciting_force_matrix(exciting_force_compiled, x, j, n_directions)` — retrieves excitation force for a specific body `x` and direction `j`.

**`Multibodies_HAMS.py`** contains the `Multibodies_HAMS` class that provides higher-level postprocessing for multi-body simulations. Key methods:

- `get_number_of_bodies()` — reads the number of bodies from `ControlFile.in`.
- `get_frequency_information()` — parses all simulation parameters from `ControlFile.in`, including frequencies, headings, body positions, water depth, and number of panels per body. Converts between frequency types (wave number, frequency, period) including a shallow-water dispersion solver.
- `hams_coefficient_conversion(body_number_1, body_number_2, dof_1, dof_2)` — converts body/DOF pairs to global DOF indices in the combined `6×NBODY` system.
- `create_hams_compiled_results(hams_output_location, run)` — reads all HAMS output files, compiles hydrodynamic coefficients and excitation forces across all DOF combinations, and writes them to individual text files in `Results_full_HAMS/`.
- `create_hams_compiled_results_zero_inf(hams_output_location, run)` — same as above but for zero and infinite frequency only.
- `compile_hams_results_per_frequency(run)` — reads the compiled results from `Results_full_HAMS/` and reorganizes them into per-frequency output files in `Results_frequency_HAMS/`, with one file per frequency containing all DOF combinations.
- `compile_hams_results_per_frequency_zero_inf(run)` — same as above for zero/infinite frequencies.

**`Obtaining_RAOs_from_HAMS_PA_array.py`** is an example script that demonstrates RAO computation for an array of point absorbers. It uses the `Multibodies_HAMS` class to compile the per-frequency coefficient and excitation matrices, then solves the coupled equation of motion:

```
[-ω²(M + A) + iωB + K] · ξ = Fₑ
```

where `M`, `A`, `B`, `K` are `6N × 6N` block-diagonal matrices (N = number of bodies) and `Fₑ` is the `6N × 1` excitation force vector. The resulting RAOs are written to a text file in WAMIT format.


### Multi-Body (Matrix-Based Variant)

**`postprocessing_functions_Mat.py`** provides an alternative set of multi-body postprocessing functions with a simplified interface tailored for matrix-based workflows (e.g. WEC-Sim conversion). The functions mirror those in `post_processing_HAMS_multibodies.py` but omit wave-direction indexing, assuming a single wave heading. Key functions:

- `get_mass_and_damping_matrices_hams_multibodies(added_mass_file_location, added_damping_file_location, dof)` — reads HAMS output files and returns the frequency list, added mass coefficients, and radiation damping coefficients for a given DOF.
- `get_external_force_hams_multibodies(external_force_file_location, column)` — reads excitation force values for a given column (real or imaginary part).
- `compile_mass_damping_matrices_multiple_bodies(n_bodies, length_lists, added_mass_final, radiation_damping_final)` — reorganizes raw coefficient lists into compiled format indexed by body pair. Returns lists of size `n_bodies²`.
- `compile_excitation_force_matrix_multiple_bodies(n_bodies, length_lists, exciting_force_final)` — reorganizes excitation force lists indexed by body. Returns a list of size `n_bodies`.
- `hydrodynamic_coefficients_exciting_forces_all_bodies(...)` — main function that reads all HAMS output files for a DOF, optionally strips zero/infinite frequency entries, constructs complex excitation forces from real and imaginary parts, and returns the compiled matrices.
- `hydrodynamic_coefficients_exciting_forces_all_bodies_zero_inf(...)` — same as above but returns only the zero and infinite frequency coefficients.
- `obtain_specific_coefficient_matrix(added_mass_compiled, radiation_damping_compiled, n_bodies, x, y)` — retrieves added mass and radiation damping for a specific body pair (`x`, `y`).
- `obtain_specific_exciting_force_matrix(exciting_force_compiled, x)` — retrieves excitation force for a specific body `x`.


## Generalized Modes

**`generalized_modes.py`** computes generalized mode shape functions for use with the generalized modes feature in HAMS-MREL. The script produces `gen_mod_{body}_{dof}.txt` files that can be placed in the input directory to override standard rigid-body normals. See the [Input Files](input_files.md#gen_mod-optional) documentation for details on how these files are read by the Fortran code.

The script contains the following functions:

- `read_pnl_mesh(filepath)` — reads a HAMS `.pnl` hull mesh file and returns node coordinates and panel connectivity.
- `extract_body_number(mesh_filepath)` — extracts the body number from the mesh filename (e.g., `HullMesh_2.pnl` → 2).
- `read_lcs_from_controlfile(control_filepath, body_number)` — reads the LCS translation offsets and z-axis rotation angle for a given body from `ControlFile.in`.
- `local_to_global(coords_local, x0, y0, z0, rot_deg)` — transforms coordinates from the local body frame to the global frame by applying z-axis rotation followed by translation.
- `compute_centroids(nodes, panels)` — computes the centroid of each panel as the arithmetic mean of its vertex coordinates.
- `compute_generalized_normal(centroids_global, nx, ny, nz)` — evaluates the generalized mode shape function at each panel centroid as `n_j = nx * X + ny * Y + nz * Z`.

The user sets the mesh file, control file, DOF, and shape function components (`nx`, `ny`, `nz`) as variables in `main()`, then runs the script to produce the output file.


## HAMS to WEC-Sim Conversion

**`HAMS_to_WECSim.py`** converts HAMS-MREL multi-body results into a WAMIT-format `.out` file that can be used as input for WEC-Sim. The script imports postprocessing functions from `postprocessing_functions_Mat.py` and contains the `HAMS_to_WAMIT` class with the following methods:

- `__init__(filename, hams_result_folder, hams_input_folder, rho)` — initializes the converter with output filename, HAMS directory paths, and water density (default 1000 kg/m³).
- `create_initial_wamit_output()` — writes the WAMIT-format header block to the output file, including MREL attribution.
- `get_number_of_bodies()` — reads the number of bodies from `ControlFile.in`.
- `get_frequency_information(g)` — parses all simulation parameters from `ControlFile.in` including frequencies, headings, body positions (LCS), rotation centers, water depth, and panel counts. Converts between frequency types (wave number, frequency, period) for all four HAMS input frequency types. Stores period, wavenumber, and angular frequency lists for WAMIT output formatting.
- `obtain_wavenumber_shallow_water(period, depth, g)` — iterative shallow-water dispersion solver that computes the wavenumber from wave period and water depth.
- `hams_to_wamit_coefficient_conversion(body_number_1, body_number_2, dof_1, dof_2)` — maps body/DOF pairs to global indices in the combined `6N` system used by WAMIT.
- `hams_to_wamit_coefficient_conversion_excitation(body_number_1, dof_1)` — maps body/DOF to a global excitation force index.
- `create_hams_compiled_results_for_wamit(hams_output_location, run)` — reads all HAMS output files, compiles hydrodynamic coefficients (added mass, radiation damping) and excitation forces across all DOF combinations, and writes them to individual text files in `Results_full_HAMS/`. Excitation forces are stored as magnitude and phase.
- `create_hams_compiled_results_for_wamit_zero_inf(hams_output_location, run)` — same as above but for zero and infinite frequency added mass only, writing to `Results_full_HAMS_inf_zero/`.
- `compile_hams_results_per_frequency_for_wamit(run)` — reads the compiled results from `Results_full_HAMS/` and reorganizes them into per-frequency files in `Results_frequency_HAMS/`, with one coefficient file and one excitation file per frequency.
- `compile_hams_results_per_frequency_for_wamit_zero_inf(run)` — same as above for zero/infinite frequencies, writing to `Results_frequency_HAMS_inf_zero/`.
- `add_input_file_description()` — writes the WAMIT input file description block listing hull mesh files.
- `add_frequency_information()` — writes the frequency table, gravity, water depth, solver parameters, and symmetry information in WAMIT format.
- `add_hydrostatic_data(volumes, hydrostatic, center_of_gravity, center_of_bouyancy, g)` — writes per-body hydrostatic data including displaced volumes, center of buoyancy, hydrostatic restoring coefficients (non-dimensionalized by `ρg`), center of gravity, and radii of gyration.
- `add_added_mass_zero_inf()` — writes the zero and infinite frequency added mass coefficient blocks, non-dimensionalized by `ρ`.
- `add_hydrodynamic_coefficients_per_frequency(omega_number)` — writes the added mass (non-dimensionalized by `ρ`) and radiation damping (non-dimensionalized by `ρω`) for a given frequency.
- `add_excitation_forces_per_frequency(omega_number, g)` — writes the diffraction exciting force magnitude (non-dimensionalized by `ρg`) and phase for a given frequency and wave heading.
- `add_hydrodynamic_coefficients_excitation_forces()` — loops over all frequencies and writes the combined coefficient and excitation force blocks.
- `add_initial_output_info()` — writes the WAMIT FORCE run header.

The script includes a two-stage execution workflow controlled by flags at the bottom:

1. **`create_data_from_HAMS`** — reads raw HAMS output, compiles all DOF combinations, and produces intermediate per-frequency result files.
2. **`create_final_WECSIM_files`** — reads the intermediate files and assembles the final WAMIT `.out` file with user-specified hydrostatic data, volumes, and centers of gravity/buoyancy.

## Mesh Generation for Second-Order Cases

**`SecondOrder/make_mesh.py`** writes a `HullMesh.pnl` (format in [Input Files](input_files.md#hullmeshpnl)) and a matching `Hydrostatic.in` for a surface-piercing vertical cylinder or box. It generated the meshes of the `test/data/drift_*` integration cases.

```
python make_mesh.py cylinder --radius 10 --draft 10 --ntheta 32 --nz 6 --nr 6 --out HullMesh.pnl
python make_mesh.py box --lx 20 --ly 10 --draft 5 --nx 10 --ny 5 --nz 4 --out HullMesh.pnl
```

- The hull has no lid. The script stops if any panel normal, computed as HAMS does, points into the fluid.
- `--x0` and `--y0` shift the body in the horizontal plane. `--half-y` keeps the panels with y > 0 and sets the y symmetry flag. `--graded` spaces the cylinder rows by a cosine, finer at the waterline and the bilge. `--hydro-out` sets the path of the hydrostatic file, which by default is `Hydrostatic.in` next to the mesh.
- The hydrostatics are those of a freely floating, uniform solid prism on the discretized waterplane polygon: mass $\rho V$, centre of gravity and of buoyancy at half draft, about the point $(x_0, y_0, 0)$. Give that point as the rotation centre in `ControlFile.in`. External damping and stiffness are zero. These are whole-body values, also with `--half-y`. The script stops if a metacentric height is not positive, and prints the values it writes.
- `python make_mesh.py oc4 --multiple 1|2 [--out-dir D]` writes the columns of the OC4 semisubmersible (1224 or 4608 panels), its `Hydrostatic.in` (platform mass and linearized mooring), and the control surface and free surface as `CSMesh_1.pnl` and `FSMesh_1.pnl` (for `cs_shape = file`), after checking them, and prints the checks. The sources of the model are given in the script.

**`SecondOrder/make_drift_fixtures.py`** regenerates all the `test/data/drift_*` inputs, including `ControlFile.in` and `SecondOrder.in`. See `test/data/drift_cyl/README.md`.

**`SecondOrder/drift_convergence.py`** generates and runs sets of mean drift and QTF cases of a cylinder, a box, the OC4 semisubmersible or an array of bottom-mounted cylinders over a series of meshes (HAMS-MREL, and Nemoh as a separate program), writes Slurm scripts for them, and collects the results into a table. Its `collect` step needs `scipy`, for the MacCamy-Fuchs value.

```
python drift_convergence.py list [--tier full|small] [--markdown] [--only PREFIX ...]
python drift_convergence.py generate ROOT [--case NAME | --group G]
python drift_convergence.py run ROOT --case NAME --hams HAMS_EXE [--threads N] [--dry-run]
python drift_convergence.py slurm ROOT --repo CLUSTER_CHECKOUT [--account A] [--only PREFIX ...] [--batch NAME] [--limit H:MM:SS] [--cpus N] [--memory GB]
python drift_convergence.py collect ROOT [--out FILE]
python drift_convergence.py check
```

- Every command takes the set of cases as a flag: none (mean drift), `--control` (the control-surface method), `--qtf`, `--qtf-potential` (the QTF potential part, the sum frequency and the second-order motion) or `--linton-evans` (the four-cylinder array). `list` prints the cases of a set with their meshes and estimated time and memory.
- `run --dry-run` runs nothing: it writes the case into a temporary folder and prints `ControlFile.in` and `SecondOrder.in` (HAMS), or `Nemoh.cal` and the head of `FSmesh.dat` (Nemoh), then deletes the folder.
- `slurm` writes one array script per size class. `--only` keeps the cases whose name starts with a prefix, `--batch NAME` puts all the chosen HAMS cases in one array (`NAME.txt`, `NAME.sbatch`), `--limit` sets the time limit of every array, and `--cpus` (default 32) and `--memory` (default 64 GB) set the resources of each run. The estimate of a `qtf_terms = full` case is that of `body` (the free-surface passes are not modelled), so its time limit is 4 times the estimate, at least 2 h.
- The external damping and stiffness written in `Hydrostatic.in` of some cases (`SHARP_DAMPING`, `SHARP_MOTION_K11` and the other constants) are derived beside the constants in the script; on a hull without them a case that needs them raises an error.

**`SecondOrder/truncated_cylinder_reference.py`** computes a semi-analytical reference (eigenfunction matching, no panels) for a sharp truncated cylinder: A, B, X, RAOs, far-field F_x, and near-field F_x, F_z and M_y with each of the six terms, fixed and free. It needs `numpy` and `scipy`. The method and its checks are in the script's docstring. Run it from its folder:

```
python -I truncated_cylinder_reference.py check [NAME ...]
python -I truncated_cylinder_reference.py reference --out tc_reference.csv
python -I truncated_cylinder_reference.py compare RESULTS --out tc_compare.csv [--reference tc_reference.csv]
```

- `check` runs the named self-checks, or all but the three slow ones (`uncertainty`, `compare_signs`, `potential_convergence`), and exits with 1 if one fails. A check that cannot run for want of data (`compare_signs` without HAMS results) prints SKIP and is not counted as passed. About 10 min on a laptop.
- `reference` writes the reference values at 11 frequencies, fixed and free, with their stated uncertainties, prints the convergence table and runs all the checks. The stated uncertainties hold for the default hull at the default resolution only; for any other cylinder or resolution they are `nan`, an empty field in the CSV.
- `qtf-potential --out CSV [--motion fixed|free|free_damped|both] [--light]` writes the quadratic part and the five body-forcing terms of the difference QTF on the grid of the QTF runs. `free_damped` is the free body with `drift_convergence.SHARP_DAMPING`, for the `q_sharp_free_damp_*` cases. `check damped_motions` checks that this damping enters the motions as −iωB (under 1 s).
- `qtf-potential ... --kernel sum` (default `diff`, whose output is unchanged) writes the sum-frequency QTF T+ at Omega = omega_i + omega_j: rows of kernel `sum` for every pair omega_i <= omega_j, the diagonal included, with the same columns, so the quadratic part plus the five terms (plus T_FS with `--free-surface`) is the reference of `QTF_Direct.12s`. The waves of the reference are 1/(-i) times unit waves, so every second-order value is multiplied by -1 to match the unit-wave convention of the HAMS files. The free-surface term uses the sum kernel with the second state unconjugated and three phase classes (k1 + k2 + K, -k1 + k2 + K, k1 - k2 + K), all positive, so the rotated contour needs no new rule. The new checks are `incident2_sum` (one wave with itself gives the second-order Stokes potential, 1e-10), `potential_sum_symmetry` (T+_ij = T+_ji, 1e-12), `fs_sum_bvp` (the direct solve of phi_D at the sum frequency equals ID + FS of the indirect terms, 1e-5, at N = 4000), `fs_sum_forcing_vs_incident` (the sum forcing of incident-only fields equals the sum-frequency incident wave's own, 1e-10) and `fs_sum_radius_invariance` (T_FS+ at R = 15 and 20 m, 1e-6).
- `qtf-potential ... --free-surface [--fs-radius R] [--fs-check]` adds the free-surface integral T_FS of the second-order potential, the columns `Fx_FS_re` … `My_FS_im`, so that the quadratic part plus all six terms is the complete difference QTF (the reference of `qtf_terms = full`). The integral over z = 0 is split at R (default 20 m): a polar grid inside, angular modes outside, and beyond about 2 km a rotated contour that sums the oscillating tail exactly. `--fs-check` also splits at 1.5 R and writes the relative change in `fs_check`. The sign of the term comes from a direct solve of the second-order problem (`check fs_sign_from_bvp`). The method and its accuracy are in the script's docstring, section "QTF free-surface integral". On a laptop at the default resolution it adds about 25 min per motion on the `hams` grid (50 min with `--fs-check`); run the full grid on a cluster.
- `qtf-motion --out CSV [--kernel diff|sum] [--k-ext K11] [--b-ext B11] [--grid motion|nemoh|hams] [--light] [--free-surface] [--from CSV]` writes the second-order motion xi_1, xi_3, xi_5 of the free cylinder (Phase 6c, the reference of `QTF_Motion.*`): for each pair of the kernel it takes T (the quadratic part, the five body-forcing terms and, with `--free-surface`, T_FS, computed as `qtf-potential` does or read from a `qtf-potential` CSV with `--from`) and solves [-Ω² (M + A) - iΩ (B + B_ext) + C + K_ext] ξ = ρ g T with A and B at |Ω| from the radiation potentials of T, M and C of the exact cylinder, and the surge spring and damper `--k-ext`, `--b-ext` (default the 100 s spring and 5 % damper of the motion cases, `SHARP_MOTION_K11` and `_B11` of `drift_convergence.py`, [second-order motion](input_files.md#second-order-motion); sway and yaw are not solved, so their springs are not needed), on the grid of the motion cases by default (`--grid motion`, `drift_convergence.QTF_MOTION_FREQUENCIES`). Ω = ω_i − ω_j for the difference (the matrix conjugated for Ω < 0, as in HAMS) and ω_i + ω_j for the sum; at Ω = 0 the equation is static. A pair whose matrix has rcond < 1e-12 gets empty `xi` fields, as HAMS leaves it out. The CSV has the columns `big_omega`, `rcond`, the total T (`Fx_total_re` …) and `xi1_re` … `xi5_im`. The cost is that of `qtf-potential` on the same grid; with `--from` it takes seconds.
- `check motion_ref_one_dof_limit motion_ref_ab_match` checks the motion solve (under 1 s; part of a bare `check`): heave alone against the one-degree-of-freedom closed form ξ = F/(c − Ω²(m + a) − iΩb) of the unit test `motion_one_dof`, at Ω = 2π/100, ±0.30, 1.40 and 0 rad/s, with surge and pitch zero (tolerance 1e-12), and A, B at Ω against the reference's first-order A, B at ω = Ω (1e-12).
- `check fs_forcing_vs_incident fs_radius_invariance fs_outer fs_sign_from_bvp` runs the four checks of the free-surface integral (about 1 min); they are part of a bare `check`.
- `compare` reads the `a_sharp_*` case folders of `drift_convergence.py` under `RESULTS` (`output/Wamit_format/`) and prints, per case and frequency, the percentage difference of HAMS from the reference. Folders that are not such cases, or have no output, are skipped with a note. 93 s, nearly all of it the reference rows; `--reference` reads them from a CSV of `reference` instead (under 1 s). The uncertainty shown after each drift percentage is the one the reference rows carry.

**`SecondOrder/linton_evans_reference.py`** computes the exact solution for wave scattering by an array of bottom-mounted vertical circular cylinders (Linton & Evans 1990), without panels: the exciting forces and moments per cylinder, the elevation at any point, the near-field mean drift per cylinder, the far-field mean drift of the array, and the mean drift on a vertical control circle. It needs `numpy` and `scipy` (and `matplotlib` for `plot`). The method and its checks are in the script's docstring. Run it from its folder:

```
python -I linton_evans_reference.py check [NAME ...]
python -I linton_evans_reference.py reference --out le_reference.csv
python -I linton_evans_reference.py compare RESULTS [--out le_compare.csv] [--reference le_reference.csv]
python -I linton_evans_reference.py plot [RESULTS] --out FOLDER [--reference le_reference.csv]
```

- `check` runs the 13 self-checks, or the named ones, and exits with 1 if one fails (about 4 min on a laptop).
- `reference` writes the rows of the HAMS case (`drift_convergence.py --linton-evans`: four cylinders of radius 5 m, depth 10 m, centres (±10, ±10) m, headings 0 and 45 deg, 17 frequencies, 16 field points) in the conventions of the HAMS files: kind `body` (X_1, X_2, X_4, X_5 as `ExcForce.3` writes them, and the near-field drift per body, moments about the body centre), `array` (far-field F_x, F_y and M_z about the origin) and `field` (the total elevation of `PressureElevation.6p`).
- `compare` reads the `l_square_*` case folders under `RESULTS` (`output/Wamit_format/`) and prints the largest error of each quantity per case: X forces and moments (complex, over the largest of each kind for that body), near-field F and M (over the largest load of the array at that frequency and heading), far-field F, far-field M_z and the elevation. `--out` writes the rows per frequency, heading and body.
- `plot` writes PNGs of the reference (the wave field, the loads per cylinder against one cylinder alone, Kim (1993), Fig. 3, and Linton & Evans (1990), Figs. 2 and 3 redrawn with the digitised values of `linton_evans_digitised.csv`, which is not in the repository: set `LINTON_EVANS_DIGITISED` to its path, or these two figures are skipped); with `RESULTS`, also HAMS in those figures and the elevation error at the field points (`plot_field_error`), as maps when the `l_square_field_*` cases are there; with `RESULTS`, also |X_1|, |X_2| and the near-field and far-field drift of HAMS against the reference, per heading, and the largest errors against the panel count.

**`SecondOrder/qtf_compare.py`** compares the QTFs of the `q_*` cases of `drift_convergence.py` with the exact reference of `truncated_cylinder_reference.py` and with Nemoh. The measure is |value − reference| / max |reference| per pair, and the table gives its largest and median over the pairs, per case. Run it from its folder:

```
python qtf_compare.py check
python qtf_compare.py RESULTS --reference qtf_reference.csv [--nemoh FOLDER] [--max-dw 0.3]
python qtf_compare.py RESULTS --potential qtf_potential.csv [--nemoh FOLDER] [--max-dw 0.3] [--exclude LO HI]
python qtf_compare.py RESULTS --direct [qtf_potential.csv] [--max-dw 0.3] [--exclude LO HI] [--by-omega]
python qtf_compare.py RESULTS --direct-sum qtf_potential_sum.csv [--exclude LO HI]
python qtf_compare.py RESULTS --indirect-sum qtf_potential_sum.csv [--exclude LO HI]
python qtf_compare.py RESULTS --motion qtf_motion.csv [qtf_motion_sum.csv] [--exclude LO HI]
python qtf_compare.py RESULTS --against-nemoh FOLDER [--hull oc4] [--reading as-is|conj|-1|-conj] [--max-dw 0.3] [--exclude LO HI] [--kernel diff|sum]
python qtf_compare.py RESULTS --against-wamit FILE [--hull oc4w] [--reading as-is|conj|-1|-conj] [--max-dw 0.3] [--max-omega 1.90] [--exclude LO HI]
```

- Run `qtf_compare.py` without `-I`: it imports `qtf_tools` from its own folder, which `-I` leaves off the path.
- `check` (or `--self-test`) runs the self-checks of `--potential` (also on a `qtf_terms = full` fixture, `fs_columns`), `--against-nemoh` and `--against-wamit` on temporary fixtures, those of `--indirect-sum` (`indirect_sum_reader`) and `--against-nemoh --kernel sum` (`nemoh_sum_reader`), and that of `--motion` (`motion_compare`).
- `--direct` compares the direct solve with the indirect one, case by case, for the `q_*direct*` cases (`qtf_direct = 1`): T_dir (`QTF_Direct.12d` − `.10d`, or the parts of `QTF_DirectTerms.txt`) against T_pot (`QTF_Indirect.11d` − `.10d`) of the same run, as |T_dir − T_pot| / max |T_pot| over the forces and over the moments (max and median), and the number of difference frequencies kept by `--max-dw`. `--by-omega` adds the largest error at each difference frequency, where an irregular frequency shows as a spike. With the CSV of `qtf-potential` the sharp cases also get `.12d` and `.11d` against the exact reference. For the `q_*_directfs_*` cases (`qtf_terms = full`, three parts D_I2, D_D, D_FS in `QTF_DirectTerms.txt`, seven in `QTF_PotentialTerms.txt`) the reference must come from `qtf-potential --free-surface` (else the case is noted `no FS reference`), and a row `FS` compares FS_inner + FS_outer of the indirect terms with the reference column FS; the free-surface columns of the direct file are D_FS. The number of parts is read from the file, with a note if it differs from the case name. `drift_convergence.py` makes these cases (`qtf_directfs_cases`: `q_round_directfs_free_m*` and `q_sharp_directfs_free_m*`).
- `--direct-sum REF` compares the sum-frequency direct solve of the `q_*directsum*` cases (`qtf = 3`, `qtf_direct = 1`, `qtf_terms = full`) with the CSV of `qtf-potential --kernel sum --free-surface`: `QTF_Direct.12s` against the exact total, `.12s` − `.10s` against the exact potential part, and from `QTF_DirectTermsSum.txt` D_FS against the exact free-surface part and P_BODY + P_FS against the potential part (the file holds both orders of an unequal pair, read in either order). A last table is the in-run check, D = D_I2 + D_D + D_FS against P = P_BODY + P_FS row by row. `check` tests it on a fixture. `drift_convergence.py` makes the cases (`qtf_directsum_cases`: `q_round_directsum_free_m*` and `q_sharp_directsum_free_m*`, IRSP = 1, differencing).
- `--reference` compares `QTF_Quadratic.10d` and `.10s` with the CSV of `truncated_cylinder_reference.py qtf`; `--potential` compares `QTF_Indirect.11d`, `.11d` − `.10d` and each term of `QTF_PotentialTerms.txt` with the CSV of `qtf-potential`. For a `q_sharp_full_*` case (`qtf_terms = full`) the reference total and potential part include the free-surface integral, and a row `FS` compares `FS_inner` + `FS_outer` of the 19-column terms file with the `Fx_FS_*` … columns of a CSV made with `qtf-potential --free-surface` (the split between inner and outer depends on the mesh and R, so only the sum is compared; the scale is the larger of FS and the potential part, as for the other terms). Without those columns the case prints a note. `--nemoh` adds Nemoh's `OUT_QTFM_N.dat` against the same reference.
- `--indirect-sum REF` compares `QTF_Indirect.11s` of the `q_*indsum*` cases (`qtf = 3`, `qtf_terms = full`, `qtf_direct = 0`) with the CSV of `qtf-potential --kernel sum --free-surface`: `.11s` against the exact total, `.11s` − `.10s` against the exact potential part, and from `QTF_PotentialTermsSum.txt` (FK, ID, B1, B2, B3, FS_inner, FS_outer; both orders of an unequal pair) the sum of the terms against the potential part and FS_inner + FS_outer against the exact free-surface part, as errors over the pairs ω_i ≤ ω_j relative to the largest exact value. A case that also has `QTF_Direct.12s` gets the row `.11s-12s`, |.11s − .12s| over the largest |.12s|. `--exclude LO HI` leaves a band out, as for the other tables.
- `--motion REF [REF]` compares `QTF_Motion.11d` and `.12d` (with the `diff` rows) and `.11s` and `.12s` (with the `sum` rows) of the `q_*motion*` cases with the CSVs of `truncated_cylinder_reference.py qtf-motion`, surge, heave and pitch, as errors over the pairs ω_i ≤ ω_j relative to the largest exact value, with the in-run row `.11-12`, |.11 − .12| over the largest |.12|. The files hold conj(ξ) L and conj(ξ) L² (L = 1 in these cases) after `#` header lines, and are conjugated back. The pairs a file leaves out (an Ω whose matrix has rcond < 1e-12, [the static diagonal and mooring](input_files.md#second-order-motion), as the difference diagonal of a hull moored in surge only, where sway and yaw have no restoring) are not compared: the cell shows (compared/pairs), and one line per file counts them with its `# Left out` lines. `--exclude LO HI` leaves a band out, off by default; use it for the undamped first-order surge-pitch resonance of the free cylinder (0.75 to 0.80 rad/s). The reference keeps the M and C of the exact circle, while a HAMS polygon has slightly lower M and C, which shifts the response near the 100 s resonance. `check` tests it (`motion_compare`).
- `--against-nemoh FOLDER --kernel sum` compares `QTF_Indirect.11s` of each `q_<hull>_indsum_free_m<k>` case with `OUT_QTFP_N.dat` of `nemoh_<hull>_free_fs_m<k>` (the full sum-frequency QTF) over every pair ω_i ≤ ω_j, the diagonal included. Nemoh writes the lower triangle of a symmetric matrix for the sum frequency, so the other half is the same value, with no conjugate. The reading of Nemoh's value is `as-is` unless `--reading` says otherwise, never chosen by agreement.
- `--against-nemoh FOLDER` compares `QTF_Indirect.11d` of each `q_<hull>_<motion>_m<k>` case with Nemoh's `nemoh_<hull>_<motion>_m<k>` run under FOLDER, when there is no exact reference. `--hull` names the hull (default `oc4`). `--reading` sets how Nemoh's value is read: `as-is` (default), its conjugate, or either with the sign changed. The reading is never chosen by agreement with HAMS; the table warns when another reading would agree much better.
- `--against-wamit FILE` compares `QTF_Indirect.11d` and `QTF_Quadratic.10d` of each `q_<hull>_*` case with a WAMIT numeric `.12d` (difference QTF), such as `marin_semi.12d` of the OpenFAST r-test. `--hull` names the hull (default `oc4w`). The `.12d` is read as written, with the periods turned into ω and the value conjugated back like the HAMS files; `--reading` is as for `--against-nemoh`, and is never chosen by agreement. It prints, per case and file, the largest and median error over the pairs with ω ≤ `--max-omega` for three sets: all off-diagonal pairs, those with 0 < ω_j − ω_i ≤ `--max-dw`, and the diagonal (the mean drift), with the median ratio Re(HAMS) / Re(WAMIT) on the diagonal. The diagonal fixes the sign but not the conjugation, since it is real.
- `--max-omega` limits the pairs to ω_i and ω_j ≤ max-omega (default 1.90 rad/s, the top of the `q_oc4w` grid).
- `--exclude LO HI` leaves out the pairs with a frequency in [LO, HI] rad/s, such as a resonance band, and keeps the scale of all pairs.
- `--max-dw` limits the difference pairs to ω_j − ω_i ≤ max-dw (default 0.3 rad/s).

## Partial Second-Order Field Maps

**`SecondOrder/qtf_field_plots.py`** maps `QTF_FieldPartial.txt` (format in [Output Files](output_files.md#wamit-format)) and prints grids of field points for `ControlFile.in`. It needs `numpy` and `matplotlib`. The field in the file is partial: it omits the second-order diffracted wave and, for a free body, the wave radiated by its second-order motion. Near the bodies, where the free-surface forcing drives the second-order scattered wave, the omitted part is not small, so the maps say PARTIAL in their title. The complete field is a later phase (see [Partial second-order field](input_files.md#partial-second-order-field-at-the-field-points)).

```
python -I qtf_field_plots.py plot Output/Wamit_format/QTF_FieldPartial.txt --out maps [--kernel diff|sum|both] [--type pressure|elevation|both]
python -I qtf_field_plots.py grid --x X0 X1 NX --y Y0 Y1 NY [--z 0] [--exclude XC YC R]...
python -I qtf_field_plots.py check
```

- `plot` writes one PNG per kernel, frequency pair, heading pair, type and depth (pressure points are grouped by z), named `field_<diff|sum>_<PER1>_<PER2>_<beta1>_<beta2>_<elevation|pressure_z<z>>.png`, with `m` for a minus sign. Each has two panels over x and y, in the units of the file. The modulus is contoured between the points, or shown as coloured markers when there are fewer than 3 points off one line, two points at the same (x, y), or one value. The phase (degrees, cyclic colour map) is always shown as markers at the points, since a contour would interpolate across the jump from 180 to −180. The number of PNGs grows with the pairs; `--kernel` and `--type` select. The title gives the first two columns as PER1 and PER2, as in the file; they are periods or frequencies as `Output_frequency_type` says.
- `grid` prints the field-point block: a `Number_of_field_points` line, then one `x y z` line per point of the NX by NY grid, leaving out the points within R of each (XC, YC). Make the discs larger than the bodies, or points on or inside a hull stop the run.
- `check` writes a synthetic file, reads, plots and grids it, and exits with 1 if one of its 9 checks fails.
