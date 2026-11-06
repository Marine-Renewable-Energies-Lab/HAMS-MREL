# HAMS-MREL

HAMS-MREL (Hydrodynamic Analysis of Marine Structures – Marine Renewable Energies Lab) is an open-source numerical tool for studying wave–structure interactions involving multiple floating bodies. It uses a Boundary Integral Equation Method to compute hydrodynamic coefficients, wave excitation forces, and related quantities for diffraction and radiation problems. The software includes OpenMP-based parallelization to improve computational efficiency and has been validated against WAMIT for a range of multi-body test cases.

HAMS-MREL is developed as an extension of the original [HAMS solver](https://github.com/YingyiLiu/HAMS). While HAMS focuses on single-body hydrodynamic analysis, HAMS-MREL expands the framework to handle interactions among multiple bodies. Other features include free surface elevation and pressure fields, suppression of irregular frequencies, global symmetry, generalized modes, an optional GMRES iterative solver with H-matrix acceleration for large multi-body problems, and integration with WEC-Sim.

Current version: 2.3.3 (The repository has migrated from https://github.com/vaibhavraghavan/HAMS-MREL to keep all open source tools from MREL under one unbrella. For v2.0 to v2.2, got to the old repository.)

The binary for Windows is available in the [Releases](https://github.com/Marine-Renewable-Energies-Lab/HAMS-MREL/releases). This was compiled using OneAPI 2022 (ifort) and Visual Studio 2019. 

## Release notes

v2.3.3 corrects these errors in first-order results:

- **Free-surface image term.** Now integrated over the source panel instead of taken at its centroid. Results change a little on all meshes; both forms converge to the same limit. See [the image term section](docs/developer_documentation.md#free-surface-image-term).
- **Triangular panels.** A triangle with an obtuse angle at node 2 is now integrated over the correct shape. Other meshes are not affected. See [triangle panels](docs/developer_documentation.md#triangle-panels).
- **Generalized modes in DOF 1 to 3.** The panel direction was read from the columns the mode overrides, which changed all modes of all bodies in v2.0 to v2.3.2. Runs without a mode file do not change. See [gen_mod](docs/input_files.md#gen_mod-optional).
- **Finite depth, panels closer than 2% of the depth.** Now computed with Linton's method up to ν·h = 50, which removes jumps in the Green function. Infinite depth is not affected.
- **Finite depth at large kh.** Results above kh = 355 are no longer NaN. Below kh = 300 the output is unchanged.

It also adds two checks:

- A run with NaN or Inf results says so and exits with code 1.
- A hull panel with its centroid on z = 0 stops the run with a message.

## Installation and Usage

Installation instructions are provided for:
- [Windows](docs/installation_windows.md)
- [Linux (Ubuntu)](docs/installation_linux.md)
- [High Performance Computing Clusters](docs/installation_linux.md)

HAMS-MREL requires input files to be placed in an input directory before execution. The program can be run in two ways, depending on whether the input and output directory locations are specified explicitly.

### 1. Default mode (Input/Output directories next to the executable)

In this mode, the program expects two directories named **Input** and **Output** to be located in the same folder as the executable.

**Windows (Command Prompt or Power Shell)**
```
hams-mrel.exe
```

**Linux (Terminal)**
```
./hams-mrel
```

### 2. Specifying input and output directory paths

The user may also provide paths to the input and output directories. These directories may be located anywhere on the system. In this mode, the program expects exactly two arguments: the input directory and the output directory. Providing only one argument is not allowed and will result in an error.

**Windows**
```
hams-mrel.exe <path-to-input-dir> <path-to-output-dir>
```

**Linux**
```
./hams-mrel <path-to-input-dir> <path-to-output-dir>
```

## Input and Output Files

### Input

HAMS-MREL requires a set of input files to run. These files contain the general setting for the simulation (`ControlFile.in`), hydrostatic data (`Hydrostatic[_i].in`) and meshes for the submerged bodies (`HullMesh[_i].pnl`) and inner water plane for each body (`WaterPlaneMesh[_i].pnl`). The exact number of files depends on the number of bodies and on whether irregular frequencies are analyzed.  

**A full description of each inout file is provided in [docs/input_files.md](docs/input_files.md)**.


### Output

HAMS-MREL produces output files in multiple established formats to support post-processing and compatibility with external tools. Every simulation outputs data in the HAMS and WAMIT formats (plus Hydrostar for single-body cases), covering added mass, damping, excitation forces, motions, pressure/elevation, and hydrostatics, along with an _ErrorCheck.txt_ file summarizing key system parameters. Output formats are fixed and always generated.

**A full description of each file is provided in [docs/output_files.md](docs/output_files.md)**.

## Developer Documentation

Additional developer-focused documentation is available in [`docs/developer_documentation.md`](docs/developer_documentation.md). This includes details on the test suite, the GitHub Actions continuous integration pipeline and results from memory profiling analysis.

## How to cite

If you use HAMS‑MREL in your research, please cite:

> Raghavan, V., Loukogeorgaki, E., Mantadakis, N., Metrikine, A. V., & Lavidas, G. (2024). _HAMS‑MREL, a new open‑source multiple body solver for marine renewable energies: Model description, application and validation._ Renewable Energy, 237, 127857. [https://doi.org/10.1016/j.renene.2024.121577]

> Raghavan, V., Metrikine, A.V. & Lavidas, G. Theory and validation of the new features in BIEM solver HAMS-MREL. J. Ocean Eng. Mar. Energy 12, 53–71 (2026). [https://doi.org/10.1007/s40722-025-00431-8]

## License

HAMS-MREL is available under the Apache 2.0 license. 

## Contact

HAMS‑MREL was developed by Vaibhav Raghavan. For questions, please contact the author at vaibhav.raghavan@gmail.com. Users are encouraged to open a new issue in this repository to report bugs, request new features, or suggest improvements (Using Discussions).

This work is part of the TU Delft Marine Renewable Energies Lab (MREL) led by George Lavidas (G.Lavidas@tudelft.nl). Learn more about the laboratory [here](https://www.tudelft.nl/en/ceg/about-faculty/departments/hydraulic-engineering/sections/offshore-engineering/research/marine-renewable-energies-lab-mrel). 

HAMS-MREL was supported by Yasel Quintero Lares from the Digital Competence Centre, Delft University of Technology.  

## Publications utilizing HAMS-MREL

> Raghavan, V., Loukogeorgaki, E., Mantadakis, N., Metrikine, A. V., & Lavidas, G. (2024). _HAMS‑MREL, a new open‑source multiple body solver for marine renewable energies: Model description, application and validation._ Renewable Energy, 237, 127857. [https://doi.org/10.1016/j.renene.2024.121577]

> Raghavan, V., Metrikine, A.V. & Lavidas, G. Theory and validation of the new features in BIEM solver HAMS-MREL. J. Ocean Eng. Mar. Energy 12, 53–71 (2026). [https://doi.org/10.1007/s40722-025-00431-8]

> Raghavan, V., Simonetti, I., Metrikine, A. V., Lavidas, G., & Cappietti, L. (2024). A new numerical modelling framework for fixed oscillating water column wave energy conversion device combining BEM and CFD methods: Validation with experiments. Ocean Engineering, 301, 117543. [https://doi.org/10.1016/j.oceaneng.2024.117543]

> Lavidas, G., Mezilis, L., Alday G., M., Baki, H., Tan, J., Jain, A., Engelfried, T., & Raghavan, V. (2025). Marine renewables in Energy Systems: Impacts of climate data, generators, energy policies, opportunities, and untapped potential for 100% decarbonised systems. Energy, 336, 138359. [https://doi.org/10.1016/j.energy.2025.138359]

> Asiikkis, A. T., Raghavan, V., Grigoriadis, D. G. E., Metrikine, A. V., Lavidas, G., & Vakis, A. I. (2025). Experimental validation of time domain simulations with HAMS-MREL and comparative analysis of linear and weakly nonlinear models for dense WEC arrays. Ocean Engineering, 341, 122553. [https://doi.org/10.1016/j.oceaneng.2025.122553]

> Alday, M., Raghavan, V., & Lavidas, G. (2025). Effects of wave spectrum representation on power production estimations from point absorbers. Applied Ocean Research, 161, 104626. [https://doi.org/10.1016/j.apor.2025.104626]

> Engelfried, T., Alday, M., Raghavan, V., & Lavidas, G. (2026). Geospatial analysis of the life cycle global warming impacts from marine renewables. Renewable and Sustainable Energy Reviews, 226, 116338. [https://doi.org/10.1016/j.rser.2025.116338]

> Raghavan, V., Metrikine, A. V., & Lavidas, G. (2026). Mixed Wave Energy arrays, re-thinking wave array configurations: A hydrodynamics based approach. Renewable Energy, 272, 126005. [https://doi.org/10.1016/j.renene.2026.126005]

Publications to follow...
