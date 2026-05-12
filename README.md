# MD LAMMPS Workshop

Introductory molecular-dynamics workshop using **LAMMPS**.

The workshop starts from the basic structure of a LAMMPS input script and gradually builds toward simple but useful molecular-dynamics examples: timestep stability, ensembles, computes, phase behavior, free films, and polymer-chain dynamics.

Repository:

```text
https://github.com/ArisSgouros/md_lammps_workshop/
```

---

## Goals

By the end of the workshop, participants should be able to:

- read and modify basic LAMMPS input scripts;
- run simple MD simulations from the command line;
- understand the role of common commands such as `units`, `read_data`, `pair_style`, `velocity`, `fix`, `compute`, `dump`, and `run`;
- compare NVE, NVT, NPT, and anisotropic pressure-control simulations;
- compute and interpret basic observables such as MSD, RDF, density profiles, and pressure-tensor quantities;
- connect simple model systems to more realistic MD applications.


---

## Repository structure

```text
md_lammps_workshop/
├── LICENSE
├── README.md
├── examples/
├── slides/
└── utilities/
```

Main folders:

```text
examples/    Runnable LAMMPS examples
slides/      Presentation material
utilities/   Scripts for postprocessing atomistic trajectories
```

---

## Examples

The examples are organized as separate folders:

```text
00_lammps_intro_core
01_lammps_intro_command_structure
02_timestep
03_computes
04_nve_nvt_npt_ensembles
05_phase_behavior
06_free_film
07_rouse_chain
```

| Folder | Topic |
|---|---|
| `00_lammps_intro_core` | Minimal LAMMPS run |
| `01_lammps_intro_command_structure` | Command IDs, groups, fixes, computes, dumps |
| `02_timestep` | NVE timestep stability and energy conservation |
| `03_computes` | Basic observables and LAMMPS computes |
| `04_nve_nvt_npt_ensembles` | NVE, NVT, NPT, and anisotropic pressure control |
| `05_phase_behavior` | LJ fluids, binary mixtures, RDFs, and demixing trends |
| `06_free_film` | Density profile and surface tension of a free LJ film |
| `07_rouse_chain` | Polymer-chain dynamics and Rouse-model benchmarks |

Each example is intended to be self-contained.

---

### Windows

LAMMPS provides precompiled Windows executables:

<https://docs.lammps.org/Install_windows.html>

A practical alternative is to use **Windows Subsystem for Linux (WSL)**. In that case, install LAMMPS through the Linux/Conda route:

<https://docs.lammps.org/Install_conda.html>

### Source build

LAMMPS can also be built from source. This is useful when custom packages, additional modules, or code development are required.

<https://docs.lammps.org/Install_linux.html> <br>

<https://docs.lammps.org/Build.html>


---

## Running an example

Move into an example folder:

```bash
cd examples/02_timestep
```

Run the main input file:

```bash
lmp -in in.main
```

For an MPI-enabled build:

```bash
mpirun -np 4 lmp -in in.main
```

---

## Naming convention

Each runnable example uses:

```text
in.main
```

Typical files inside an example folder may include:

```text
README.md
in.main
in.data
log.main
dump.main.lammpstrj
data.final
restart.final
```

Analysis files use descriptive names, for example:

```text
msd.dat
rdf.dat
density_profile_z.dat
surface_tension.dat
```

---

## Suggested workflow

For each example:

1. Read the local `README.md`.
2. Inspect `in.main`.
3. Run the example.
4. Check the thermo output.
5. Visualize the trajectory if a dump file is written.
6. Modify one parameter.
7. Rerun and compare.

Change one parameter at a time.

---

## Recommended tools

Required:

- LAMMPS
- Python 3

Useful:

- VMD for trajectory visualization

---

## References

Useful starting points:

- LAMMPS documentation: https://docs.lammps.org/
- LAMMPS website: https://www.lammps.org/
- M. P. Allen and D. J. Tildesley, *Computer Simulation of Liquids*
- D. Frenkel and B. Smit, *Understanding Molecular Simulation*
- M. Doi and S. F. Edwards, *The Theory of Polymer Dynamics*
- M. Rubinstein and R. H. Colby, *Polymer Physics*

---

## License

See `LICENSE`.

This repository is intended as teaching material for introductory molecular-dynamics simulations with LAMMPS.

