# LAMMPS workshop examples

This folder contains the example scripts for the workshop. The examples are ordered from the smallest possible LAMMPS input to more realistic simulation and analysis workflows.

## Organization

| Folder | Topic |
|---|---|
| `00_lammps_intro_core` | Minimal complete MD workflow |
| `01_lammps_intro_command_structure` | LAMMPS command grammar, groups, fixes, and dumps |
| `02_timestep` | Timestep choice and NVE energy conservation |
| `03_computes` | MSD, diffusion, RDF, and averaged output |
| `04_nve_nvt_npt_ensembles` | NVE, NVT, NPT, and anisotropic barostatting |
| `05_phase_behavior` | Binary LJ mixture and partial RDFs |
| `06_free_film` | Liquid slab, density profile, and pressure-tensor surface tension |
| `07_rouse_chain` | Brownian Rouse-chain benchmark |

## Commenting style

The first examples are intentionally commented heavily because they introduce the script structure and command grammar. Later examples use the same conventions but avoid repeating explanations that have already been introduced.

## ID convention

Short capital IDs are used consistently:

- `G...` for groups, for example `G1`, `GA`, `GB`
- `F...` for fixes, for example `FNVE`, `FNVT`, `FRDF`
- `C...` for computes, for example `CMSD`, `CRDF`, `CPRESS`
- `D...` for dumps, for example `DALL`, `DEQ`, `DPROD`

## File naming convention

Input files use the `in.*` prefix. Generated outputs use the `o.*` prefix. Logs are written as `log.<tag>.lammps`.

Typical output names are:

- `o.dump.all.lammpstrj`
- `o.dump.prod.lammpstrj`
- `o.data.final`
- `o.msd.dat`
- `o.rdf.dat`

## Running an example

From inside an example folder:

```bash
./run.sh
```

To remove generated output:

```bash
./clean.sh
```

The scripts assume that the LAMMPS executable is available as `lmp_mpi`. Change `run.sh` if your executable has a different name.
