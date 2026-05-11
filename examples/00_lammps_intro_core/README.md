# Example 00: Minimal LAMMPS core workflow

This example introduces the smallest complete LAMMPS molecular-dynamics workflow used in the workshop.

## Learning goals

- Recognize the standard blocks of a LAMMPS input file.
- Read a data file with `read_data`.
- Define a Lennard-Jones pair potential.
- Assign masses and velocities.
- Run a short NVE simulation.
- Write a trajectory dump.

## Files

| File | Purpose |
|---|---|
| `in.main` | Heavily commented teaching version |
| `in.main_raw` | Compact version of the same workflow |
| `in.data` | Initial two-component LJ configuration |
| `run.sh` | Runs `in.main` with LAMMPS |
| `clean.sh` | Removes generated output files |
| `preproc/make_data.py` | Regenerates `in.data` if needed |

## Run

```bash
./run.sh
```

Expected outputs include `log.intro_core.lammps` and `o.dump.all.lammpstrj`.
