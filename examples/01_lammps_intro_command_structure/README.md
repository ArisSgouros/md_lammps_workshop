# Example 01: LAMMPS command structure

This example emphasizes LAMMPS command grammar, object IDs, and how groups, fixes, and dumps are connected.

## Learning goals

- Understand the pattern `command ID group-ID style arguments`.
- Define atom groups by type.
- Write separate dumps for the full system and for one group.
- Cleanly remove fixes and dumps before writing the final data file.

## Files

| File | Purpose |
|---|---|
| `in.main` | Main teaching script |
| `in.data` | Initial two-component LJ configuration |
| `run.sh` | Runs the example |
| `clean.sh` | Removes generated output files |
| `preproc/make_data.py` | Regenerates `in.data` if needed |

## Run

```bash
./run.sh
```

Expected outputs include `log.command_structure.lammps`, `o.dump.all.lammpstrj`, `o.dump.g1.lammpstrj`, and `o.data.final`.
