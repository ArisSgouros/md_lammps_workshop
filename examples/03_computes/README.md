# Example 03: Computes and averaged output

This example demonstrates LAMMPS computes and time-averaged output using an LJ fluid at `T* = 1.0` and `rho* = 0.8`.

## Learning goals

- Use `compute msd` to calculate mean-square displacement.
- Estimate self-diffusion from `MSD(t) ≈ 6Dt`.
- Use `compute rdf` and `fix ave/time` to write `g(r)`.
- Compare on-the-fly analysis with postprocessing utilities.

## Run

```bash
./run.sh
```

Main outputs:

- `o.msd.dat`
- `o.rdf.dat`
- `o.dump.prod.lammpstrj`
- `o.data.init`

The subfolders contain optional postprocessing workflows.
