# Example 07: Brownian Rouse-chain benchmark

This example simulates noninteracting Gaussian chains with overdamped Brownian dynamics.

## Learning goals

- Use `atom_style molecular` and harmonic bonds.
- Run overdamped Brownian dynamics with `fix brownian`.
- Avoid velocity-based diagnostics when velocities are not physically meaningful.
- Compare simulation results with Rouse-model reference values.

## Files

| File or folder | Purpose |
|---|---|
| `in.main` | Main Brownian dynamics input |
| `in.data` | Initial chain configuration |
| `preproc/` | Generator for `in.data` |
| `postproc_full_from_dump/` | Full trajectory analysis from the dump |
| `postproc_msd_from_dump/` | Alternative MSD postprocessing route |

## Run

```bash
./run.sh
```

Then run the full postprocessing workflow:

```bash
cd postproc_full_from_dump
./run.sh
```
