# Example 02: Timestep and energy conservation

This example shows how the integration timestep affects energy conservation in an NVE simulation.

## Learning goals

- Run the same LJ fluid with different timesteps.
- Monitor total-energy drift, `dE = E(t)-E(0)`.
- Identify when a timestep is conservative, marginal, or too large.

## Suggested exercise

Change this line in `in.main`:

```lammps
variable        dt equal 0.005
```

Try `0.001`, `0.005`, `0.02`, and `0.05`, then compare the logs.

## Run

```bash
./run.sh
```

Generated files include `log.dt<dt>.lammps` and `o.dump.dt<dt>.lammpstrj`.
