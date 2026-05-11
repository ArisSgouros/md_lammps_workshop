# Example 05: Binary LJ mixture and partial RDFs

This example uses a binary Lennard-Jones mixture to illustrate how interaction parameters affect local structure and mixing tendencies.

## Learning goals

- Create a random binary mixture using `set type/fraction`.
- Define different `A-A`, `A-B`, and `B-B` LJ interactions.
- Compute partial RDFs: `g_AA(r)`, `g_AB(r)`, and `g_BB(r)`.

## Suggested exercise

Change:

```lammps
variable        eps_ab equal 1.0
```

to a smaller value, for example `0.5`, and compare the RDFs and trajectory.

## Run

```bash
./run.sh
```

Main outputs are `o.dump.all.lammpstrj`, `o.rdf.dat`, and `o.data.final`.
