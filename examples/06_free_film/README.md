# Example 06: LJ free film

This example builds a liquid slab in a larger periodic box and estimates density profiles and surface tension.

## Learning goals

- Build a slab by filling only part of the simulation box.
- Average a number-density profile along `z` using `chunk/atom` and `ave/chunk`.
- Estimate mechanical surface tension from the pressure tensor.

## Main formula

For a slab with two interfaces:

```text
gamma = Lz/2 * (P_N - P_T)
P_N = Pzz
P_T = (Pxx + Pyy)/2
```

## Run

```bash
./run.sh
```

Main outputs are `o.density_profile_z.dat` and `o.surface_tension_pressure_tensor.dat`.
