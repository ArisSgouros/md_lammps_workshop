#!/usr/bin/env python3
"""Generate Gaussian initial conformations for the Brownian Rouse-chain example."""

import numpy as np

# -------------------------
# Parameters
# -------------------------

n_chains = 128
n_beads = 8
b = 1.0
box = 100.0
margin = 8.0
seed = 12345

out = "in.data"

rng = np.random.default_rng(seed)

atoms = []
bonds = []

atom_id = 1
bond_id = 1

for mol_id in range(1, n_chains + 1):

    # Gaussian random walk.
    # Each bond vector has <r^2> = b^2.
    steps = rng.normal(
        loc=0.0,
        scale=b / np.sqrt(3.0),
        size=(n_beads - 1, 3)
    )

    coords = np.zeros((n_beads, 3))
    coords[1:] = np.cumsum(steps, axis=0)

    # Center the chain near zero, then place it inside the box.
    coords -= coords.mean(axis=0)

    low = -coords.min(axis=0) + margin
    high = box - coords.max(axis=0) - margin

    if np.any(high <= low):
        raise RuntimeError("Box too small or margin too large for generated chain.")

    shift = rng.uniform(low, high)
    coords += shift

    for i in range(n_beads):
        if i == 0:
            atom_type = 2       # first bead
        elif i == n_beads - 1:
            atom_type = 3       # last bead
        else:
            atom_type = 1       # internal bead

        x, y, z = coords[i]
        atoms.append((atom_id, mol_id, atom_type, x, y, z))

        if i > 0:
            bonds.append((bond_id, 1, atom_id - 1, atom_id))
            bond_id += 1

        atom_id += 1

n_atoms = len(atoms)
n_bonds = len(bonds)

with open(out, "w") as f:
    f.write("Rouse chains: Gaussian initial conformations\n\n")

    f.write(f"{n_atoms} atoms\n")
    f.write(f"{n_bonds} bonds\n")
    f.write("\n")

    f.write("3 atom types\n")
    f.write("1 bond types\n")
    f.write("\n")

    f.write(f"0.0 {box:.8f} xlo xhi\n")
    f.write(f"0.0 {box:.8f} ylo yhi\n")
    f.write(f"0.0 {box:.8f} zlo zhi\n")
    f.write("\n")

    f.write("Masses\n\n")
    f.write("1 1.0\n")
    f.write("2 1.0\n")
    f.write("3 1.0\n")
    f.write("\n")

    f.write("Atoms\n\n")
    for a in atoms:
        atom_id, mol_id, atom_type, x, y, z = a
        f.write(f"{atom_id} {mol_id} {atom_type} {x:.8f} {y:.8f} {z:.8f}\n")

    f.write("\nBonds\n\n")
    for bnd in bonds:
        bond_id, bond_type, i, j = bnd
        f.write(f"{bond_id} {bond_type} {i} {j}\n")

print(f"Wrote {out}")
print(f"Chains: {n_chains}")
print(f"Beads per chain: {n_beads}")
print(f"Atoms: {n_atoms}")
print(f"Bonds: {n_bonds}")
