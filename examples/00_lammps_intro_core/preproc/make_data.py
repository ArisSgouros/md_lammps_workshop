import numpy as np

rho = 0.8
ncell = 6
xB = 0.5
seed = 12345

rng = np.random.default_rng(seed)

# FCC basis in fractional coordinates
basis = np.array([
    [0.0, 0.0, 0.0],
    [0.5, 0.5, 0.0],
    [0.5, 0.0, 0.5],
    [0.0, 0.5, 0.5],
])

# FCC number density: rho = 4 / a^3
a = (4.0 / rho) ** (1.0 / 3.0)
L = ncell * a

atoms = []
atom_id = 1

for i in range(ncell):
    for j in range(ncell):
        for k in range(ncell):
            origin = np.array([i, j, k], dtype=float) * a
            for b in basis:
                pos = origin + b * a
                atom_type = 2 if rng.random() < xB else 1
                atoms.append((atom_id, atom_type, *pos))
                atom_id += 1

with open("in.data", "w") as f:
    f.write("Binary LJ intro data file\n\n")
    f.write(f"{len(atoms)} atoms\n\n")
    f.write("2 atom types\n\n")
    f.write(f"0.0 {L:.8f} xlo xhi\n")
    f.write(f"0.0 {L:.8f} ylo yhi\n")
    f.write(f"0.0 {L:.8f} zlo zhi\n\n")

    f.write("Masses\n\n")
    f.write("1 1.0\n")
    f.write("2 1.0\n\n")

    f.write("Atoms # atomic\n\n")
    for atom_id, atom_type, x, y, z in atoms:
        f.write(f"{atom_id} {atom_type} {x:.8f} {y:.8f} {z:.8f}\n")

print("Wrote data.lj_intro")
print(f"Number of atoms: {len(atoms)}")
print(f"Box length: {L:.6f}")
