# Postprocessing: RDF from dump

This folder computes the radial distribution function from `../o.dump.prod.lammpstrj` using the external `pyrdf.py` utility.

Run after completing the parent example:

```bash
./run.sh
```

The file `in.data_init_mod.dat` provides the corresponding structure information required by the utility.

## Form-factor note

The file `form_factors.dat` follows the convention used by the RDF/XRD utility. The original note in this folder stated that the form factors were taken from the International Tables for Crystallography and are in electron/Angstrom units.
