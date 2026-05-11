# Postprocessing: full Rouse-chain analysis from one dump

This folder analyzes the trajectory written by the parent Brownian Rouse-chain example.

It reads `../o.dump.all.lammpstrj` and computes:

- chain center-of-mass MSD using multiple time origins
- radius of gyration
- end-to-end vector size
- end-to-end vector autocorrelation
- fitted chain diffusion coefficient
- fitted longest relaxation time

Run after completing the parent example:

```bash
./run.sh
```
