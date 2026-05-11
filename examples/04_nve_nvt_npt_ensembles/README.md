# Example 04: NVE, NVT, NPT, and anisotropic barostatting

This example compares common ensembles using the same LJ fluid.

## Learning goals

- Distinguish NVE, NVT, and NPT behavior.
- Observe how thermostatting affects energy conservation.
- Observe how isotropic NPT changes all box lengths together.
- Use `z` barostatting to keep `Lx` and `Ly` fixed while controlling `Pzz`.

## Run

```bash
./run.sh
```

Watch the thermo columns `pxx`, `pyy`, `pzz`, `density`, `lx`, `ly`, and `lz`.
