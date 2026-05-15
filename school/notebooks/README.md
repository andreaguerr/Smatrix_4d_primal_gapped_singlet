# School Notebooks

This folder contains a small Z2-even 4D S-matrix bootstrap example for the
school material.

## Files

- `00_run_Z2even_simple.wl`: laptop runner for the example.
- `Z2even_Bootstrap_4D_simple.wl`: simple bootstrap definitions.
- `Z2even_Bootstrap_4D_simple_threshold_option.wl`: variant with the threshold
  bound-state option.
- `SDPB.m`: Mathematica exporter used to write SDPB input.

## Integral Cache

The precomputed integral cache is not stored in git. The runner can download it
from the public release when required:

```text
https://github.com/andreaguerr/Smatrix_4d_primal_gapped_singlet/releases/tag/integrals-refined-4d-v2026-05-15
```

The cache is installed locally as `integral_storage_bin/`, which is ignored by
git together with generated SDPB run output.
