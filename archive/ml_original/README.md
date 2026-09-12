# Recovered original machine-learning sources

These source and aggregate files are preserved byte-for-byte. Their SHA256
hashes are recorded in `sha256.json`; the active convention-adapted versions
are in `scripts/ml/`.

| Archived file | Recovered source |
| --- | --- |
| `run_ml_recovery.py` | User attachment `run_ml_recovery(1).py`, byte-identical to the original Library file |
| `run_class6_ml_discovery.py` | `class56_lax_discovery(1).zip` |
| `class5_lax.py` | Shared Class 5/6 helper in the same ZIP |
| `run_ml_discovery.py` | Earlier Class 5 search in the same ZIP |
| `class6_ml_summary.json` | Historical aggregate in the same ZIP |
| `class5_ml_aggregate.json` | `class5_class6_landau_lifshitz_complete_calculation_notebook_provenance_audited.zip`, `results/class5_ml_aggregate.json` |

The uploaded ZIP has SHA256
`bf4c98435001c3cbddb14be51b99a35ebe0fff5870e93bc34059b15008a1d77a`
and matches the recovered Library archive byte-for-byte. The organized
provenance package changed comments, output labels and an import path in the
training scripts; its numerical searches agree with these originals.

Class 5 optimizes twelve temporal coefficients plus `beta` with a uniform
sparsity penalty and an off-shell residual. Class 6 optimizes eight
coefficients with on-shell curvature. Both use Adam followed by L-BFGS.
The original Class 6 on-shell objective is retained as an option in the active
code; its default off-shell extension is exactly equivalent for every choice
of coefficients in this ansatz.

The convention verification imports these originals directly, compares
residuals/losses/gradients at arbitrary coefficients, and independently checks
the final matrices against the manuscript. The active scripts preserve the
original optimizer settings and spatial sampling. See `REPRODUCIBILITY.md`.
