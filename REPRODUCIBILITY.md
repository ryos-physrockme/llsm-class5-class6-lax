# Reproducibility guide

All commands below are run from the repository root after installing `requirements.txt`.

## Class 5: quantum-to-classical time-Lax limit

Exact symbolic derivation:

```bash
python scripts/class5/verify_coherent_symbolic.py
```

This checks the fundamental operator identities, Sutherland relation, coherent-state product correction, corrected time-Lax matrix, Hamiltonian flow, off-shell EOM factorization, and the undeformed XXX limit. Results are written to `results/class5/symbolic_verification.json` and `results/class5/symbolic_expressions.json`.

Finite-spacing and finite-lattice checks:

```bash
python scripts/class5/verify_coherent_numeric.py
```

Outputs are `results/class5/numeric_verification.json` and `results/class5/finite_spacing.csv`.

Independent classical `r`-matrix / monodromy comparison:

```bash
python scripts/class5/verify_rmatrix_comparison.py
```

The script checks the classical Poisson algebra, projector recursion, and agreement of the independently generated time-Lax matrix with the quantum-lattice continuum result up to a field-independent scalar matrix.

## Class 6: quantum-to-classical time-Lax limit

Exact symbolic derivation:

```bash
python scripts/class6/verify_coherent_symbolic.py
```

This retains the second spatial coefficient required in Class 6, verifies the shared-site correction including its commutator contribution, checks the Hamiltonian flow and off-shell Lax factorization, and includes a negative control showing that the Class-5-only total-derivative correction fails in Class 6.

Finite-spacing, Yang-Baxter, finite-lattice zero-curvature, and independent off-shell checks:

```bash
python scripts/class6/verify_coherent_numeric.py
```

Outputs are written under `results/class6/`.

## Machine-learning discovery path

Class 5 finite feature-library search:

```bash
python scripts/ml/class5_ml_recovery.py
```

The script samples local spin jets subject only to the unit-spin kinematic constraints, fits the off-shell `F = Q E` residual for an over-complete feature library, and writes `results/ml/class5_reproduced_summary.json`.

Class 6 nilpotent-polynomial search:

```bash
python scripts/ml/class6_ml_recovery.py
```

The script fits the coefficients in the finite matrix library `{I, N, N^2}`, verifies held-out off-shell residuals, and checks the nilpotent matrix-square-root relation. It writes `results/ml/class6_reproduced_summary.json`.

`results/ml/historical_summary.json` records the aggregate errors from the historical runs described in the complete research note. The reconstructed scripts are intentionally separated from these historical metrics.

## Build the documents

Manuscript:

```bash
make paper
```

Complete Japanese research note:

```bash
make note
```

CI rebuilds both PDFs from source and uploads them as workflow artifacts. The repository treats the TeX sources as the durable record.
