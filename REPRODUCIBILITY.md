# Reproducibility guide

All commands below are run from the repository root. Install `requirements.txt` for analytic checks, or `requirements-ml.txt` for the complete suite including CPU PyTorch and pandas.

## Class 5: quantum-to-classical time-Lax limit

Exact symbolic derivation:

```bash
python scripts/class5/verify_coherent_symbolic.py
```

This checks the fundamental operator identities, Sutherland relation, coherent-state product correction, corrected time-Lax matrix, Hamiltonian flow, off-shell EOM factorization, and the undeformed XXX limit. Results are written to `results/class5/symbolic_verification.json` and `results/class5/symbolic_expressions.json`.

The Section 5 comparison is checked against the same direct operator contractions: the single-site Pauli product, the common expression for the connected source in terms of the difference of second moments, and the isotropic source and time correction. The Class 6 script below checks the corresponding common expression and isotropic formulas independently.

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

This retains the second spatial coefficient required in Class 6, verifies the shared-site correction including its commutator contribution, checks the Hamiltonian flow and off-shell Lax factorization, and includes a negative control showing that the Class-5-only total-derivative correction fails in Class 6. It also derives the instantaneous spin derivative from the two adjacent quantum bonds and checks its continuum limit against the Hamiltonian flow used in manuscript section 4.3.

Finite-spacing, Yang-Baxter, finite-lattice zero-curvature, and independent off-shell checks:

```bash
python scripts/class6/verify_coherent_numeric.py
```

Outputs are written under `results/class6/`.

## Machine-learning discovery path

The active scripts derive directly from the recovered originals in
`archive/ml_original/`; that directory records their file hashes and origins.
The candidate spaces and numerical protocols are preserved. Install and run:

```bash
pip install -r requirements-ml.txt
make verify-ml
```

`make verify-ml` first runs `scripts/ml/verify_conventions.py`. It compares
the original and adapted samplers, residuals, losses and gradients at arbitrary
coefficients. Independent symbolic formulas also compare the complete recovered
`U` and `V` matrices against the explicit manuscript entries. The checks are
recorded in `results/ml/convention_verification.json`.

### Convention dictionary

Let `y` and `tau` denote the original code's field and time, and let `S`, `t`,
`alpha5`, `alpha6`, `lambda` denote the manuscript quantities. The manuscript
matrices `M5` and `N6` are defined in `scripts/ml/paper_conventions.py`.

| Quantity | Class 5 | Class 6 |
| --- | --- | --- |
| Proper spin-basis rotation | `R5 = diag(-1, 1, -1)` | `R6 = diag(1, -1, -1)` |
| Spin | `S(x,t) = exp(-alpha5*x*M5/2) R5 y(x,-2t)` | `S(x,t) = R6 y(x,-2t)` |
| Old coupling | `c = alpha5/2` | `alpha = alpha6` |
| Old leading spectral coefficient | `z = i/(2*lambda)` | `a0 = i/(2*lambda)` |
| Spatial differentiation | `D_x = d_x + alpha5*M5/2` | `d_x` |
| Constant auxiliary vector rotation | `O5 = exp(alpha5*lambda*M5)` | identity |
| Training residual in paper variables | `-O5^(-1)(F5-Q5 E5)/2` | default `-(F6-A6 E6)/2`; original `-F6/2` with `E6=0` |

The traceless matrix map is `iota(v) = -i v.sigma/2`. The time Lax vector
has the factor `-2` from `tau=-2t`. Fixed scalar matrix parts complete the
representatives used in the manuscript. The Class 5 coefficient map is
`Q5 = O5 [i/(2*lambda) I + beta m5 m5^T]`, where `m5=(i,-1,0)`;
it shares `beta` with the spatial ansatz. Because `O5` is complex orthogonal,
retaining its inverse in the residual is needed to preserve the original
positive numerical norm. The comparison script verifies this for arbitrary
coefficients, including the sparsity contribution and parameter gradients.

The original coefficients are real. The default numerical runs retain that
domain by taking positive imaginary `lambda`. Exact symbolic identities use
general nonzero complex `lambda`.

### Preserved search protocols

| Setting | Class 5 | Class 6 |
| --- | --- | --- |
| Candidate space | 12 temporal vectors plus spatial `beta` | 8 coefficients in `{I,N6,N6^2}` |
| Training loss | off-shell `F=QE`, plus uniform smooth L1 penalty on all 13 coefficients | off-shell `F=A6 E`, equivalent to the original on-shell curvature |
| Training / validation points | 768 / 4096 | 768 / 3072 |
| Adam steps / rate / gradient clip | 600 / 0.03 / 100 | 1300 / 0.02 / 10 |
| L-BFGS max iterations / rate | 50 / 0.8 | 120 / 0.7 |
| L-BFGS gradient / change tolerance | `1e-12` / `1e-14` | `1e-13` / `1e-15` |
| Initial coefficients | independent Gaussian, standard deviation 0.25 | undeformed point plus Gaussian noise, standard deviation 0.08 |
| Independent runs | 3 parameter pairs × 3 seeds | 3 parameter pairs × 2 seeds |

Both L-BFGS optimizers use the strong-Wolfe line search. Class 5 uses penalty
weight `1e-9` during Adam and `5e-11` during L-BFGS, with smoothing `1e-20`
inside each square root. The Class 6 selection of the Adam iterate passed to
L-BFGS is retained from the original. All calculations use double precision.
Exact coefficients are evaluated after training.

The scripts can also be run individually:

- `python scripts/ml/class5_ml_recovery.py` writes coefficients, training history
  and aggregate metrics to `results/ml/class5/`.
- `python scripts/ml/class6_ml_recovery.py` writes coefficients and aggregate
  metrics to `results/ml/class6/`.
- `python scripts/ml/class6_ml_recovery.py --training-mode on-shell --out results/ml/class6_onshell`
  selects the original formulation and writes to a separate output directory.

Class 6 now samples a tangent time derivative independently of the EOM in its
default off-shell mode. The spatial samples are unchanged. Since the spatial
ansatz is `u6=A6 S` with field-independent `A6`, the time-derivative terms cancel
in `F6-A6 E6` for arbitrary coefficients. The loss and its gradients therefore
equal those of the original on-shell objective. This extension adds no new
coefficient constraints or parameters; the convention check compares both
formulations directly against the archived original code.

The recorded rerun used Python 3.12.14 and PyTorch 2.10.0+cpu. Class 5's maximum
physical coefficient error was `3.61e-7`, maximum extra coefficient `8.10e-8`,
and maximum validation MSE `4.94e-14`. Class 6's off-shell maximum coefficient error was
`1.12e-11`, and maximum validation RMS `2.92e-11`. The original on-shell
formulation gave `4.33e-12` and `1.12e-11`, respectively. Numerical optimization may
reach slightly different stopping points across environments; CI checks
recovery tolerances and the exact convention identities.

`results/ml/historical_summary.json` retains the rounded historical metrics
quoted in the research note. The recovered original aggregate files are also
in `archive/ml_original/`. The previous SciPy reconstructions had fewer Class 5
features, omitted its sparsity penalty, and used off-shell training for Class 6.
Their outputs are retained in `archive/ml_reconstructed/`; their code remains
accessible in commit `22ae4b88cb6e96af5dab2484ebf8156faaede0d9`.

The manuscript Introduction summarizes the discovery route. Appendix A gives
the original ansatz spaces expressed in paper variables, the two distinct
losses, sampling, optimization settings and rerun metrics. The final off-shell
certificates and lattice calculations remain in the main text.

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
