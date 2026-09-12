# Superseded least-squares reconstructions

These are the outputs of the earlier SciPy reconstructions, preserved for
provenance. Their source is available at commit
`22ae4b88cb6e96af5dab2484ebf8156faaede0d9` in `scripts/ml/`.

Class 5 used eight temporal features, omitted the uniform sparsity penalty,
and used least squares. Class 6 also used least squares and an off-shell
residual. These implementations differed from the historical training scripts.

The active code now derives directly from the recovered original PyTorch
scripts. Class 5 restores twelve temporal features, thirteen learned real
parameters, sparsity, and Adam/L-BFGS. Class 6 retains the original eight
parameters and optimizer, with an exactly equivalent off-shell extension.
See `archive/ml_original/` and `REPRODUCIBILITY.md`.
