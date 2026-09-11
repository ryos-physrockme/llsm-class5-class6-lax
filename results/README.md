# Results and provenance

`results/class5/` and `results/class6/` contain representative machine-readable outputs of the analytic and finite-spacing verification scripts.

`results/ml/` distinguishes two kinds of records:

- `historical_summary.json`: aggregate metrics copied from the historical ML runs documented in the complete research note. The original raw historical training scripts were not available when this repository was initialized.
- `class5_reproduced_summary.json` and `class6_reproduced_summary.json`: outputs of the clean reconstructed searches in `scripts/ml/`, implemented from the documented feature libraries and off-shell equations.

CI does not rely on committed result files to decide success. It reruns the executable checks from source.
