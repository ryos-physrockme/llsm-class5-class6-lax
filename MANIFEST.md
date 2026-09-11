# Repository manifest and recovery map

Use this file when recovering project context in a new session.

## Read first

1. `PROJECT_STATUS.md` — current physics conclusions, paper status, notation and style decisions.
2. `README.md` — repository purpose and quick commands.
3. `REPRODUCIBILITY.md` — verification scope and reproducibility notes.
4. `AGENTS.md` — working rules for future automated/editor sessions.

## Scientific record

- `notes/full_research_note_ja.tex` — complete self-contained research/history note.
- `notes/parts/02_ml_search.tex` — ML discovery history, off-shell certificate, Class 5/6 coefficient searches.
- `notes/parts/08a_audit_class5_quantum.tex` — independent Class 5 quantum/coherent-symbol derivation.
- `notes/parts/08b_audit_class6_quantum*.tex` — independent Class 6 derivation and Class-5-only negative control.
- `notes/parts/08c_audit_class5_rmatrix*.tex` — independent classical r-matrix/monodromy closure for Class 5.

## Manuscript

- `paper/main.tex` — modular paper root.
- `paper/sections/02*.tex` — common quantum-lattice/coherent-state formulation.
- `paper/sections/03*.tex` — completed Class 5 section through classical r-matrix comparison.
- `paper/sections/04_class6.tex` — next section to write.

## Executable checks

- `scripts/class5/verify_coherent_symbolic.py`
- `scripts/class5/verify_coherent_numeric.py`
- `scripts/class5/verify_rmatrix_comparison.py`
- `scripts/class6/verify_coherent_symbolic.py`
- `scripts/class6/verify_coherent_numeric.py`

## ML

- `scripts/ml/class5_ml_recovery.py`
- `scripts/ml/class6_ml_recovery.py`

The ML scripts above are **clean reconstructions from the equations documented in the complete research note**. They are not claimed to be byte-for-byte copies of the original historical training scripts.

- `results/ml/historical_summary.json` records aggregate numbers from the historical runs.
- `results/ml/class5_reproduced_summary.json` and `class6_reproduced_summary.json` record the current reconstructed runs.

Keep those provenance classes separate.

## Verification commands

```bash
python -m pip install -r requirements.txt
make verify
make paper
make note
```

GitHub Actions runs the physics checks, ML reconstructions, TeX style constraints, and both LaTeX builds.

## Current notation

- lattice spatial Lax operator: `L_{a,n}(u)`
- lattice time Lax operator: `A_{a,n}(u)` in the manuscript
- continuum Lax matrices: `U(x,t;lambda)` and `V(x,t;lambda)`
- lattice spacing: `epsilon`

Some independent-audit text embedded from earlier research reports uses the older calligraphic time-operator notation. Treat that as archival notation, not the current manuscript convention.
