# Repository instructions for future assistants and coding agents

Read `PROJECT_STATUS.md`, `README.md`, and `REPRODUCIBILITY.md` before changing physics formulas or manuscript text.

## Source of truth

- `paper/main.tex`: current compressed manuscript draft.
- `notes/full_research_note_ja.tex`: self-contained long-form calculation/history note.
- `scripts/`: executable verification and reconstructed ML searches.
- `results/`: current machine-readable outputs.
- `archive/`: independent reports kept for provenance, not the current manuscript source.

## Required workflow for physics changes

1. Identify the analytic claim being changed.
2. Update or add an executable check in `scripts/` when possible.
3. Run the relevant script(s), preferably `make verify` before merging substantial changes.
4. Update the research note before compressing the result into the paper if the change alters the calculation history or conventions.
5. Compile the affected LaTeX document and inspect the resulting PDF.
6. Keep generated PDFs synchronized with the committed source for material manuscript/note changes when PDFs are intentionally tracked; CI always builds artifact PDFs from source.

## Writing and notation constraints

- Every nonstandard symbol must be defined before use and, when helpful, introduced as `<physical concept> <symbol>`.
- Do not rely on parameter names or labels that make sense only inside a cited paper.
- No `\boxed`, `\fbox`, `\framebox`, or `\mbox` in `paper/` or `notes/`.
- Avoid unnecessary `\subsubsection` proliferation; use paragraph-level logical flow unless a real section boundary exists.
- Avoid short-lived helper variables that make the derivation read like source code.
- Use standard literature terminology; do not invent technical terms.
- Quantum lattice time Lax operator: `A_{a,n}`. Continuum time Lax matrix: `V`.
- Quantum lattice spatial Lax operator: `L_{a,n}`. Continuum spatial Lax matrix: `U`.

## Scientific-claim constraints

- ML is a discovery/audit tool, not the final proof.
- Genuine Lax claims require off-shell curvature/EOM equivalence, not only on-shell flatness.
- The non-multiplicativity of coherent-state lower symbols is standard; novelty claims concern the explicit finite-lattice time-Lax continuum bridge and model-specific surviving correction.
- Class 5 is related to a known Jordanian/null-like Landau-Lifshitz hierarchy.
- Class 6 is quantum-mechanically the known eleven-vertex model; do not claim a new Class 6 hierarchy.
- Distinguish local identities from claims requiring boundary conditions.

## CI expectations

The main branch should pass `.github/workflows/ci.yml`, including physics verification, reconstructed ML searches, manuscript compilation, and complete-note compilation.
