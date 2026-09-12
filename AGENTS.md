# Repository instructions for future assistants and coding agents

Read `PROJECT_STATUS.md`, `README.md`, and `REPRODUCIBILITY.md` before changing physics formulas or manuscript text.

## Source of truth

- `paper/main.tex`: current compressed manuscript draft.
- `notes/full_research_note_ja.tex`: self-contained long-form calculation/history note.
- `scripts/`: executable verification and original-derived ML searches in the manuscript conventions.
- `results/`: current machine-readable outputs.
- `archive/`: independent reports kept for provenance, not the current manuscript source.

## Required workflow for physics changes

1. Identify the analytic claim being changed.
2. Update or add an executable check in `scripts/` when possible.
3. Run the relevant script(s), preferably `make verify` before merging substantial changes.
4. Update the research note before compressing the result into the paper if the change alters the calculation history or conventions.
5. Use the pull-request CI build to compile the affected LaTeX document and inspect the resulting PDF.
6. Deliver the review PDFs as described below. After merging, confirm that CI has synchronized the tracked PDFs on `main` with the committed sources.

## PDF review and publication

- For manuscript or note changes, obtain the affected PDFs from the successful pull-request CI run's `compiled-documents` artifact. Check that the run corresponds to the current source revision; the PDF tracked on a working branch may be stale.
- Before requesting review in ChatGPT, save these PDFs to the active session's Library and provide direct file links. Include the source commit SHA and the workflow run URL in the handoff so the reviewed version is identifiable. A CI link alone is not a PDF handoff.
- Keep generated PDF changes out of working-branch commits. The existing `publish-pdfs` job updates the tracked PDFs only after successful verification and compilation on a push to `main`.
- Do not merge a manuscript pull request solely to make its review PDF available.
- After merging, check that `publish-pdfs` succeeded and that the tracked PDF contains the expected manuscript changes; successful LaTeX compilation alone does not confirm publication.

## Writing and notation constraints

- Every nonstandard symbol must be defined before use and, when helpful, introduced as `<physical concept> <symbol>`.
- Do not rely on parameter names or labels that make sense only inside a cited paper.
- No `\boxed`, `\fbox`, `\framebox`, or `\mbox` in `paper/` or `notes/`.
- Avoid unnecessary `\subsubsection` proliferation; use paragraph-level logical flow unless a real section boundary exists.
- Avoid short-lived helper variables that make the derivation read like source code.
- Use standard literature terminology; do not invent technical terms.
- Quantum lattice time Lax operator: `A_{a,n}`. Continuum time Lax matrix: `V`.
- Quantum lattice spatial Lax operator: `L_{a,n}`. Continuum spatial Lax matrix: `U`.
- In the manuscript, use `κ_5, α_5` for the Class 5 lattice and continuum deformation parameters, and `κ_6, α_6` for Class 6. Their scalings are `κ_5 = ε α_5` and `κ_6 = ε² α_6`.

## Scientific-claim constraints

- ML is a discovery/audit tool, not the final proof.
- Genuine Lax claims require off-shell curvature/EOM equivalence, not only on-shell flatness.
- The non-multiplicativity of coherent-state lower symbols is standard; novelty claims concern the explicit finite-lattice time-Lax continuum bridge and model-specific surviving correction.
- Class 5 is related to a known Jordanian/null-like Landau-Lifshitz hierarchy.
- Class 6 is quantum-mechanically the known eleven-vertex model; do not claim a new Class 6 hierarchy.
- Distinguish local identities from claims requiring boundary conditions.

## CI expectations

The main branch should pass `.github/workflows/ci.yml`, including physics verification, original-code convention checks and ML searches, manuscript compilation, and complete-note compilation.

The original ML sources are preserved in `archive/ml_original/`. Preserve their candidate spaces, sampling, losses, and optimizer settings when changing conventions: Class 5 uses twelve temporal features plus one spatial coefficient, a uniform sparsity penalty, and off-shell training; Class 6 originally used eight coefficients and on-shell training. The current Class 6 default uses the user-approved, exactly equivalent off-shell `F-A6 E` objective; `--training-mode on-shell` retains the original formulation. Both final analytic certificates are off-shell. See `REPRODUCIBILITY.md` for the exact dictionary and commands. The paper describes mathematical setup, numerical settings and results; keep execution commands and file organization in the repository guides.
