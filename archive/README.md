# Archive and provenance

This directory records the provenance of calculations that were originally produced as separate reports during the research process.

The current repository does **not** duplicate every historical PDF/TeX report as a second source of truth. Their derivations were consolidated, without removing the research-history commentary, into:

- `notes/full_research_note_ja.tex` — the self-contained complete calculation/history note;
- `notes/parts/08a_audit_class5_quantum.tex` — independent Class 5 finite-lattice/coherent-symbol audit;
- `notes/parts/08b_audit_class6_quantum*.tex` — independent Class 6 audit;
- `notes/parts/08c_audit_class5_rmatrix*.tex` — independent Class 5 classical r-matrix/monodromy audit.

The executable versions of the calculations are under `scripts/`, with machine-readable outputs under `results/`.

Historical standalone reports that led to these sections included:

1. `class5_lax_coherent_symbol_ja` — Class 5 finite-lattice time-Lax/coherent-symbol derivation;
2. `class6_lax_coherent_symbol_ja` — Class 6 independent derivation and failure boundary for the Class-5-only correction;
3. `class5_sts_rmatrix_comparison_ja` — independent Class 5 classical r-matrix/monodromy comparison.

The complete note is the durable textual record; CI-generated PDFs are build artifacts rather than primary sources.
