# Class 5 Lax correspondence: reproduction and scope

The self-contained derivation is in
`notes/parts/12_class5_generic_lax_comparison.tex`. It is included in the full
Japanese note after the eleven-vertex comparison. A standalone review copy is
built from `notes/class5_lax_comparison_ja.tex`.

From the repository root, with SymPy installed:

```bash
python scripts/class5/verify_lax_correspondence.py
```

The deterministic report is written to
`results/class5/lax_correspondence/verification_results.json`.
Use `--output PATH` to keep a separate rerun. A nonzero residual raises an
exception and gives a nonzero process exit code.

The 22 exact checks cover the reduction of the general Class-5 R matrix and
local Hamiltonian to the manuscript representative, the reconstructed quantum
L operator, reduction to the fundamental R matrix, leading and subleading
differences from the printed expression, the RLL relation in spin 1/2 and spin
1, the field transformation, both Kameyama--Yoshida Lax components, the time
convention, the Hamiltonian flow, and the scalar part of the connection. The
two trace checks additionally verify the scalar coefficients directly from the
manuscript U and V formulas.

Sources: de Leeuw--Fontanella--Nieto Garcia, arXiv:2506.13598v2, equations
(1.1), (1.5), (1.8), (1.11), (1.12), (2.14), (2.23); Kameyama--Yoshida,
arXiv:1405.4467v2, equations (4.1)--(4.9), (A.1), (A.4). The script states the
complete parameter and field conventions and performs no network requests.

These are local algebraic checks on the complexified unit-spin constraint.
The reconstruction is not an author-confirmed erratum. Testing RLL in two
representations is not a proof for arbitrary spin. Local gauge equivalence
neither equates real forms nor fixes periodic monodromy sectors, Poisson
algebras, scattering states, or finite-time quantum dynamics.

The initial local rerun used Python 3.13.5 and SymPy 1.14.0. The dedicated CI
workflow uses Python 3.11 with SymPy 1.14.0 and compares a fresh report with the
committed result. It also builds the full note and the standalone comparison;
review PDFs are artifacts, not working-branch commits.
