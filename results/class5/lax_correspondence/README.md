# Class 5 Lax correspondence: reproduction and scope

The self-contained derivation is in
`notes/parts/12_class5_generic_lax_comparison.tex`. It is included in the full
Japanese note after the eleven-vertex comparison. A standalone review copy is
built from `notes/class5_lax_comparison_ja.tex`.

From the repository root, with SymPy installed:

```bash
python scripts/class5/verify_lax_correspondence.py
python scripts/class5/verify_hamiltonian_correspondence.py
```

The deterministic reports are written to
`results/class5/lax_correspondence/verification_results.json` and
`results/class5/lax_correspondence/hamiltonian_verification.json`.
Use `--output PATH` to keep a separate rerun. A nonzero residual raises an
exception and gives a nonzero process exit code.

The first script contains 22 exact checks covering the reduction of the
general Class-5 R matrix and local Hamiltonian to the manuscript representative,
the reconstructed quantum L operator, reduction to the fundamental R matrix,
leading and subleading differences from the printed expression, the RLL
relation in spin 1/2 and spin 1, the field transformation, both
Kameyama--Yoshida Lax components, the time convention, the Hamiltonian flow,
and the scalar part of the connection.

The second script verifies independently that the spatially dependent field
rotation removes the derivative interaction from the Class 5 Hamiltonian
density, that the resulting density maps to the null-like warped
Landau--Lifshitz Hamiltonian density with the stated parameter normalization,
and that the same field and time dictionary maps the equations of motion.
These checks establish the dynamical comparison before the Lax-pair gauge
transformation is considered.

Sources: de Leeuw--Fontanella--Nieto Garcia, arXiv:2506.13598v2, equations
(1.1), (1.5), (1.8), (1.11), (1.12), (2.14), (2.23); Kameyama--Yoshida,
arXiv:1405.4467v2, equations (4.1)--(4.9), (A.1), (A.4). The scripts state the
parameter and field conventions and perform no network requests.

These are local algebraic checks on the complexified unit-spin constraint.
The reconstruction is not an author-confirmed erratum. Testing RLL in two
representations is not a proof for arbitrary spin. The field and gauge
transformations are not generically periodic, so the comparison does not
identify real forms or global periodic sectors.

The initial local rerun used Python 3.13.5 and SymPy 1.14.0. The dedicated CI
workflow uses Python 3.11 with SymPy 1.14.0 and compares fresh reports with the
committed results. It also builds the full note and the standalone comparison;
review PDFs are artifacts, not working-branch commits.
