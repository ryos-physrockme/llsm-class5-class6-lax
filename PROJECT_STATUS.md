# Project status

Last updated: 2026-09-15.

This file restores the state of the project without relying on chat history.

## Current physics result

The project studies the quantum-to-classical correspondence for the **time** component of the Lax pair in the Class 5 and Class 6 spin chains and their Landau--Lifshitz continuum limits.

The spatial continuum Lax matrix follows from the normalized lattice Lax operator by a standard long-wavelength/coherent-state expansion. The local lattice time equation is subtler because

```text
A_{a,n+1} L_{a,n},   L_{a,n} A_{a,n}
```

contain factors sharing a physical site. Their coherent-state lower symbols do not factor into products of separate lower symbols. In the fixed-spin Landau--Lifshitz scaling the surviving connected contribution occurs at the same order as the continuum time evolution.

Coherent-state nonmultiplicativity, continuum Lax constructions, quantum time-Lax hierarchies, and quantum-to-classical Lax limits are treated as known. The result specific to this work is the direct continuum limit of the **local quantum lattice time-Lax equation** and the organization of its shared-site contribution.

## Common spin-1/2 identity

Use

```text
Lhat = I + epsilon X + epsilon^2 Y + ...
h    = 2 P + epsilon h1 + ...
d_u Lhat = epsilon^2 Z + ...
```

and define

```text
U  = X^downarrow
Xi = (X^2)^downarrow - U^2.
```

Let

```text
B = [2P, X2 X1 + Y2 + Y1] + [h1, X2 + X1] - Z1 + Z2
```

be the order-`epsilon^2` Sutherland residual. For physical spin `1/2`, rank-one product lower symbols, and a first Hamiltonian correction scalar on the triplet sector,

```text
C = 2 i (d_x Xi + [Xi,U]) - i <X1 B>.
```

For the lattice models considered here, `B=0`, so

```text
C = 2 i (d_x Xi + [Xi,U]),
Delta V = 2 i Xi + f(lambda) I.
```

For both Class 5 and Class 6, the triplet condition follows from the exchange-even part of the Sutherland relation because the three Pauli coefficients of the leading spatial operator are linearly independent for `lambda != 0`; the relevant rank-three minor is `i/(32 lambda^3)`.

The explicit second spatial coefficient `Y` must be kept in intermediate products but cancels from the final common formula. The analytic proof allows arbitrary finite auxiliary dimension but uses a two-dimensional physical spin space. Higher physical spin, higher Hamiltonian flows and finite-time quantum convergence are not established.

Primary research record:

- `notes/parts/10_shared_site_generalization.tex`
- `scripts/generalization/verify_shared_site_identity.py`
- `scripts/generalization/verify_shared_site_proof.py`
- `results/generalization/shared_site_identity.json`
- `results/generalization/shared_site_proof.json`

## Class 5

- The continuum spatial matrix is obtained from the normalized Class 5 quantum Lax operator.
- The finite-lattice time operator is taken directly to the continuum while retaining shared-site products.
- The Class 5 source reduces to a total spatial derivative.
- The corrected time matrix is off-shell equivalent to the Hamiltonian equations for generic spectral parameter.
- The undeformed XXX limit still requires the shared-site correction.
- An independent classical `r`-matrix/monodromy construction gives the same time matrix up to a field-independent scalar matrix.
- The traceless continuum pair is locally related, after a complex field redefinition and local auxiliary-space gauge transformation, to the null-like warped `SL(2)` Landau--Lifshitz pair of Kameyama--Yoshida.
- This local relation does not identify real sections, periodic monodromies, scattering data, Poisson structures, or global boundary sectors; the transformations are not generically periodic.

Detailed comparison:

- `notes/parts/12_class5_generic_lax_comparison.tex`
- `scripts/class5/verify_lax_correspondence.py`
- `results/class5/lax_correspondence/verification_results.json`

## Class 6

- The Class 6 quantum `R`-matrix and normalized lattice Lax operator are expanded through the required continuum order.
- The second spatial coefficient is required in intermediate products but is not independent final data in the common formula.
- The source contains a spatial derivative and a nonzero commutator contribution; a spatially uniform spin gives a negative control against a total-derivative-only correction.
- The corrected time matrix is off-shell equivalent to the Hamiltonian equations for generic spectral parameter.
- Class 6 is already the eleven-vertex model; no new hierarchy is claimed.
- With

```text
beta^2 = -2 alpha6
z      = beta lambda
k      = -4/beta
tau    = (i beta^2/2) t
```

the traceless spatial matrix, time matrix, zero-curvature convention and spin equation agree with the known rational eleven-vertex Landau--Lifshitz structure, up to auxiliary-space transpose and the field-independent scalar part of `U6`.

Detailed comparison:

- `notes/parts/11_class6_eleven_vertex_comparison.tex`
- `scripts/class6/verify_eleven_vertex_comparison.py`
- `results/class6/eleven_vertex_comparison.json`

## Novelty and literature scope

Audits:

- `research/shared_site_literature_audit.md`
- `research/novelty_reassessment_2026-09-14.md`
- `paper/literature_review.md`

The manuscript does **not** claim a new Class 5 or Class 6 classical hierarchy. It treats coherent-state lower-symbol nonmultiplicativity, coherent-state continuum Lax constructions, quantum time-Lax hierarchies, and quantum-to-classical Lax limiting procedures as established.

The main distinction from the global lattice-sum arguments in Avan--Doikou--Sfetsos is that the local zero-curvature equation has no suppression by the number of independent site sums. Under `t = epsilon^2 t_lat`, the first spatial Taylor term of the quadratic overlap and the coincident cubic term survive at the dynamical order. A targeted literature audit has not found the same local fixed-spin formula in the checked primary sources, but literature absence or priority is not claimed.

## Machine-learning role

The ML searches identify coefficient patterns within stated finite ansatz spaces; low loss is not used as proof of integrability. Final claims are established by analytic operator identities, off-shell curvature/EOM factorization, finite-lattice checks and literature comparisons. The recovered original search code is under `archive/ml_original/`; active scripts and results are under `scripts/ml/` and `results/ml/`. See `REPRODUCIBILITY.md`.

## Manuscript and Git status

The full research and manuscript revision has been integrated into `main` in this order:

1. PR #12 -- shared-site generalization, research-note integration, XXZ control and Class 6/eleven-vertex validation -- normal merge.
2. PR #13 -- Class 5 normalization and Kameyama--Yoshida correspondence -- normal merge.
3. PR #14 -- manuscript restructuring and final novelty framing -- squash merge.

The final manuscript source is `paper/main.tex`. Before PR #14 was merged, CI passed all physics checks, ML checks, manuscript compilation and complete Japanese research-note compilation. The resulting 36-page manuscript PDF was rendered page-by-page and visually reviewed.

Current paper structure:

```text
1 Introduction
2 Quantum lattice setup and continuum limit
  2.x Common form of the shared-site correction
3 Class 5 model
  3.x Relation to the null-like warped SL(2) Landau--Lifshitz pair
4 Class 6 model
  4.x Comparison with the known eleven-vertex Landau--Lifshitz pair
5 Comparison of Class 5 and Class 6
6 Summary and discussion
A Machine-learning searches for the Lax pairs
```

XXZ remains a research/control example rather than a third full manuscript model.

## Notation and authoring rules

- Quantum lattice pair: `L_{a,n}(u)`, `A_{a,n}(u)`.
- Continuum pair: `U(x,t;lambda)`, `V(x,t;lambda)`.
- `u` is the quantum spectral parameter; `lambda` is the continuum spectral parameter.
- `kappa_5 = epsilon alpha_5`; `kappa_6 = epsilon^2 alpha_6`.
- The generic Section 2 quantities `X`, `Y`, `Z`, `B`, `Xi` are defined there and do not replace model-labelled coefficients elsewhere.
- Symbols are introduced only after their physical meaning is stated.
- Do not present low ML loss as proof of integrability.
- Do not claim coherent-state symbol nonmultiplicativity as new.
- Do not claim a new Class 5 or Class 6 classical hierarchy.
- Keep literature input, convention changes, new derivations and novelty assessment logically separate.

The paper narrative is

```text
known spatial continuum construction + known quantum time hierarchy
        -> local fixed-spin time-equation continuum limit
        -> conditional shared-site identity
        -> Class 5 / Class 6 realizations and independent classical checks.
```
