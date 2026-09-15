# Project status

Last updated: 2026-09-15.

This file is intended to restore the state of the project without relying on chat history.

## Current physics question

The main problem is the quantum-to-classical correspondence for the **time** component of the Lax pair in the Class 5 and Class 6 spin chains and their Landau--Lifshitz continuum limits.

The spatial continuum Lax matrix follows from the normalized lattice Lax operator by the usual long-wavelength/coherent-state expansion. The time component is subtler because the finite-lattice zero-curvature equation contains products such as

```text
A_{a,n+1} L_{a,n},   L_{a,n} A_{a,n},
```

whose factors share a physical lattice site. The coherent-state lower symbol therefore does not factor into the ordinary product of the separate lower symbols. Under the Landau--Lifshitz continuum scaling this shared-site contribution survives at the same order as the continuum time-Lax equation.

The non-multiplicativity of coherent-state symbols, coherent-state continuum Lax constructions, quantum time-Lax hierarchies, and quantum-to-classical Lax limiting procedures are established ingredients. The project-specific question is the **local** passage from the finite-lattice quantum time equation to the fixed-spin long-wavelength continuum time matrix and the operator products that survive at that order.

## Common spin-1/2 shared-site identity

Use the normalized operator expansions

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

Let `B` denote the order-`epsilon^2` Sutherland residual,

```text
B = [2P, X2 X1 + Y2 + Y1] + [h1, X2 + X1] - Z1 + Z2.
```

For physical spin `1/2`, rank-one product lower symbols, and a first Hamiltonian correction that is scalar on the spin-triplet sector, the local shared-site source obeys

```text
C = 2 i (d_x Xi + [Xi,U]) - i <X1 B>.
```

Hence `B=0` gives

```text
C = 2 i (d_x Xi + [Xi,U]),
Delta V = 2 i Xi + f(lambda) I,
```

with `f(lambda)` independent of the fields and spacetime coordinates.

For Class 5 and Class 6, the triplet-scalar condition is not an extra model assumption: the three Pauli coefficients of the leading spatial operator are linearly independent for `lambda != 0`, and the exchange-even part of the order-two Sutherland relation implies the required total-spin commutant condition. The common rank-three minor is `i/(32 lambda^3)` in both models.

Important interpretation:

- the coincident order-`epsilon^2` connected products cancel;
- their first spatial Taylor coefficients and the coincident order-`epsilon^3` products survive at the continuum time-evolution order;
- the explicit second spatial coefficient `Y` must be retained in intermediate products but cancels from the final common formula;
- in Class 6, `L6^(2)` is a useful model-specific representation of `Xi`, not independent final input;
- for spin `1/2`, `Xi` can be reconstructed from the full spin dependence of `U` by the Pauli algebra;
- the analytic proof allows arbitrary finite auxiliary dimension but uses a two-dimensional physical spin space;
- higher physical spin, higher Hamiltonian flows and finite-time quantum convergence are not established.

Research record and checks:

- `notes/parts/10_shared_site_generalization.tex`
- `scripts/generalization/verify_shared_site_identity.py`
- `scripts/generalization/verify_shared_site_proof.py`
- `results/generalization/shared_site_identity.json`
- `results/generalization/shared_site_proof.json`

The second proof file records 34 exact checks, including Class 5/Class 6/XXZ rank tests and rank-deficient negative controls.

## Class 5 results

- The continuum spatial Lax matrix is obtained from the normalized Class 5 quantum lattice Lax operator.
- The finite-lattice time Lax operator has been taken directly to the continuum while retaining shared-site operator products.
- The Class 5 source reduces to a total spatial derivative.
- The corrected time matrix reproduces the Hamiltonian equations off shell, with curvature and EOM mutually reconstructible for generic spectral parameter.
- The undeformed XXX limit still requires the shared-site correction; the direct time symbol alone has the wrong coefficient in the standard isotropic time matrix.
- An independent classical `r`-matrix / monodromy construction gives the same time matrix up to a field-independent scalar matrix.
- The full Class 5 normalization and Lax comparison has been checked independently. The traceless continuum pair is locally related, after a complex field redefinition and a local auxiliary-space gauge transformation, to the null-like warped `SL(2)` Landau--Lifshitz pair of Kameyama--Yoshida.
- This local relation does **not** identify the real sections, periodic monodromies, scattering data, Poisson structures, or global boundary sectors. The field and gauge transformations are not generically periodic.

Detailed Class 5 comparison:

- `notes/parts/12_class5_generic_lax_comparison.tex`
- `scripts/class5/verify_lax_correspondence.py`
- `results/class5/lax_correspondence/verification_results.json`

The dedicated comparison workflow and symbolic checks are green.

## Class 6 results

- The Class 6 quantum `R`-matrix and normalized lattice Lax operator have been expanded through the order required by the time-Lax limit.
- The second spatial coefficient must be retained in the intermediate shared-site products; the common derivation shows that it is not independent final data.
- The shared-site source contains both a spatial derivative and a nonzero commutator contribution. A spatially uniform spin configuration provides a negative control against a Class-5-only total-derivative correction.
- The corrected time matrix reproduces the Hamiltonian equations and admits an off-shell curvature/EOM factorization with an invertible coefficient map for generic spectral parameter.
- Class 6 is already the eleven-vertex model at the quantum level; no new hierarchy is claimed.
- The continuum Class 6 pair has been matched directly to the known rational eleven-vertex Landau--Lifshitz pair. With

```text
beta^2 = -2 alpha6
z      = beta lambda
k      = -4/beta
tau    = (i beta^2/2) t
```

the traceless spatial matrix, time matrix, zero-curvature convention and spin equation agree, up to auxiliary-space transpose and the field-independent scalar part of `U6`.

Detailed Class 6 comparison:

- `notes/parts/11_class6_eleven_vertex_comparison.tex`
- `scripts/class6/verify_eleven_vertex_comparison.py`
- `results/class6/eleven_vertex_comparison.json`

## Literature positioning and novelty

The current audits are

- `research/shared_site_literature_audit.md`
- `research/novelty_reassessment_2026-09-14.md`
- `paper/literature_review.md`

The manuscript explicitly treats the following as known:

- coherent-state lower-symbol non-multiplicativity;
- coherent-state continuum Lax constructions;
- quantum time-Lax hierarchies;
- quantum-to-classical Lax limiting procedures;
- the classical structures to which the final Class 5 and Class 6 pairs are locally related.

The manuscript does **not** claim a new Class 5 or Class 6 classical hierarchy. The result specific to this work is the direct fixed-spin long-wavelength limit of the **local quantum lattice time-Lax equation**, including the shared-site operator-product terms that survive at the continuum time-evolution order, together with the conditional local identity above.

The literature audit distinguishes this from the global lattice-sum power counting of Avan--Doikou--Sfetsos and from the quantum time-Lax hierarchy of Doikou--Findlay. A targeted search has not found the same local fixed-spin shared-site formula in the checked sources, but literature absence or priority is not established.

## Machine-learning discovery path

Machine learning selected local structures and coefficient relations for subsequent analytic verification.

- Earlier PCM and symmetric-coset tests recovered spectral-parameter families.
- A low-loss `T^{1,1}` candidate was shown analytically to be a fake Lax connection, motivating the off-shell `F = Q E` certificate and rank test.
- Class 5 used twelve temporal structures plus a spatial coefficient, with an off-shell residual and a uniform sparsity penalty.
- Class 6 used eight coefficients in the finite nilpotent matrix library `{I,N,N^2}`; the restored implementation includes an exactly equivalent off-shell objective.
- Final Lax claims are established by analytic factorization, symbolic verification, finite-lattice checks and literature comparison, not by low loss alone.

The recovered original search code is preserved under `archive/ml_original/`. Active scripts and results are under `scripts/ml/` and `results/ml/`. See `REPRODUCIBILITY.md` for exact commands and convention dictionaries.

## Git / manuscript status

The research results have been integrated into `main` in the following order:

1. PR #12: shared-site generalization, research-note integration, XXZ control and Class 6 / eleven-vertex validation -- merged with a normal merge commit.
2. PR #13: Class 5 normalization and Kameyama--Yoshida correspondence -- merged with a normal merge commit.

The remaining manuscript integration is draft PR #14 from `paper/shared-site-framework` to `main`. The paper branch has been synchronized with the merged research results. PR #14 contains the manuscript/editorial changes only:

- common shared-site formula and rank-three argument in Section 2;
- concise Class 5 comparison with the known null-like warped pair;
- concise Class 6 comparison with the known eleven-vertex pair;
- revised Abstract, Introduction, model comparison and Summary;
- updated literature/novelty framing and reproducibility documentation.

The current paper structure is

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

PR #14 CI rebuilds the physics checks, ML checks, manuscript and complete Japanese research note. Merge only after all jobs are green and the generated manuscript PDF has been visually checked.

## Notation decisions

- Quantum lattice pair: `L_{a,n}(u)` and `A_{a,n}(u)`.
- Continuum pair: `U(x,t;lambda)` and `V(x,t;lambda)`.
- `u` is the quantum spectral parameter; `lambda` is the continuum spectral parameter.
- Class 5 lattice/continuum couplings are `kappa_5`, `alpha_5`, with `kappa_5 = epsilon alpha_5`.
- Class 6 lattice/continuum couplings are `kappa_6`, `alpha_6`, with `kappa_6 = epsilon^2 alpha_6`.
- In the common Section 2 derivation, `X`, `Y`, `Z`, `B`, and `Xi` are local generic quantities defined there and are not replacements for the model-labelled coefficients elsewhere.
- Repeated Cartesian indices are summed once this convention has been declared.
- Symbols are introduced only after their physical meaning has been stated.
- Equation-of-motion residuals carry model superscripts `E^(5)` and `E^(6)`.
- Retain the explicit matrix elements of `V_5` and `V_6`; no generator-basis rewrite has been adopted.

## Authoring rules

- No `\\boxed`, `\\fbox`, `\\framebox`, or `\\mbox` in the manuscript or research note.
- Avoid context-dependent labels imported from another paper unless the underlying object has first been defined self-containedly.
- Prefer explicit derivations over short-lived helper notation introduced only to compress algebra.
- Do not present low ML loss as proof of integrability.
- Do not claim coherent-state symbol non-multiplicativity as new.
- Do not claim a new Class 5 or Class 6 classical hierarchy.
- Keep literature input, convention changes, new derivations and novelty assessment logically separate.

## Editorial scope

The analytic quantum-to-classical time-Lax construction is the paper narrative:

```text
known spatial continuum construction + known quantum time hierarchy
        -> local fixed-spin time-equation continuum limit
        -> conditional shared-site identity
        -> Class 5 / Class 6 realizations and independent classical checks.
```

The Introduction gives a short ML overview, and Appendix A supplies the mathematical search formulation and numerical settings. Detailed proof history, controls, normalization audits and global caveats remain in the research record rather than being copied wholesale into the manuscript.
