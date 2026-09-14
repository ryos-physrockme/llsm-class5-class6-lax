# Project status

Last updated: 2026-09-14.

This file is intended to restore the state of the project without relying on chat history.

## Current physics question

The main problem is the quantum-to-classical correspondence for the **time** component of the Lax pair in the Class 5 and Class 6 spin chains and their Landau--Lifshitz continuum limits.

The spatial continuum Lax matrix follows from the normalized lattice Lax operator by the usual long-wavelength/coherent-state expansion. The time component is subtler because the finite-lattice zero-curvature equation contains products such as

```text
A_{a,n+1} L_{a,n},   L_{a,n} A_{a,n},
```

whose factors share a physical lattice site. The coherent-state lower symbol therefore does not factor into the ordinary product of the separate lower symbols. Under the Landau--Lifshitz continuum scaling this shared-site contribution survives at the same order as the continuum time-Lax equation.

The non-multiplicativity of coherent-state symbols, coherent-state continuum Lax constructions, and quantum time-Lax hierarchies are all established ingredients. The project-specific question is the **local** passage from the finite-lattice time equation to the fixed-spin long-wavelength continuum time matrix and the operator products that survive at that order.

## Common spin-1/2 shared-site identity

For the normalized spatial operator and Hamiltonian expansion

```text
Lhat = I + epsilon X + epsilon^2 Y + ...
h    = 2 P + epsilon h1 + ...
d_u Lhat = epsilon^2 Z + ...
```

define the continuum spatial matrix and the lower-symbol second-moment difference by

```text
U  = X^downarrow
Xi = (X^2)^downarrow - U^2.
```

Let `B` be the order-`epsilon^2` residual of the Sutherland relation,

```text
B = [2P, X2 X1 + Y2 + Y1] + [h1, X2 + X1] - Z1 + Z2.
```

For physical spin `1/2`, rank-one product-state lower symbols, and the sufficient condition that the first Hamiltonian correction is scalar on the spin-triplet sector, the local shared-site source obeys

```text
C = 2 i (d_x Xi + [Xi,U]) - i <X1 B>.
```

Thus an integrable lattice Lax operator satisfying the Sutherland relation at this order (`B=0`) gives

```text
C = 2 i (d_x Xi + [Xi,U]),
Delta V = 2 i Xi + f(lambda) I,
```

where the scalar function does not affect the curvature.

Important interpretation:

- the coincident order-`epsilon^2` connected products cancel;
- their first spatial Taylor coefficients and the coincident order-`epsilon^3` products survive at the continuum time-evolution order;
- the explicit second spatial coefficient `Y` must be retained in the intermediate products but cancels from the final common formula;
- in Class 6, `L6^(2)` is therefore a useful model-specific representation of `Xi`, not independent final input;
- for spin `1/2`, `Xi` can be reconstructed from the complete spin dependence of `U` by the Pauli algebra;
- the analytic proof does not fix the finite auxiliary-space dimension, but it does use a two-dimensional physical spin space;
- higher physical spin, higher Hamiltonian flows and finite-time quantum convergence are not established.

The detailed derivation is in `notes/parts/10_shared_site_generalization.tex`. Independent checks are in `scripts/generalization/verify_shared_site_identity.py` and `scripts/generalization/verify_shared_site_proof.py`.

## Established Class 5 results

- The Class 5 continuum spatial Lax matrix is recovered from the normalized quantum lattice Lax operator.
- The finite-lattice time Lax operator has been taken directly to the continuum while retaining shared-site operator products.
- The connected contribution reduces to the Class-5-specific total-derivative form proportional to `d_x(U_5^2)`.
- The corrected continuum time Lax matrix reproduces the Class 5 Hamiltonian equations of motion off shell; curvature and EOM are mutually reconstructible for generic spectral parameter.
- The undeformed XXX / isotropic Landau--Lifshitz limit provides a control and still requires the shared-site correction. The direct time symbol alone has the wrong coefficient in the standard isotropic time matrix.
- An independent classical `r`-matrix / monodromy construction gives the same time-Lax matrix up to a field-independent scalar matrix.
- Ref. de Leeuw--Fontanella--Nieto García constructs the Class 5 spatial Lax matrix but explicitly leaves the compatible time/companion matrix open; the finite-lattice derivation fills this concrete gap.

Primary verification code: `scripts/class5/`.

A separate branch `research/class5-generic-lax-comparison` has been reserved for additional Class 5 work in another session. Changes from that branch must be reconciled with the manuscript branch before merging overlapping paper sources.

## Established Class 6 results

- The Class 6 quantum `R`-matrix and normalized lattice Lax operator have been expanded through the order required by the time-Lax limit.
- The second spatial coefficient must be retained in the intermediate shared-site products. In the model-specific identity `(L6^(1))^2 = I/(4 lambda^2) + 2 L6^(2)`, it makes the second moment `Xi_6` and its commutator transparent; the common derivation shows that it is not independent final data.
- The shared-site contribution contains both a spatial derivative and a non-vanishing commutator term. A spatially uniform spin configuration gives an explicit negative control showing that the Class-5-only total-derivative correction fails.
- The corrected Class 6 continuum time Lax matrix reproduces the Hamiltonian equations of motion and admits an off-shell Lax/EOM factorization with invertible coefficient map for generic spectral parameter.
- Class 6 is already known at the quantum level as the eleven-vertex model; the project does not claim a new hierarchy.
- The manuscript Class 6 pair has now been compared directly with the known rational eleven-vertex Landau--Lifshitz pair. With
  `beta^2 = -2 alpha6`, `z = beta lambda`, `k = -4/beta`, and
  `tau = (i beta^2/2) t`, the traceless spatial matrix, the time matrix,
  the zero-curvature convention, and the spin equation agree exactly, up to
  the auxiliary-space transpose and the field-independent scalar part of `U6`.

Primary verification code: `scripts/class6/`, including `verify_eleven_vertex_comparison.py`.
The detailed convention dictionary is recorded in `notes/parts/11_class6_eleven_vertex_comparison.tex`.

## Literature positioning of the common identity

The current literature audit is in `research/shared_site_literature_audit.md`.

- Avan--Doikou--Sfetsos (2010) already develops coherent-state continuum Lax constructions, discusses noncommutativity of expectation values with nonlinear operations, and tracks coincident/overlapping site indices in global lattice sums. Their global suppression argument uses the reduction in the number of independent site sums.
- Doikou--Findlay (2017/2020) already constructs quantum time-Lax hierarchies and discusses coherent-state dynamics.
- Therefore the manuscript must not claim that coherent-state nonmultiplicativity, overlapping sites in general, or the existence of quantum time-Lax operators is new.
- The distinction relevant here is that the **local** lattice zero-curvature equation has no global-sum suppression. Under `t = epsilon^2 t_lat`, the first spatial Taylor term and coincident cubic product survive precisely at the continuum time-evolution order.
- A targeted literature audit has not found the same local fixed-spin formula `C = 2 i (d_x Xi + [Xi,U])` in the checked primary sources, but literature absence or priority is not established.

## Machine-learning discovery path

Machine learning selected local structures and coefficient relations for subsequent analytic verification.

- Earlier PCM and symmetric-coset tests recovered spectral-parameter families.
- A low-loss `T^{1,1}` candidate was shown analytically to be a fake Lax connection, motivating the off-shell `F = Q E` certificate and rank test.
- Class 5 used twelve temporal structures plus a spatial coefficient, with an off-shell residual and a uniform sparsity penalty. Nine historical runs selected the physical coefficients and suppressed eight extra terms.
- Class 6 used eight coefficients in the finite nilpotent matrix library `{I, N, N^2}`, with on-shell curvature training. Six historical runs exposed the finite matrix-square-root coefficient pattern. The restored implementation also supports an exactly equivalent off-shell `F-A6 E` objective, now the default; the time-derivative terms cancel for arbitrary fitted coefficients.
- Final claims are established by analytic factorization, symbolic verification, finite-lattice checks, and literature comparison.

Historical aggregate metrics are in `results/ml/historical_summary.json`. The original training scripts have been recovered from the user attachments and Library files; matching hashes and unchanged sources are preserved in `archive/ml_original/`. The active scripts derive directly from them, using the manuscript spin, time, deformation and spectral conventions. They retain the original candidate spaces, spatial sampling, initialization, sparsity and Adam/L-BFGS settings. `scripts/ml/verify_conventions.py` checks losses and gradients at arbitrary coefficients, the on/off-shell equivalence for Class 6, and exact agreement with the manuscript Lax matrices. Current numerical results are in `results/ml/class5/` and `results/ml/class6/`, with the original on-shell Class 6 rerun in `results/ml/class6_onshell/`. Earlier least-squares reconstruction outputs have been moved to `archive/ml_reconstructed/`.

## Manuscript status

`main` still contains the previously audited Class 5/Class 6 manuscript. Two draft branches are currently open for the new physics and framing:

1. `research/shared-site-generalization`, draft PR #12: analytic common identity, research-note integration, literature audit, XXZ control and Class 6 eleven-vertex comparison. The integrated CI (physics, ML, paper build and full research-note build) is green.
2. `paper/shared-site-framework`, draft PR #14: manuscript restructuring stacked on the research result. It adds the common spin-`1/2` formula to Section 2, sharpens the Introduction, adds a concise Class 6 eleven-vertex comparison, updates Section 5 and the Summary, and reframes the Abstract.

The paper-branch structure is currently:

```text
1 Introduction
2 Quantum lattice setup and continuum limit
  2.x Common form of the shared-site correction
3 Class 5 model
4 Class 6 model
  4.4 Comparison with the known eleven-vertex Landau--Lifshitz pair
5 Comparison of Class 5 and Class 6
6 Summary and discussion
A Machine-learning searches for the Lax pairs
```

The common Section 2 formula is intentionally shorter than the complete research-note proof: the rank-three derivation of the triplet condition and the rank-deficient controls remain in the research note. XXZ remains a research/control example rather than a third full manuscript model.

The first complete CI build of the paper restructuring succeeded and its PDF was visually inspected. The new common subsection occupies roughly one and a half pages, the Class 6 known-pair comparison roughly one page, and no page overflow was observed. Subsequent edits only clarified first-use notation in the Abstract, Introduction and common subsection; current CI should be checked before merge.

## Notation decisions

- Quantum lattice pair: `L_{a,n}(u)` and `A_{a,n}(u)`.
- Continuum pair: `U(x,t;lambda)` and `V(x,t;lambda)`.
- `u` is the quantum spectral parameter; `lambda` is the continuum Lax spectral parameter when that convention is used.
- In manuscript section 3, `kappa_5` and `alpha_5` denote the Class 5 lattice and continuum deformation parameters, with `kappa_5 = epsilon alpha_5`.
- In manuscript section 4, `kappa_6` and `alpha_6` denote the Class 6 lattice and continuum deformation parameters, with `kappa_6 = epsilon^2 alpha_6`.
- `K_6` is the four-dimensional two-site deformation operator, and `L_6^(j)` denotes the coefficient of `epsilon^j` in the normalized lattice spatial operator.
- In the common Section 2 derivation, `X`, `Y`, `Z`, `B`, and `Xi` are local generic quantities defined there and are not replacements for the model-labelled coefficients elsewhere.
- Repeated Cartesian indices are summed once this convention has been declared.
- The Class 5 two-site deformation operator is `K_5`. Model subscripts are retained on `A^(j)` and the overlapping-product differences in both model sections.
- Symbols are introduced only after their physical meaning has been stated.
- Equation-of-motion residuals carry model superscripts, `E^(5)` and `E^(6)`, matching the curvature notation.
- Retain the explicit matrix elements of `V_5` and `V_6`; a generator-basis rewrite has not been adopted.
- Bond lower symbols retain their lattice indices, with the continuum density evaluated at `(x_n,t)` for model `m = 5,6`.

## Authoring rules

- No `\boxed`, `\fbox`, `\framebox`, or `\mbox` in the manuscript or research note.
- Avoid context-dependent labels imported from another paper unless the underlying object has first been defined self-containedly.
- Prefer explicit derivations over introducing short-lived helper symbols solely to compress algebra.
- Do not present a low ML loss as a proof of integrability.
- Do not claim that coherent-state symbol non-multiplicativity is new.
- Do not claim a new Class 6 / eleven-vertex hierarchy.
- Keep literature input, convention changes, new derivations, and novelty assessment logically separate.

## Editorial scope

The analytic quantum-to-classical time-Lax construction remains the main paper narrative. The stronger current framing is:

```text
known continuum/spatial construction + known quantum time hierarchy
        -> local fixed-spin time-equation continuum limit
        -> conditional shared-site identity
        -> Class 5 / Class 6 realizations and independent classical checks.
```

The Introduction gives an ML overview, and Appendix A supplies the mathematical search formulation and numerical settings. Code execution instructions belong in the repository guides. The detailed proof history, rank-condition analysis, controls and literature audit remain in the research record rather than being copied wholesale into the paper.
