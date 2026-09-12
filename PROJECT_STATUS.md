# Project status

Last updated: 2026-09-12.

This file is intended to restore the state of the project without relying on chat history.

## Current physics question

The main problem is the quantum-to-classical correspondence for the **time** component of the Lax pair in the Class 5 and Class 6 spin chains and their Landau-Lifshitz continuum limits.

The spatial continuum Lax matrix follows from the normalized lattice Lax operator by the usual long-wavelength/coherent-state expansion. The time component is subtler because the finite-lattice zero-curvature equation contains products such as

```text
A_{a,n+1} L_{a,n},   L_{a,n} A_{a,n},
```

whose factors share a physical lattice site. The coherent-state lower symbol therefore does not factor into the ordinary product of the separate lower symbols. Under the Landau-Lifshitz continuum scaling this shared-site contribution survives at the same order as the continuum time-Lax equation.

The non-multiplicativity of coherent-state symbols is standard. The project-specific result is the explicit finite-lattice-to-continuum evaluation of this contribution for Class 5 and Class 6 and its absorption into a local correction of the continuum time-Lax matrix.

## Established Class 5 results

- The Class 5 continuum spatial Lax matrix is recovered from the normalized quantum lattice Lax operator.
- The finite-lattice time Lax operator has been taken directly to the continuum while retaining shared-site operator products.
- The connected contribution reduces to the Class-5-specific total-derivative form proportional to `∂_x(U_5^2)`.
- The corrected continuum time Lax matrix reproduces the Class 5 Hamiltonian equations of motion off shell; curvature and EOM are mutually reconstructible for generic spectral parameter.
- The undeformed XXX / isotropic Landau-Lifshitz limit provides a control and still requires the shared-site correction.
- An independent classical `r`-matrix / monodromy construction gives the same time-Lax matrix up to a field-independent scalar matrix.

Primary verification code: `scripts/class5/`.

## Established Class 6 results

- The Class 6 quantum `R`-matrix and normalized lattice Lax operator have been expanded through the order required by the time-Lax limit.
- Unlike Class 5, the second spatial coefficient is essential.
- The shared-site contribution contains both a spatial derivative and a non-vanishing commutator term. A spatially uniform spin configuration gives an explicit negative control showing that the Class-5-only total-derivative correction fails.
- The corrected Class 6 continuum time Lax matrix reproduces the Hamiltonian equations of motion and admits an off-shell Lax/EOM factorization with invertible coefficient map for generic spectral parameter.
- Class 6 is already known at the quantum level as the eleven-vertex model; the project does not claim a new hierarchy. The result is an explicit bridge/dictionary and the nontrivial time-Lax continuum construction.

Primary verification code: `scripts/class6/`.

## Machine-learning discovery path

Machine learning is a discovery/audit tool rather than the final proof.

- Earlier PCM and symmetric-coset tests recovered spectral-parameter families.
- A low-loss `T^{1,1}` candidate was shown analytically to be a fake Lax connection, motivating the off-shell `F = Q E` certificate and rank test.
- Class 5 used a symmetry-informed over-complete local feature library. Nine historical runs selected the correct physical coefficients and suppressed distractors.
- Class 6 used the finite nilpotent matrix library `{I, N, N^2}`. Six historical runs exposed the finite matrix-square-root coefficient pattern.
- Final claims are established by analytic factorization, symbolic verification, finite-lattice checks, and literature comparison.

Historical aggregate metrics are in `results/ml/historical_summary.json`. The executable scripts in `scripts/ml/` are clean reconstructions from the documented equations, because the original historical scripts were not available as raw files when this repository was initialized.

## Manuscript status

`paper/main.tex` is the current draft.

Written in the current draft:

1. common quantum-lattice setup and continuum limit;
2. Class 5 quantum model and spatial Lax matrix;
3. Class 5 finite-lattice derivation of the time Lax matrix;
4. Class 5 zero-curvature/EOM equivalence;
5. Class 5 independent classical `r`-matrix comparison.
6. Class 6 section 4.1: quantum Hamiltonian and `R`-matrix, continuum
   scaling, exact spatial-operator expansion, and continuum spatial Lax matrix.
7. Class 6 section 4.2: finite-lattice time-operator expansion, overlapping
   lower-symbol products, local correction involving the second spatial
   coefficient, and the corrected time Lax matrix with explicit entries.

Next technical section to write: **Section 4.3, Class 6 zero-curvature/EOM
equivalence**, including the independent Hamiltonian flow and the off-shell
factorization. After that, write the Class 5/Class 6
comparison, scope, Introduction, and Conclusions.

Current planned paper structure:

```text
1 Introduction
2 Quantum lattice setup and continuum limit
3 Class 5 model
4 Class 6 model
5 Comparison of Class 5 and Class 6
6 Conclusions
Appendices
```

## Notation decisions

- Quantum lattice pair: `L_{a,n}(u)` and `A_{a,n}(u)`.
- Continuum pair: `U(x,t;λ)` and `V(x,t;λ)`.
- `u` is the quantum spectral parameter; `λ` is the continuum Lax spectral parameter when that convention is used.
- In manuscript section 3, `κ_5` and `α_5` denote the Class 5 lattice and
  continuum deformation parameters, with `κ_5 = ε α_5`. They correspond to
  `a` and `alpha` in the Class 5 verification scripts. The added model
  subscripts are a notation change only; the normalizations are unchanged.
- In manuscript section 4, `κ_6` and `α_6` denote the Class 6 lattice and
  continuum deformation parameters, with `κ_6 = ε² α_6`. They correspond to
  `a` and `alpha` in the existing Class 6 verification scripts and research
  note; this is a notation change only. `K_6` is the four-dimensional
  two-site deformation operator, and `L_6^(j)` denotes the coefficient of
  `ε^j` in the normalized lattice spatial operator.
- Repeated Cartesian indices are summed once this convention has been declared.
- Symbols are introduced only after their physical meaning has been stated.

## Authoring rules

- No `\boxed`, `\fbox`, `\framebox`, or `\mbox` in the manuscript or research note.
- Avoid context-dependent labels imported from another paper unless the underlying object has first been defined self-containedly.
- Prefer explicit derivations over introducing short-lived helper symbols solely to compress algebra.
- Do not present a low ML loss as a proof of integrability.
- Do not claim that coherent-state symbol non-multiplicativity is new.
- Keep literature input, convention changes, new derivations, and novelty assessment logically separate.

## Open editorial decision

The physics calculation is currently strong enough to remain the main paper narrative. The ML part should at minimum be mentioned as the discovery path and preserved in the research note/repository. After the Class 6 paper section is written, decide whether the ML details belong in an appendix/supplement or deserve a more central role.
