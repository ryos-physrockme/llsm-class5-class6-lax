# Manuscript consistency audit

Review baseline: `f96342319be848892f8c3eff0e1203cbdfaf5821`, PR #11.

The review covers the complete manuscript, with particular attention to the
overlapping-product derivations, the classical monodromy comparison, and
Appendix A. It compares the paper with the analytic verification scripts,
the active ML implementation and recorded results, the detailed research
note, and the primary references used for the convention and monodromy
formulas.

Primary-source comparisons use
[de Leeuw, Fontanella and Nieto García](https://arxiv.org/html/2506.13598v2)
for the Class 5 Hamiltonian and time convention, and
[Doikou and Karaiskos](https://arxiv.org/pdf/1105.5042) for the generating
time matrix and local monodromy expansion. The supplied
Krippendorf–Lüst–Syvaeri paper provides the optimization context.

## Findings and revisions

| Location | Finding | Revision |
| --- | --- | --- |
| Section 2.2 | The scaling sentence called the time scaling model dependent, although both models use the same time convention. | State the common time scaling and refer to the model sections for the coupling scalings. |
| Section 3.3 | The Hamiltonian variation was less explicit than in Section 4.3. The comparison with the reference used arrows that left the direction of the time change ambiguous. | State the unconstrained Cartesian variation with periodic boundary conditions. Define the reference field, coupling and time and give their explicit relations. |
| Section 3.4 | The passage from the exact normalized monodromy to a rank-one projector omitted the formal asymptotic qualification. | Explain the dominant exponential branch and the omission of the other branch beyond all orders of the local power series. Display the normalized monodromy whose local expansion is used. |
| Section 3.4 | The generating-flow identity was implicit, and the Hamiltonian flow was identified before the time-matrix comparison. | Display the Poisson bracket with the logarithm of the transfer matrix and the two equations fixing the first projector coefficient. Identify the Hamiltonian vector field after the two time matrices agree up to a scalar. |
| Appendix A.1 | The Class 5 scalar matrix components were prescribed without explaining their origin or their role in the fit. | Relate them to the local conservation law for the transverse spin component. State that the numerical fit determines the traceless sector. |
| Appendix A.1–A.2 | The input time derivative was listed, but its exact cancellation in both residuals was unstated. | Explain the cancellation for arbitrary coefficients and what remains to be fitted. This also makes the scope of the off-shell numerical objective explicit. |
| Appendix A.3 | The statement that analytic coefficients enter only after optimization was too broad: the normalization and the undeformed Class 6 initial point are supplied. Finite parameter samples were not clearly distinguished from general spectral dependence. | Separate supplied values, coefficient fitting at each fixed parameter pair, and analytic verification for general nonzero spectral parameter. |
| Abstract, Introduction, Summary | The numerical contribution needed the same scope throughout the paper. | Describe recovery at tested parameter values and the separate analytic verification. |
| Notation and prose | The ML auxiliary representation and complex orthogonal matrix could be stated more precisely; one equation ended with a comma before a new sentence. | Make the definitions explicit and correct the punctuation. |

## Mathematical checks

The analytic checks cover the quantum operator identities, both continuum
time matrices, Hamiltonian equations, curvature/EOM factorization and its
inverse, the common isotropic limit, and the independent Class 5 classical
Poisson algebra and projector recursion. The reference convention change is
also checked directly against the transformed Class 5 vector equation.

The overlapping-product calculations retain the required order of
operations: expand the quantum coefficients, form the products on the
shared physical site, take their lower symbols, and expand the neighbouring
spin arguments. The coincident quadratic terms cancel; their first spatial
Taylor terms and the coincident cubic terms survive. The Class 6 second
spatial coefficient contributes to both the derivative and commutator
parts. Its spatially uniform example tests the latter independently.

The ML checks compare losses and gradients at arbitrary coefficients with
the archived searches, and compare the complete recovered matrices with the
main-text entries. The audit adds explicit checks of the Class 5 time-data
cancellation and its parameter gradients, the scalar conservation identity,
and invertibility of both coefficient maps before fitting. Existing checks
cover the Class 6 time-data cancellation and gradient equivalence.

Local verification passed 34 Class 5 symbolic checks, 58 Class 6 symbolic
checks, four independent classical r-matrix checks, and 31 ML convention
and fit-structure checks. The Class 6 uniform-spin negative control rejected
the incomplete total-derivative correction as expected.

## Scope of the conclusions

- The monodromy calculation concerns the formal local hierarchy. It does
  not assert convergence of the power series or control of exponentially
  small terms.
- The sampled losses test coefficients in the specified finite candidate
  spaces. General spectral dependence and curvature/EOM equivalence are
  established analytically. The equations of motion and the candidate
  tensor structures are inputs to the searches.
- The lattice-to-continuum claim concerns instantaneous quantum
  commutators evaluated on product coherent states. The discussion already
  identifies finite-time correlation and error estimates as further work.
- The correction and the two model examples are not presented as a theorem
  for arbitrary quantum chains or spin representations.

The PR CI supplies the complete verification run and the review PDF for the
revised source. Build provenance is included in the review handoff.
