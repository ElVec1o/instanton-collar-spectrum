# CORE / EXTRA split

The paper is split into a tightly-scoped journal submission (`CORE/`) and an
arXiv-only supplementary companion (`EXTRA/`). The split is the result of a
deliberate "include in CORE only what is 100% proven; demote everything else to
the companion" decision made during drafting.

## `CORE/`, the journal submission

**`CORE/core.tex`**, *Asymptotic product structure and collar essential spectrum on SU(2) instanton moduli.* 25 pages.

Three theorems, two lemmas, one self-contained weighted-Sobolev appendix. All content is fully proved (or, for Theorem 3.2, a clean folklore-saturation packaging with explicit radius bookkeeping). No "proof sketch" labels; no citation-chain-dependent statements.

Contents:

- §3 Theorem (M_1 McKean saturation): λ_0(M_1(S^4_r)) = 4/r²
- §4.1 Lemma (Schwinger cross-term closed form): exact elementary closed form for the BPST scale-derivative L²-pairing, equal to 24π²(cosh d sinh d − d)/sinh³d with d the hyperbolic distance in the M₁ ≅ H⁵ identification; machine-precision verified ~2.5 ulp via independent Rust adaptive Gauss–Kronrod cubature
- §4.2 Lemma (off-diagonal cross-block bounds): scale-position closed form, position-position IR-divergence diagnosis, gauge-fixed cross-block C s₁ s₂/R²(1+|log|), bubble-background s_i²/R
- §4.3 Theorem 4.3 (codim-j collar essential bottom = 4j/r²) + Corollary 4.5 (two-sided rate Cε|log ε|/r²): rigorous via operator-norm metric comparison (1 ± Cε|log ε|) g_prod, Weyl quasi-modes upper bound, min-max lower bound
- §4.4 Theorem 4.8 (Mourre spectral type with C¹·¹ regularity): no σ_sc above threshold, locally finite σ_pp, AC spectrum on (4j/r² + Cε|log ε|, ∞)
- Remark 4.10 (SU(N) extension), end of §4.4
- §5 Remarks and open questions
- Appendix A (~6 pages): self-contained Bartnik β = −3/2 weighted-Sobolev proof of Lemma 4.2(iv)-(v); scale-uniform Coulomb estimate via conformal rescaling; gauge-equivalence + Bianchi for the tail cancellation; three-region cross-block split with explicit constant tracking; bubble-background s²/R via parity + Sobolev; Δ_g vs Δ_prod bookkeeping with Hardy on hyperbolic cusp; iterated-commutator bounds (Lemmas A.7, A.8) closing the Mourre C¹·¹ regularity gap.

## `EXTRA/`, the arXiv-only supplementary companion

**`EXTRA/extra.tex`**, *AHS cohomology and the Yang-Mills mass gap: supplementary observations.* 30 pages.

A heterogeneous collection of observations and conditional structural results that did not meet the fully-proven bar required for CORE:

- An observation that no Yang-Mills mass-gap formula can depend on the AHS Euler characteristic alone (bundle-dependence + perturbative independence).
- A McKean-Singer ζ-identity on the AHS complex, constraining Singer-style regularization schemes.
- A KKN-curved-surface Hessian formula with explicit Polyakov-Wiegmann curvature-mass coefficient α = (N²−1)/9 (derived via three independent routes; convention-dependent).
- The global λ_0(M_2(S^4_r)) question as an **open question** with sharp diagnostic identifying the orbit-volume function W(ρ) at ρ → ∞ as the deeper analytical obstruction.
- Two **conditional propositions** explicitly labeled as such, with "proof outline":
    - SO(5)-isotypic reduction of M_2 to the invariant subspace (conditional on the V_true cusp asymptote)
    - Mazzeo-Melrose finiteness of σ_pp(0, 4/r²) + meromorphic resolvent (conditional on stratified iterated-edge calculus adapted to Uhlenbeck)
- A T⁴ vacuum-structure discussion (Witten twisted sectors) plus a structural Bakry-Émery lattice-continuum bridge.

The companion is *not* part of the journal submission. It is intended for arXiv deposit alongside CORE.

## Why the split

The criterion:

1. CORE contains only fully proved content.
2. Material that rests on additional inputs (numerical Bargmann robustness; stratified iterated-edge 0-calculus citation chains) is placed in the companion and labeled as conditional.

## Versioned history

`paper/old/` contains every numbered draft v4–v15 plus the v16 omnibus that was split into CORE + EXTRA. These are kept for the record only; none are current.
