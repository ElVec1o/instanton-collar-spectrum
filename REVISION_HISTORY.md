# Revision history

Substantive changes to the manuscript, in chronological order. All internal review was LLM-based; no human expert has reviewed this work. The changes below are to the mathematics and its presentation, and stand on their own.

1. **Split.** The original omnibus draft was split into CORE (only fully proved content; the journal submission) and EXTRA (conditional and exploratory material; arXiv-only companion).

2. **Independent verification of Lemma 4.1.** A standalone Rust adaptive Gauss-Kronrod cubature was added, computing the original 4-dimensional integral and the 1-dimensional Schwinger closed form independently on a 106-point grid spanning four orders of magnitude in s/R. Maximum relative discrepancy: about 2.5 ulp. See `verification/lemma_ct_rust/RESULTS.md`.

3. **Mourre regularity closed.** Appendix A was extended with iterated dilation-derivative and iterated-commutator bounds (Lemmas A.7, A.8) establishing the C^{1,1} hypothesis needed for the Sahbani-type Mourre conclusion, rather than assuming it.

4. **Metric comparison and Green's-kernel estimates made explicit.** The operator-norm comparison between the moduli metric and the product metric on the collar, and the horizontal-projection Green's-kernel estimate, were written out with explicit constants in Appendix A instead of being asserted.

5. **Commutator identity corrected (conjugated frame).** The identity [-Delta_{H^5_r}, iA] = 2(-Delta - 4/r^2), as stated in an earlier draft with A the radial dilation generator, is false as an operator identity on H^5_r: direct computation leaves a transverse term 2i r^{-2} y^2 Delta_x. The Mourre argument was rewritten in the Froese-Hislop framework: spherical-harmonic (respectively Fourier-mode) decomposition, reduction of each sector to a radial Schrodinger operator, conjugation by y^{(n-1)/2} and the substitution u = log y to the model H_conj = -d^2/du^2 + 4/r^2 on L^2(R), Mourre estimate there with the dilation generator, transfer back by direct sum. Both computations are recorded in `verification/mourre_estimate_cusp.sage`.

6. **Revision 6 (final consolidations).**
   - The Mourre residual-factor argument is stated at the identity level: A_j acts trivially on the residual factor M_{k-j}, so [-Delta_{M_{k-j}} otimes 1, iA_j] = 0, and the full commutator lives on the cusp factors.
   - Two uncited bibliography entries removed.
   - The Lemma 4.2(v) / Proposition A.10 constant is restated as the scale-invariant Sobolev embedding H^2(B_R) -> L^infty(B_R) in dimension 4, with constant C R^{-2}.
   - The horosphere transverse spectrum is treated by Fourier-mode decomposition on R^4 (continuous transverse spectrum on [0, infinity)), replacing an incorrect "discrete after fibering" phrasing; the zero mode sets the 4/r^2 threshold and positive modes shift the bottom up.
   - Theorem 4.3 states lim rather than liminf, since the upper and lower bounds match.

7. **Elementary and hyperbolic closed form for Lemma 4.1; rate corollary (v1.2).** The 1-dimensional Schwinger integral of Lemma 4.1 is evaluated in elementary terms: the cross-term equals 24 pi^2 (cosh d sinh d - d)/sinh^3 d, where d is the hyperbolic distance between the two bubble points in the upper-half-space model of H^5 (the M_1 identification of Proposition 3.1). The d -> 0 limit reproduces the diagonal norm 16 pi^2 continuously; the d -> infinity limit recovers the 48 pi^2 s_1 s_2 / R^2 asymptotic; the logarithm in the (s/R)^4 remainder is identified as the distance d itself. Verified symbolically and to 50 significant digits (`verification/lemma_ct_elementary.py`). A quantitative two-sided rate corollary, |lambda_0^ess(U_eps^(j)) - 4j/r^2| <= C eps |log eps| / r^2, is added after Theorem 4.3. CORE is now 25 pages.
