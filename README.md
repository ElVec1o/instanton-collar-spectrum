# instanton-collar-spectrum

[![DOI](https://zenodo.org/badge/1248611476.svg)](https://doi.org/10.5281/zenodo.20370060)

LaTeX sources, PDFs, and verification code for

> V. Bonfioli, *Asymptotic product structure and collar essential spectrum on SU(2) instanton moduli*, 2026. 25 pages. [`paper/CORE/core.pdf`](paper/CORE/core.pdf)

and the supplementary companion (arXiv-only, not part of the submission)

> *AHS cohomology and the Yang-Mills mass gap: supplementary observations*, 2026. 30 pages. [`paper/EXTRA/extra.pdf`](paper/EXTRA/extra.pdf)

## Results

M_k(S^4_r) is the charge-k SU(2) instanton moduli space over the round 4-sphere of radius r, with the L^2 metric.

- **Theorem 3.2.** lambda_0(M_1(S^4_r)) = 4/r^2. Habermann / Doran-Murchadha-Marolf identification of M_1 with hyperbolic 5-space of curvature -1/r^2; McKean's bound is saturated.
- **Lemma 4.1.** The BPST scale-derivative L^2 cross-term in elementary closed form: 24 pi^2 (cosh d sinh d - d)/sinh^3 d, where d is the hyperbolic distance between the bubble points in the M_1 identification. Limits: 16 pi^2 at d -> 0 (the diagonal norm), 48 pi^2 s_1 s_2 / R^2 asymptotic at d -> infinity.
- **Lemma 4.2.** Off-diagonal cross-block bounds for the L^2 metric on the multi-instanton collar.
- **Theorem 4.3.** On the codimension-j Uhlenbeck collar, the bottom of the essential spectrum equals 4j/r^2, with two-sided rate C eps |log eps| / r^2 (Corollary 4.5). Operator-norm metric comparison (1 +/- C eps |log eps|) g_prod; Weyl quasi-modes; min-max.
- **Theorem 4.8.** Mourre estimate in the Froese-Hislop conjugated radial framework with C^{1,1} regularity: no singular continuous spectrum above threshold, locally finite point spectrum, absolutely continuous spectrum above 4j/r^2 + C eps |log eps|.

The global bottom lambda_0(M_k(S^4_r)) for k >= 2 is open; the paper proves the collar-localized statement and says so. The companion contains conditional and exploratory material.

## Checking the paper

Suggested order for a first pass:

1. Lemma 4.1 (self-contained; one Schwinger computation). Its closed form is verified independently twice: `verification/lemma_ct_rust/` (4-D integral vs 1-D form, 106 points, max rel. error ~2.5 ulp; table in `RESULTS.md`) and `verification/lemma_ct_elementary.py` (elementary/hyperbolic form, symbolic identities plus 50-digit quadrature).
2. Theorem 4.3 (uses Lemmas 4.1, 4.2; the metric comparison is eq. 3.3, the appendix carries Lemma 4.2).
3. Theorem 4.8 (Mourre; the conjugated-frame commutator identities are in `verification/mourre_estimate_cusp.sage`, which also shows why the Froese-Hislop conjugated frame is required, and `verification/mourre_iterated_commutator.sage` carries the C^{1,1} input).

## Verification

| Artifact | Checks | Run |
|---|---|---|
| `verification/lemma_ct_rust/` | Lemma 4.1: 4-D integral vs closed form, 106-point grid | `cargo run --release` |
| `verification/lemma_ct_elementary.py` | Lemma 4.1: elementary/hyperbolic form, symbolic + 50-digit grid | `python3 lemma_ct_elementary.py` |
| `verification/lemma_ct_montecarlo.py` | Lemma 4.1: Monte Carlo | `python3 lemma_ct_montecarlo.py` |
| `verification/lemma_ct_higher_order.sage` | Lemma 4.1: next-order coefficient | `sage lemma_ct_higher_order.sage` |
| `verification/lemma_od_position_terms.sage` | Lemma 4.2: scale-position terms | `sage lemma_od_position_terms.sage` |
| `verification/mourre_estimate_cusp.sage` | Theorem 4.8: commutator identities | `sage mourre_estimate_cusp.sage` |
| `verification/mourre_iterated_commutator.sage` | Lemmas A.7, A.8 | `sage mourre_iterated_commutator.sage` |
| `verification/lemma_ct_SUN.sage` | Remark 4.10: SU(N) | `sage lemma_ct_SUN.sage` |
| `verification/alpha_*.sage`, `m2_*` | companion only | as above |

## Building

```
cd paper/CORE  && pdflatex core.tex  && pdflatex core.tex
cd paper/EXTRA && pdflatex extra.tex && pdflatex extra.tex
```

Expected: CORE 25 pages, EXTRA 30 pages, no undefined references.

## Provenance

The author is an independent researcher. The manuscript was produced with extensive LLM assistance (drafting, internal adversarial review, verification tooling); the author takes responsibility for the content. It has not been reviewed by a human expert. The numerical verification is machine-checked and independent of the derivations. Corrections and counterexamples: vicobonfioli@gmail.com.

## License and citation

Paper and companion: CC-BY 4.0 (`LICENSE-paper`). Code: MIT (`LICENSE-code`).

```
Bonfioli, Vico. Asymptotic product structure and collar essential spectrum
on SU(2) instanton moduli. 2026. DOI: 10.5281/zenodo.20370060
https://github.com/ElVec1o/instanton-collar-spectrum
```
