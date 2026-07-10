# instanton-collar-spectrum

[![DOI](https://zenodo.org/badge/1248611476.svg)](https://doi.org/10.5281/zenodo.20370060)

LaTeX sources, PDFs, and verification code for

> V. Bonfioli, *Asymptotic product structure and collar essential spectrum on SU(2) instanton moduli*, 2026. 25 pages.

and the supplementary companion

> *AHS cohomology and the Yang-Mills mass gap: supplementary observations*, 2026. 30 pages.

## Results

Let M_k(S^4_r) denote the charge-k SU(2) instanton moduli space over the round 4-sphere of radius r, with the L^2 metric.

- **Theorem 3.2.** lambda_0(M_1(S^4_r)) = 4/r^2. (Habermann / Doran-Murchadha-Marolf identification of M_1 with hyperbolic 5-space of curvature -1/r^2; McKean's bound is saturated.)
- **Lemma 4.1.** The BPST scale-derivative L^2 cross-term in elementary closed form: equal to 24 pi^2 (cosh d sinh d - d)/sinh^3 d, where d is the hyperbolic distance between the two bubble points in the upper-half-space model of H^5 (the M_1 identification). The d -> 0 limit recovers the diagonal norm 16 pi^2; the d -> infinity limit gives the 48 pi^2 s_1 s_2 / R^2 asymptotic with next-order coefficient and remainder bound.
- **Lemma 4.2.** Off-diagonal cross-block bounds for the L^2 metric on the multi-instanton collar: scale-position asymptotics, position-position IR behavior, gauge-fixed bound, bubble-background bound.
- **Theorem 4.3.** On the codimension-j Uhlenbeck collar, the bottom of the essential spectrum equals 4j/r^2. (Operator-norm metric comparison (1 +/- C eps |log eps|) g_prod; Weyl quasi-modes; min-max.)
- **Theorem 4.4.** Mourre estimate in the Froese-Hislop conjugated radial framework with C^{1,1} regularity: no singular continuous spectrum above the threshold, locally finite point spectrum, absolutely continuous spectrum above 4j/r^2 + C eps |log eps|.
- **Remark 4.5.** SU(N) extension.

The global bottom lambda_0(M_k(S^4_r)) for k >= 2 is open; the paper proves the collar-localized statement only. The companion contains conditional and exploratory material and is not part of the submission.

## Layout

| Path | Contents |
|---|---|
| `paper/CORE/` | Submission paper (`core.tex`, `core.pdf`) |
| `paper/EXTRA/` | Supplementary companion (`extra.tex`, `extra.pdf`) |
| `paper/SPLIT.md` | Criterion for the CORE/EXTRA split |
| `external-review/` | Single-file bundle (paper + companion + code, one markdown file), reviewer guide, revision history |
| `verification/` | Verification code (Rust, SageMath, Python) |

## Verification

| Artifact | Checks | Run |
|---|---|---|
| `verification/lemma_ct_rust/` | Lemma 4.1 closed form against the original 4-D integral; 106-point grid over four orders of magnitude in s/R; max relative error ~2.5 ulp (`RESULTS.md`) | `cargo run --release` |
| `lemma_ct_elementary.py` | Lemma 4.1 elementary/hyperbolic closed form: symbolic identities plus 50-digit quadrature grid | `python3 lemma_ct_elementary.py` |
| `lemma_5_2_schwinger.py` | Lemma 4.1, Monte Carlo | `python3 lemma_5_2_schwinger.py` |
| `lemma_ct_higher_order.sage` | Lemma 4.1 next-order coefficient | `sage lemma_ct_higher_order.sage` |
| `lemma_od_position_terms.sage` | Lemma 4.2 scale-position terms | `sage lemma_od_position_terms.sage` |
| `mourre_estimate_cusp.sage` | Theorem 4.4 commutator identities (unconjugated failure and conjugated-frame identity) | `sage mourre_estimate_cusp.sage` |
| `mourre_iterated_commutator.sage` | Lemmas A.7, A.8 | `sage mourre_iterated_commutator.sage` |
| `lemma_ct_SUN.sage` | Remark 4.5 | `sage lemma_ct_SUN.sage` |
| `alpha_*.sage`, `m2_*.sage`, `m2_closure_fast.py` | Companion-only material | as above |

## Building

```
cd paper/CORE  && pdflatex core.tex  && pdflatex core.tex
cd paper/EXTRA && pdflatex extra.tex && pdflatex extra.tex
```

Expected: CORE 25 pages, EXTRA 30 pages, no undefined references.

## Provenance

The author is an independent researcher. The manuscript was produced with extensive LLM assistance (drafting, internal adversarial review, verification tooling); the author takes responsibility for the content. It has not been reviewed by a human expert. The numerical verification above is machine-checked and independent of the derivations. Corrections and counterexamples: vicobonfioli@gmail.com.

## License and citation

Paper and companion: CC-BY 4.0 (`LICENSE-paper`). Code: MIT (`LICENSE-code`).

```
Bonfioli, Vico. Asymptotic product structure and collar essential spectrum
on SU(2) instanton moduli. 2026. DOI: 10.5281/zenodo.20370060
https://github.com/ElVec1o/instanton-collar-spectrum
```

`CITATION.cff` is included.
