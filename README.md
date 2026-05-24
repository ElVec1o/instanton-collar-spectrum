# instanton-collar-spectrum

[![DOI](https://zenodo.org/badge/1248611476.svg)](https://doi.org/10.5281/zenodo.20370060)

Asymptotic product structure and collar essential spectrum on SU(2) instanton moduli.

Author: Vico Bonfioli (vicobonfioli@gmail.com)

## Contents

- `paper/CORE/core.tex`, `core.pdf`, 24-page submission paper. Targeted at Annales Henri Poincaré or J. Geom. Phys.
- `paper/EXTRA/extra.tex`, `extra.pdf`, 30-page supplementary companion (arXiv-only).
- `paper/SPLIT.md`, rationale for the CORE/EXTRA split.
- `external-review/REVIEWER_BUNDLE.md`, single-file reviewer bundle combining paper, companion, split rationale, revision log, and all verification artifacts.
- `external-review/REVIEWER_GUIDE.md`, short guide for reviewers.
- `verification/`, computational verification scripts:
  - `lemma_ct_rust/`, independent Rust adaptive Gauss-Kronrod cubature for Lemma 4.1 (Schwinger closed form), machine-precision verified to ~2.5 ulp.
  - Sage / Python scripts for off-diagonal cross-blocks, the Mourre commutator, the SU(N) extension, the Polyakov-Wiegmann alpha, and the M_2 supplementary analyses.

## Main results

- **Theorem 3.2.** lambda_0(M_1(S^4_r)) = 4/r^2, via the Habermann / Doran-Murchadha-Marolf identification of M_1(S^4_r) with hyperbolic 5-space and McKean saturation.
- **Lemma 4.1.** Exact Schwinger closed form for the BPST scale-derivative L^2 cross-term, with explicit next-order coefficient and remainder. Machine-precision verified.
- **Lemma 4.2.** Off-diagonal cross-block bounds (scale-position, position-position, gauge-fixed, bubble-background).
- **Theorem 4.3.** Codim-j Uhlenbeck collar essential bottom = 4j/r^2, via operator-norm metric comparison (1 +/- C eps |log eps|) g_prod with Weyl quasi-modes upper bound and min-max lower bound.
- **Theorem 4.4.** Mourre spectral type with C^{1,1} regularity (Froese-Hislop conjugated radial Schrodinger framework): no singular continuous spectrum above threshold, locally finite point spectrum, absolutely continuous spectrum above 4j/r^2 + C eps |log eps|.
- **Remark 4.5.** SU(N) extension with explicit framing-overlap scalar.

## Status

Final pre-submission revision (revision 6). Six rounds of external blind review have converged on minor-revision acceptance at the AHP / JGP tier.

## Licenses

- Paper, companion, and bundle: CC-BY 4.0 (`LICENSE-paper`). Reuse requires attribution.
- Verification code: MIT (`LICENSE-code`).

## Citation

If you use this work, please cite as:

```
Bonfioli, Vico. Asymptotic product structure and collar essential spectrum
on SU(2) instanton moduli. 2026.
DOI: 10.5281/zenodo.20370060
https://github.com/ElVec1o/instanton-collar-spectrum
```

A `CITATION.cff` file is included for automated tools.
