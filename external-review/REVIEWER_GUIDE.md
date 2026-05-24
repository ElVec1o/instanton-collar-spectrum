# Reviewer's guide

This bundle contains everything needed to evaluate the paper

> **"Asymptotic product structure and collar essential spectrum on SU(2) instanton moduli"** by V. Bonfioli

as well as an arXiv-only supplementary companion.

## What's in the bundle

```
external-review/
├── REVIEWER_GUIDE.md          (this file)
├── paper/
│   ├── core.pdf               <- PRIMARY: 20-page submission paper
│   ├── core.tex               <- LaTeX source of the primary paper
│   ├── extra.pdf              <- SUPPLEMENTARY: 30-page arXiv-only companion
│   ├── extra.tex              <- LaTeX source of the companion
│   └── SPLIT.md               <- one-page rationale for primary vs supplementary
└── verification/              <- computational artifacts cited in the paper
    ├── lemma_ct_rust/         <- Rust adaptive Gauss–Kronrod (Lemma 4.1)
    │   ├── Cargo.toml
    │   ├── src/main.rs
    │   └── RESULTS.md         <- machine-precision grid table (~2.5 ulp)
    ├── lemma_5_2_schwinger.py     <- Python Monte Carlo (Lemma 4.1 sanity)
    ├── lemma_ct_higher_order.sage <- Sage next-order Schwinger expansion
    ├── lemma_od_position_terms.sage <- Sage scale-position + IR-divergence diagnosis
    ├── lemma_ct_SUN.sage          <- Sage SU(N) extension (Remark 4.5)
    ├── mourre_estimate_cusp.sage  <- Symbolic H^5 commutator identity (Theorem 4.4)
    ├── mourre_iterated_commutator.sage <- Iterated commutator (Lemmas A.7, A.8)
    ├── alpha_first_principles.sage <- Sage alpha = (N^2-1)/9 (companion only)
    ├── alpha_PW_S2.sage           <- Sage Polyakov-Wiegmann check (companion only)
    ├── m2_bargmann_true_V.sage    <- Honest M_2 robustness diagnostic (companion)
    ├── m2_closure_higher_order.sage <- Closed-form Phi(mu) (companion)
    ├── m2_closure_fast.py         <- Numerical V_true(s) evaluation (companion)
    └── m2_isotypic_decomposition.sage <- Conditional SO(5) proposition (companion)
```

## Where to start

1. **Read `paper/core.pdf`** (20 pages). This is the only document being submitted to the journal. It contains three theorems and two lemmas, plus a six-page self-contained appendix.

2. (Optional) **Skim `paper/extra.pdf`** for the supplementary structural observations, explicitly labeled as conditional or arXiv-only. The two structural propositions there are *not* claims; they are recorded with explicit "proof outline" labelling and honest caveats.

3. (Optional) **Spot-check the verification.** The most consequential artifact is `verification/lemma_ct_rust/RESULTS.md`, which records an independent Rust adaptive Gauss–Kronrod cubature against the original 4-dimensional integral on a 106-point grid spanning four orders of magnitude in `s/R`, with maximum relative error ≈ 2.5 ulp (≈ 5.5 × 10⁻¹⁶, machine precision in double). If you want to rebuild it, the `lemma_ct_rust/` directory is a standalone Cargo project (`cargo run --release` after `cd verification/lemma_ct_rust/`).

## The primary paper's claims at a glance

| Result | Status | Verification artifact |
|---|---|---|
| Theorem 3.2: λ_0(M_1(S^4_r)) = 4/r² | Standard via Habermann/DMM + McKean (folklore packaging) | (none needed) |
| Lemma 4.1: exact Schwinger closed form for the BPST scale-derivative cross-term | New | `lemma_ct_rust/` (machine precision); `lemma_5_2_schwinger.py` (MC); `lemma_ct_higher_order.sage` (next-order coefficient and log) |
| Lemma 4.2: off-diagonal cross-block bounds | New (proved in Appendix A) | `lemma_od_position_terms.sage` |
| Theorem 4.3: codimension-j Uhlenbeck collar essential bottom = 4j/r² | New, quantitative upgrade of Taubes/DK folklore | (chain through Lemma 4.1 + 4.2) |
| Theorem 4.4: Mourre spectral type (no σ_sc above threshold, locally finite σ_pp, AC spectrum) | New, journal-grade | `mourre_estimate_cusp.sage`; `mourre_iterated_commutator.sage` |
| Remark 4.5: SU(N) extension | Real but trivial (recorded as Remark) | `lemma_ct_SUN.sage` |

## What is explicitly open

The global L²-spectral bottom λ_0(M_k(S^4_r)) for k ≥ 2 is open. The paper proves only the collar-localized essential-spectrum statement and disclaims the global bottom explicitly (Remark "Domain monotonicity bounds the global bottom" and the corresponding Open Question in §5). Two structural reductions of the k=2 case are recorded in the supplementary companion (extra.pdf) as *conditional propositions* rather than theorems.

## Author's self-assessment of the package

The author's view is that the primary paper is journal-acceptable at the tier of Annales Henri Poincaré, Journal of Geometry and Physics, Annals of Global Analysis and Geometry, Letters in Mathematical Physics, or SIGMA. The supplementary companion is intended for arXiv only.

A reviewer is welcome to verify, dispute, or recalibrate this assessment.
