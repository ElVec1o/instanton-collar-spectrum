# Reviewer's guide

Contents and verification map for

> **"Asymptotic product structure and collar essential spectrum on SU(2) instanton moduli"** by V. Bonfioli

and its supplementary companion.

## Files

| Path | Description |
|---|---|
| `paper/CORE/core.pdf` (`core.tex`) | Submission paper, 24 pages |
| `paper/EXTRA/extra.pdf` (`extra.tex`) | Supplementary companion, 30 pages, not part of the submission |
| `paper/SPLIT.md` | Rationale for the CORE/EXTRA split |
| `verification/lemma_ct_rust/` | Rust adaptive Gauss-Kronrod cubature (Lemma 4.1); `RESULTS.md` holds the 106-point grid table |
| `verification/lemma_5_2_schwinger.py` | Python Monte Carlo check of Lemma 4.1 |
| `verification/lemma_ct_higher_order.sage` | Next-order Schwinger expansion |
| `verification/lemma_od_position_terms.sage` | Scale-position terms and IR-divergence diagnosis (Lemma 4.2) |
| `verification/lemma_ct_SUN.sage` | SU(N) extension (Remark 4.5) |
| `verification/mourre_estimate_cusp.sage` | Symbolic H^5 commutator identity (Theorem 4.4) |
| `verification/mourre_iterated_commutator.sage` | Iterated commutators (Lemmas A.7, A.8) |
| `verification/alpha_first_principles.sage`, `alpha_PW_S2.sage` | alpha = (N^2-1)/9 (companion only) |
| `verification/m2_*.sage`, `m2_closure_fast.py` | M_2 diagnostics (companion only) |

## Reading order

1. `paper/CORE/core.pdf` (24 pages): three theorems, two lemmas, one self-contained appendix. This is the only document intended for journal submission.
2. Optional: `paper/EXTRA/extra.pdf` for the supplementary material. Its two structural propositions are labeled "proof outline" and are conditional, not claimed as theorems.
3. Optional: the verification. The main artifact is `verification/lemma_ct_rust/RESULTS.md`: an independent Rust adaptive Gauss-Kronrod cubature of the original 4-dimensional integral against the 1-dimensional closed form of Lemma 4.1, on a 106-point grid spanning four orders of magnitude in s/R; maximum relative error is about 2.5 ulp (about 5.5e-16). To rebuild: `cd verification/lemma_ct_rust && cargo run --release`.

## Claims and their verification artifacts

| Result | Status | Artifact |
|---|---|---|
| Theorem 3.2: lambda_0(M_1(S^4_r)) = 4/r^2 | Follows from Habermann/DMM + McKean; stated with explicit radius bookkeeping | none needed |
| Lemma 4.1: Schwinger closed form for the BPST scale-derivative cross-term | New | `lemma_ct_rust/`; `lemma_5_2_schwinger.py`; `lemma_ct_higher_order.sage` |
| Lemma 4.2: off-diagonal cross-block bounds | New; proved in Appendix A | `lemma_od_position_terms.sage` |
| Theorem 4.3: codimension-j Uhlenbeck collar essential bottom = 4j/r^2 | New; quantitative version of the qualitative product picture in Groisser-Parker and Donaldson-Kronheimer §7.3 | via Lemmas 4.1, 4.2 |
| Theorem 4.4: Mourre spectral type (no sigma_sc above threshold, locally finite sigma_pp, AC spectrum) | New | `mourre_estimate_cusp.sage`; `mourre_iterated_commutator.sage` |
| Remark 4.5: SU(N) extension | Routine; recorded as a remark | `lemma_ct_SUN.sage` |

## What is open

The global bottom lambda_0(M_k(S^4_r)) for k >= 2 is open. The paper proves the collar-localized essential-spectrum statement only, and says so (Remark on domain monotonicity; Open Question in §5). Two structural reductions of the k = 2 case appear in the companion as conditional propositions.
