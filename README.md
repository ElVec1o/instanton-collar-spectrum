# instanton-collar-spectrum

[![DOI](https://zenodo.org/badge/1248611476.svg)](https://doi.org/10.5281/zenodo.20370060)

LaTeX source, PDF, and verification code for

> V. Bonfioli, *The L^2 spectral geometry of the charge-one SU(2) instanton moduli space on S^4*, 2026. 5 pages. [`paper/CORRECTED/corrected_core.pdf`](paper/CORRECTED/corrected_core.pdf)

## Results

M_1(S^4_r) is the charge-one SU(2) instanton moduli space over the round 4-sphere of radius r, carrying the L^2 (Weil-Petersson) metric.

By Groisser-Parker and Doi-Matsumoto-Matumoto, (M_1(S^4_r), g_L2) is the interior of a compact smooth Riemannian 5-manifold with boundary: an SO(5)-invariant cohomogeneity-one metric on the closed 5-ball, of finite volume, strictly positive non-constant sectional curvature, whose boundary is a totally geodesic round S^4 (the Uhlenbeck bubbling locus) at finite distance. The hyperbolic metric on the same 5-ball is the information (Fisher-Rao) metric of Groisser-Murray, which is a different metric.

- **Theorem 2.1.** The L^2-Laplacian on the open moduli space is symmetric and bounded below but not essentially self-adjoint. Its self-adjoint realizations are the elliptic boundary conditions on the bounding S^4, and every one has compact resolvent, hence purely discrete spectrum, with leading Weyl asymptotics and explicit constant. No continuous spectrum. (Contrast: on the Riemann moduli space the Weil-Petersson Laplacian *is* essentially self-adjoint, Ji-Mazzeo-Muller-Vasy.)
- **Theorem 3.1.** Separation by SO(5) angular momentum reduces the eigenvalue problem to an explicit family of regular-singular Sturm-Liouville problems, one per harmonic degree k, with multiplicities 1, 5, 14, 30, 55, ... The Friedrichs (Dirichlet) ground state is simple, SO(5)-invariant, and lies in [0.8203116, 0.8203177] r^-2. The Neumann spectral gap lies in [0.9988337, 0.9988386] r^-2 and is the five-fold degenerate k = 1 mode, in the vector representation of SO(5). That this multiplet rather than any radial mode is the gap is proved, not merely computed.
- **Exact geometry.** The boundary 4-volume is exactly (2/3) pi^6 r^4.
- **Lemma 3.5.** The BPST scale-derivative L^2 cross-term in elementary closed form, with leading asymptotic 48 pi^2 s_1 s_2 / R^2.

Eigenvalues are reported as certified enclosures. Upper bounds are Rayleigh-Ritz, from a conforming consistent-mass discretization. Lower bounds come from the Liouville normal form: since p w = g^4, the substitution s = arc length, v = g u sends the separated problem to -v'' + Q_k v = lam v with Q_k = k(k+3)/g + g_ss/g, which is increasing in s and in k, so per-cell left-endpoint values are an exact minorant and min-max gives rigorous lower bounds. Monotonicity in k plus the certified separations from the k = 0 overtone (>= 2.0351863) and the k = 2 ground state (>= 2.3973341) prove the gap is the k = 1 multiplet. No claim is made about closed forms: at the precision available, inverse-symbolic search returns spurious matches and is evidence in neither direction.

## Verification

| Artifact | Checks | Run |
|---|---|---|
| `verification/M1_L2_spectrum_allk.py` | Theorem 3.1: SO(5)-separated spectrum, consistent-mass FEM | `python3 M1_L2_spectrum_allk.py` |
| `verification/M1_L2_enclosures.py` | Theorem 3.1: two-sided FEM brackets | `python3 M1_L2_enclosures.py` |
| `verification/M1_L2_certified_bounds.py` | Theorem 3.1: Liouville form, monotonicity of Q_k, certified lower bounds | `python3 M1_L2_certified_bounds.py` |
| `verification/M1_L2_certified_tight.py` | Theorem 3.1: refined certified enclosures (2e5 cells) | `python3 M1_L2_certified_tight.py 200000` |
| `verification/M1_L2_volume_enclosure.py` | Certified volume enclosure via monotone Riemann sums | `python3 M1_L2_volume_enclosure.py` |
| `verification/M1_L2_geometry_weyl.py` | Volume, boundary 4-volume, Weyl constants | `python3 M1_L2_geometry_weyl.py` |
| `verification/lemma_ct_rust/` | Lemma 3.5: 4-D integral vs closed form, 106-point grid | `cargo run --release` |
| `verification/lemma_ct_elementary.py` | Lemma 3.5: elementary form, symbolic + 50-digit quadrature | `python3 lemma_ct_elementary.py` |
| `verification/lemma_ct_montecarlo.py` | Lemma 3.5: Monte Carlo | `python3 lemma_ct_montecarlo.py` |
| `verification/lemma_ct_higher_order.sage` | Lemma 3.5: next-order coefficient | `sage lemma_ct_higher_order.sage` |
| `verification/lemma_ct_SUN.sage` | SU(N) variant | `sage lemma_ct_SUN.sage` |

## Building

```
cd paper/CORRECTED && pdflatex corrected_core.tex && pdflatex corrected_core.tex
```

Expected: 6 pages, no undefined references.

## Superseded material

`paper/old/` and `verification/old/` hold versions through v1.3. Their spectral results identified the L^2 metric on M_1(S^4_r) with hyperbolic 5-space. That identification is incorrect: the hyperbolic metric on the 5-ball is the information metric, not the L^2 metric, and the L^2 metric is of finite volume and positive curvature. The results resting on it do not hold and are retained only for the record.

## Provenance

The author is an independent researcher. The manuscript was produced with extensive LLM assistance (drafting, internal adversarial review, verification tooling); the author takes responsibility for the content. It has not been reviewed by a human expert. The numerical verification is machine-checked and independent of the derivations. Corrections and counterexamples: vicobonfioli@gmail.com.

## License and citation

Paper: CC-BY 4.0 (`LICENSE-paper`). Code: MIT (`LICENSE-code`).

```
Bonfioli, Vico. The L^2 spectral geometry of the charge-one SU(2) instanton
moduli space on S^4. 2026. DOI: 10.5281/zenodo.20370060
https://github.com/ElVec1o/instanton-collar-spectrum
```
