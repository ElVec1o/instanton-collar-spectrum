#!/usr/bin/env python3
"""Monte Carlo cross-check of the Lemma 4.1 Schwinger closed form.

The cross-term integral, after the eta-tensor contraction
sum_{a,mu} eta^a_{mu nu} eta^a_{mu rho} = 3 delta_{nu rho}, reads

    I(R, s_1, s_2) = 48 s_1 s_2 * int_{R^4}  x . (x - R e_1)
                     / [(|x|^2 + s_1^2)^2 * ((x - R e_1)^2 + s_2^2)^2]  d^4 x.

Schwinger parametrization 1/A^2 = int_0^infty alpha exp(-alpha A) d alpha
(Gaussian integration over R^4, completing the square, repartitioning to
(sigma, t)) gives the 1-D form of Lemma 4.1:

    I = 48 pi^2 s_1 s_2 * int_0^1 t(1-t) [2/X(t) - t(1-t) R^2 / X(t)^2] dt,
    X(t) = t(1-t) R^2 + (1-t) s_1^2 + t s_2^2,

strictly positive on [0,1] for s_1, s_2 > 0. Leading asymptotic for
s_1, s_2 << R: I = 48 pi^2 s_1 s_2 / R^2 + O((s/R)^4 |log(s/R)|).

Sampling note: the integrand's mass is dominated by the broad region
BETWEEN the two instanton centers, not by the two basins near them, so a
correct Monte Carlo estimator must cover that middle region. Uniform-box
sampling over [-L, L]^4 with L a few multiples of max(R, s_1, s_2) does;
importance sampling concentrated at the two basins does not, and
underestimates |I| by orders of magnitude. The elementary evaluation of
the 1-D integral is in lemma_ct_elementary.py; the machine-precision
comparison against the 4-D integral is in lemma_ct_rust/.

This script checks, independently of both:
  [1] the 1-D quadrature against the leading asymptotic 48 pi^2 s1 s2/R^2;
  [2] the 1-D quadrature against uniform-box Monte Carlo of the 4-D integral;
  [3] linear vanishing in s_1 at fixed (R, s_2)  (the decoupling rate used
      by Theorem 4.3);
  [4] the R^{-2} large-separation scaling.
"""

import numpy as np
from scipy import integrate

PI = np.pi

def I_analytic(R, s1, s2):
    """The exact Schwinger-derived integral, computed as a 1-D quadrature."""
    def F(t):
        X = t * (1-t) * R**2 + (1-t) * s1**2 + t * s2**2
        return t * (1-t) * (2.0 / X - t * (1-t) * R**2 / X**2)
    val, _ = integrate.quad(F, 0, 1)
    return 48 * PI**2 * s1 * s2 * val

def I_leading(R, s1, s2):
    """Leading large-R asymptotic: 48 pi^2 s_1 s_2 / R^2."""
    return 48 * PI**2 * s1 * s2 / R**2

def MC_uniform_box(R, s1, s2, box_extent=4.0, n_samples=10_000_000, seed=0):
    """Uniform Monte Carlo over [-L, L]^4, L = box_extent * max(R, s1, s2).
    Covers the broad middle region between the two centers that carries
    the integral's mass."""
    rng = np.random.default_rng(seed)
    L = box_extent * max(R, s1, s2)
    pts = rng.uniform(-L, L, size=(n_samples, 4))
    vol = (2 * L)**4

    e1 = np.array([1.0, 0, 0, 0])
    x_dot = np.einsum('Ni,Ni->N', pts, pts)
    xR = pts - R * e1[None, :]
    xR_dot = np.einsum('Ni,Ni->N', xR, xR)
    num = np.einsum('Ni,Ni->N', pts, xR)
    denom1 = (x_dot + s1*s1)**2
    denom2 = (xR_dot + s2*s2)**2
    integrand = 48.0 * s1 * s2 * num / (denom1 * denom2)
    return vol * integrand.mean(), vol * integrand.std(ddof=1) / np.sqrt(n_samples)

if __name__ == "__main__":
    print("=" * 72)
    print("Lemma 4.1 closed form via Schwinger parametrization: MC cross-check")
    print("=" * 72)

    print("\n[1] 1-D Schwinger integral vs leading 48 pi^2 s1 s2 / R^2:")
    print(f"\n    {'R':>5} {'s1':>5} {'s2':>5} {'I_analytic':>14} {'I_leading':>14} {'ratio':>8}")
    for R, s1, s2 in [(5, 1, 1), (10, 1, 1), (20, 1, 1), (50, 1, 1), (100, 1, 1),
                       (10, 0.5, 1), (10, 0.1, 1), (10, 0.01, 1),
                       (10, 1, 0.1)]:
        Ia = I_analytic(R, s1, s2)
        Il = I_leading(R, s1, s2)
        print(f"    {R:5.1f} {s1:5.2f} {s2:5.2f} {Ia:+.5e} {Il:+.5e} {Ia/Il:8.5f}")

    print("\n[2] 1-D quadrature vs uniform-box MC of the 4-D integral, R=10, s1=s2=1:")
    print("    analytic prediction 48 pi^2 / 100 (1 + O(1/R^2)) = {:.4f}".format(I_analytic(10, 1, 1)))
    for box in [2, 3, 4, 6]:
        v, e = MC_uniform_box(10, 1, 1, box_extent=box, n_samples=5_000_000, seed=42)
        print(f"    box=[-{box*10},+{box*10}]^4:  MC = {v:.4f} +/- {e:.4f}")

    print("\n[3] s_1 -> 0 scaling at R=10, s_2=1 (linear vanishing):")
    print(f"    {'s_1':>8} {'I_analytic':>14} {'I/s_1':>14}")
    for s1 in [1.0, 0.5, 0.2, 0.1, 0.05, 0.01]:
        Ia = I_analytic(10, s1, 1)
        print(f"    {s1:8.3f} {Ia:+.5e} {Ia/s1:+.5e}")

    print("\n[4] R -> infinity scaling at s_1=s_2=1:")
    print(f"    {'R':>5} {'I_analytic':>14} {'I * R^2':>14}")
    for R in [5, 10, 20, 50, 100, 200]:
        Ia = I_analytic(R, 1, 1)
        print(f"    {R:5.1f} {Ia:+.5e} {Ia*R**2:+.5e}")

    print("\n[5] Summary:")
    print("    I(R, s_1, s_2) = 48 pi^2 s_1 s_2 / R^2 + O((s/R)^4 |log(s/R)|),")
    print("    vanishing linearly in s_1 at fixed (R, s_2), as used by Theorem 4.3.")
