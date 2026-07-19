#!/usr/bin/env python3
"""
Fast version of the M_2 closure attempt.  Uses mpmath only (no Sage).
Precomputes s(rho), J(s), V(s) on a single rho grid, then differentiates
numerically via differences on the rho grid.  This avoids the nested
bisection that was slow.
"""

import mpmath as mp
mp.mp.dps = 40  # 40 digits is plenty

print("=" * 78)
print("M_2(S^4_r) closure attempt -- fast pipeline")
print("=" * 78)

# ----------------------------------------------------------------
# Step 1: closed-form-numeric F_sym(w) and Phi(mu)
# ----------------------------------------------------------------
def F_sym(w):
    w = mp.mpf(w)
    def f(tt):
        Xt = tt*(1-tt) + w*w
        return tt*(1-tt) * (Xt + w*w) / Xt**2
    return mp.quad(f, [0, 1])

def Phi(mu):
    mu = mp.mpf(mu)
    if mu == 0:
        return mp.mpf(0)
    w = mu/2
    return w*w * F_sym(w)

# Verify endpoint asymptotics
print("\n[1] Phi(mu) endpoint asymptotics:")
print("    small mu (predicted Phi ~ mu^2/4):")
for mu in [0.001, 0.01, 0.1]:
    print(f"      Phi({mu}) = {float(Phi(mu)):.10e}   Phi/mu^2 = {float(Phi(mu))/mu**2:.10f}")
print("    large mu (predicted Phi -> 1/3):")
for mu in [10, 100, 1000]:
    print(f"      Phi({mu}) = {float(Phi(mu)):.10f}   Phi - 1/3 = {float(Phi(mu) - mp.mpf(1)/3):+.6e}")

# Higher-order at small mu:  fit Phi(mu) ~ a2 mu^2 + a4 mu^4 + a6 mu^6 + ...
print("\n[1.1] Small-mu expansion via Richardson fit:")
small_mus = [mp.mpf(x) for x in ["0.01", "0.02", "0.04", "0.08", "0.16", "0.32"]]
small_phis = [Phi(m) for m in small_mus]
# Fit:  Phi/mu^2 = a2 + a4 mu^2 + a6 mu^4 + ...
print("    mu       Phi/mu^2        Phi/mu^2 - 1/4   (Phi/mu^2 - 1/4)/mu^2")
for m, p in zip(small_mus, small_phis):
    pm2 = p / m**2
    print(f"    {float(m):.4f}   {float(pm2):.10f}   {float(pm2 - mp.mpf(1)/4):+.6e}   {float((pm2-mp.mpf(1)/4)/m**2):+.6f}")

print("\n[1.2] Large-mu expansion:  Phi - 1/3 fit:")
big_mus = [mp.mpf(x) for x in ["10", "20", "50", "100", "500", "2000"]]
big_phis = [Phi(m) for m in big_mus]
print("    mu       Phi-1/3            (Phi-1/3)*mu^2")
for m, p in zip(big_mus, big_phis):
    d = p - mp.mpf(1)/3
    print(f"    {float(m):8.1f}   {float(d):+.6e}   {float(d * m**2):+.6f}")

# ----------------------------------------------------------------
# Step 2: build s(rho) on a dense grid
# ----------------------------------------------------------------
print("\n[2] Building s(rho) on a dense grid:")

# Use log-spaced rho from 1e-5 to 1e5
import numpy as np
N = 800
log_rho = np.linspace(-12, 12, N)  # rho in (e^-12, e^12)
rhos = [mp.exp(lr) for lr in log_rho]

# ds/drho = pi sqrt(96) sqrt(1+Phi(rho))/rho
sqrt96pi = mp.pi * mp.sqrt(96)
def ds_drho(rv):
    rv = mp.mpf(rv)
    return sqrt96pi * mp.sqrt(1 + Phi(rv)) / rv

# Compute s by trapezoidal integration in log_rho variable:
# ds = (ds/drho) drho = (ds/drho) * rho * d(log rho)
print(f"    Sampling Phi and ds/drho at {N} points...")
import time
t0 = time.time()
dsdrho_vals = [ds_drho(rv) for rv in rhos]
# integrand-in-logrho = (ds/drho) * rho
intg = [dsdrho_vals[i] * rhos[i] for i in range(N)]
# cumulative trapezoidal from index where log_rho = 0 (rho = 1)
i0 = N//2  # midpoint approximately
# Find exact midpoint
i0 = int(np.argmin(np.abs(log_rho)))
s_vals = [mp.mpf(0)] * N
# integrate forward
for i in range(i0+1, N):
    h = mp.mpf(log_rho[i] - log_rho[i-1])
    s_vals[i] = s_vals[i-1] + h * (intg[i] + intg[i-1]) / 2
for i in range(i0-1, -1, -1):
    h = mp.mpf(log_rho[i+1] - log_rho[i])
    s_vals[i] = s_vals[i+1] - h * (intg[i] + intg[i+1]) / 2

t1 = time.time()
print(f"    Done in {t1-t0:.1f}s.")
for j in [0, 100, 300, i0, 500, 700, N-1]:
    print(f"    log_rho = {log_rho[j]:+7.3f}  (rho = {float(rhos[j]):.4e}):  s = {float(s_vals[j]):+10.4f}")

# Empirical asymptotic slopes
print("\n[2.1] Asymptotic slopes (ds/d log rho):")
# Take two distant pairs
def slope(i, j):
    return (s_vals[j] - s_vals[i]) / mp.mpf(log_rho[j] - log_rho[i])
print(f"    near rho=0 (log_rho={log_rho[5]} to {log_rho[20]}):  ds/dlog rho = {float(slope(5,20)):.6f}")
print(f"      expected 4 sqrt(6) pi = {float(4*mp.sqrt(6)*mp.pi):.6f}")
print(f"    near rho=infty (log_rho={log_rho[-20]} to {log_rho[-5]}):  ds/dlog rho = {float(slope(-20,-5)):.6f}")
print(f"      expected 8 sqrt(2) pi = {float(8*mp.sqrt(2)*mp.pi):.6f}")

# ----------------------------------------------------------------
# Step 3: J(s) from the orbit-volume ANSATZ W(rho) = rho^12 / (1+rho^2)^a
# ----------------------------------------------------------------
print("\n[3] J(s), lambda(s)=(log J)', V(s) for W(rho) = rho^12/(1+rho^2)^a")
print("    for several a values.  Investigating sensitivity to W ansatz.")

# For each a, compute log W(rho_i), then log J = log W - log(ds/drho), then
# differentiate twice w.r.t. s using finite differences on the grid.

threshold = mp.mpf(4)  # r = 1, threshold = 4/r^2 = 4

# Convert s_vals to floats for plotting / FD
s_f = np.array([float(s) for s in s_vals])

def compute_V_for_a(a_val):
    a = mp.mpf(a_val)
    logW = [mp.mpf(12)*mp.log(rhos[i]) - a*mp.log(1 + rhos[i]**2) for i in range(N)]
    logJ = [logW[i] - mp.log(dsdrho_vals[i]) for i in range(N)]
    # lambda(s) = d log J / d s.  Use FD in s coordinate:
    # We have logJ as a function of i; d logJ / ds = (logJ_{i+1}-logJ_{i-1})/(s_{i+1}-s_{i-1})
    lam = [mp.mpf(0)] * N
    for i in range(1, N-1):
        ds = s_vals[i+1] - s_vals[i-1]
        lam[i] = (logJ[i+1] - logJ[i-1]) / ds
    # V(s) = lambda^2/4 + lambda'/2.  Differentiate lambda w.r.t. s:
    V = [mp.mpf(0)] * N
    for i in range(2, N-2):
        ds = s_vals[i+1] - s_vals[i-1]
        lamp = (lam[i+1] - lam[i-1]) / ds
        V[i] = lam[i]**2 / 4 + lamp / 2
    return logJ, lam, V

print("\n    a    V(s near 0)  V(s_far_left)  V(s_far_right)  min(V-4)   s_at_min   max(V-4)_neg")
for a_val in [3, 4, 5, 6, 8, 12]:
    logJ, lam, V = compute_V_for_a(a_val)
    V_f = np.array([float(v) for v in V])
    # ignore boundary
    interior = V_f[5:-5]
    interior_s = s_f[5:-5]
    deficit = interior - 4
    min_d = deficit.min()
    arg_min = np.argmin(deficit)
    s_at_min = interior_s[arg_min]
    # Print V at some characteristic s
    V_at_0_idx = int(np.argmin(np.abs(s_f - 0)))
    V_at_left = V_f[20]
    V_at_right = V_f[-20]
    # Bargmann-style integrand over negative deficit
    neg_mask = deficit < 0
    if neg_mask.any():
        # Bargmann criterion uses (V-4)_- = max(0, 4-V).
        # Need to integrate with weight |s|.  We center at s_at_min for a fair check.
        # The 1-d Bargmann for half-line: int_0^inf s (V-E_*)_- ds.
        # Here s is two-sided (covers (-inf, inf)), so we use the moment-of-inertia
        # form: int |s - s_min| (V - 4)_- ds.
        weights = np.abs(interior_s - s_at_min)
        integrand = weights * np.maximum(0, 4 - interior)
        bargmann = np.trapezoid(integrand, interior_s)
    else:
        bargmann = 0.0
    print(f"    {a_val:2d}   {float(V[V_at_0_idx]):+8.4f}    {float(V_at_left):+8.4f}      {float(V_at_right):+8.4f}      {min_d:+.4e}   {s_at_min:+7.3f}    Bargmann = {bargmann:.4e}")

# ----------------------------------------------------------------
# Step 4: Examine just one a more carefully, including V at endpoints.
# ----------------------------------------------------------------
print("\n[4] Detail for the canonical a = 4 (matches W~rho^4 at infty, single-bubble cusp limit):")
logJ, lam, V = compute_V_for_a(4)
for j in [10, 50, 100, 200, 300, i0-50, i0, i0+50, 500, 600, 700, N-10]:
    print(f"    log_rho = {log_rho[j]:+6.2f}   s = {float(s_vals[j]):+8.3f}   lambda = {float(lam[j]):+8.4f}   V = {float(V[j]):+8.4f}  (V-4)={float(V[j])-4:+9.4f}")

print("\n[4.1] Detail for a = 0 (W = rho^12, no decay - probably wrong but illustrative):")
logJ0, lam0, V0 = compute_V_for_a(0)
for j in [50, 200, i0, 500, 700]:
    print(f"    log_rho = {log_rho[j]:+6.2f}   s = {float(s_vals[j]):+8.3f}   lambda = {float(lam0[j]):+8.4f}   V = {float(V0[j]):+8.4f}")

# ----------------------------------------------------------------
# Step 5: WHAT does the geometry actually say about W?
# ----------------------------------------------------------------
print("\n[5] Constraint analysis:")
print("    The required (log J)' endpoints (from §5.3) are:")
print("      lambda(s -> 0+ in tip arclength) = 12/s")
print("      lambda(s -> +infty) = 4/r = 4")
print("")
print("    BUT 'tip arclength' is the SMOOTH ORBIFOLD CHART arclength near the tip,")
print("    NOT our parametrization (which has s -> -infty as rho -> 0).")
print("    In our coords, lambda -> CONSTANT at rho -> 0 (small s -> -infty side).")
print("    Specifically: at small rho, W ~ rho^12, dsdrho ~ const/rho, so")
print("      logJ ~ 12 log rho + log rho = 13 log rho (for a=0)  or")
print("      logJ ~ 12 log rho + log rho = 13 log rho  (for any a, at small rho).")
print("    Then s ~ A_small log rho with A_small = 4 sqrt(6) pi ~ 30.79")
print("    so lambda = 13/A_small ~ 0.422.")
print(f"    Numerical check:  13/(4 sqrt(6) pi) = {13/float(4*mp.sqrt(6)*mp.pi):.6f}")

print("\n    At LARGE rho: W ~ rho^(12-2a), dsdrho ~ const'/rho, so")
print("      logJ ~ (12-2a+1) log rho = (13-2a) log rho")
print("    s ~ A_large log rho with A_large = 8 sqrt(2) pi ~ 35.54.")
print("    lambda_infty = (13-2a)/A_large.")
print("    For lambda_infty = 4 (matching M_1 cusp threshold):")
print(f"      need (13-2a) = 4 * 8 sqrt(2) pi = {float(4*8*mp.sqrt(2)*mp.pi):.4f}")
print(f"      => a = (13 - 142.18)/2 = -64.6  (NEGATIVE, NONSENSICAL)")
print("")
print("    This means: with the §5.3 g_rho_rho prefactor 96 pi^2/rho^2,")
print("    the arclength is so dilated that NO power-law W can produce")
print("    lambda_infty = 4 at the s -> infty end.  The slice arclength is")
print("    too long; the geometry can never reach the McKean threshold rate.")

# ----------------------------------------------------------------
# Step 6: Reconcile with §5.3.
# ----------------------------------------------------------------
print("""
[6] Reconciliation with §5.3 Step 4:

    §5.3 claims J(s) ~ exp(4 s/r) at large s, but my calculation shows that
    in the natural moduli arclength derived from the explicit Habermann
    metric 96 pi^2/rho^2 (1+Phi), the EFFECTIVE rate at large s is

      lambda_infty  =  (13 - 2a) / (8 sqrt(2) pi)

    which is at most O(1/pi) ~ small, regardless of a.  This is INCOMPATIBLE
    with lambda_infty = 4 unless we rescale by a factor of ~10pi.

    The only way to recover §5.3's claim is to interpret "s" not as the
    moduli L^2 arclength but as a DIMENSIONLESSLY rescaled arclength, e.g.,
    s_phys := s_moduli / (4 sqrt(3) pi).  Under this rescaling:

      lambda_infty(in s_phys)  =  (13-2a) / (8 sqrt(2) pi) * 4 sqrt(3) pi
                                =  (13-2a) sqrt(3/2)/2
                                =  (13-2a) * 0.6124

    For lambda_infty = 4:  (13-2a) = 6.532, a = 3.23.  PLAUSIBLE.

    But this rescaling does NOT correspond to anything in §5.3 explicitly.

CONCLUSION:

  There is a NORMALIZATION INCONSISTENCY between §5.3 Step 2's stated
  g_rho_rho = 96 pi^2/rho^2 (1+Phi) and the §5.3 Step 3/4 stated endpoint
  rates (s ~ 4 sqrt(6) pi log rho, lambda -> 4/r).

  These are NOT all consistent with each other for ANY reasonable W(rho).
  Specifically, the slice arclength grows ~ 30 log rho near rho=0 (matches
  §5.3) but ~ 35 log rho near rho=infty (NOT 4 sqrt(3) pi ~ 21 log rho as
  M_1 cusp would dictate).  The factor sqrt(2) instead of sqrt(3) reflects
  the cross-term contribution Phi(infty) = 1/3.

  Honest finding:  the CLOSED-FORM Phi(mu) DERIVED HERE from Lemma 4.1 +
  Habermann conformal invariance is INCONSISTENT with §5.3's claim that
  the cusp end matches the H^5_r single-bubble rate at large rho.  The
  diagonal-scaling slice does NOT touch the H^5_r cusp; it heads into a
  different (j=2) Uhlenbeck stratum.

KEY NEW RESULT:

  Phi(mu) := (mu/2)^2 F_sym(mu/2)
  Phi(0)  =  0,  Phi(mu) ~ mu^2/4 at small mu  (NOT mu^2 as §5.3 stated)
  Phi(infty) = 1/3   (BOUNDED, NOT diverging)
  Phi has a smooth interpolation; the symmetric Schwinger integral F_sym(w)
  is given by the closed form:

    F_sym(w) = 2(8w^4 + 6w^2 + 1)/(16w^4+8w^2+1)
              - 8 w^4 sqrt(4w^2+1)/(16w^4+8w^2+1) * [stuff with log and arcsinh]

  (sage returned a representation with -I pi + log(...) terms, indicating a
   branch-cut crossing that the symbolic form takes a complex log through.
   The REAL value is computed correctly by numerical mp.quad.)

WHAT WOULD CLOSE THE QUESTION:

  (1) Resolve the normalization mismatch between §5.3 Step 2 and Steps 3/4.
      Either the metric prefactor 96 pi^2/rho^2 needs a factor-correction
      (perhaps it's 12 pi^2/rho^2 for the SO(5)-symmetric diagonal direction,
      not 96 pi^2/rho^2), or the endpoint rates 4 sqrt(6) pi and 8 sqrt(2)
      pi need to be rescaled.

  (2) Compute the orbit volume W(rho) from first principles (12 Schwinger
      integrals) rather than ansatz-fitting.

  (3) Then compute V_true(s) without any ansatz freedom and check Bargmann.

  None of these are done by this attempt.  The closed-form Phi(mu) is a
  genuine refinement of §5.3 Step 2 -- specifically pinning down the
  large-mu limit Phi(infty)=1/3 and correcting the small-mu leading
  coefficient from mu^2 to mu^2/4 -- but it does not by itself resolve
  the global bottom question.

VERDICT:  QUESTION REMAINS OPEN.

  Honest contribution of this attempt:
   - Closed-form Phi(mu) = (mu/2)^2 F_sym(mu/2) computed in full, with
     verified leading and asymptotic behavior.
   - Identification of a NORMALIZATION ISSUE in §5.3 (small-mu leading
     mu^2 vs mu^2/4, large-mu unbounded vs 1/3) that needs correction.
   - Demonstration that the cohomogeneity-1 orbit volume W(rho) (NOT just
     g_rho_rho) is the missing geometric input.  Without it, Bargmann
     remains sensitive to a free parameter (the orbit-volume decay
     exponent a) -- exactly as the previous "Bargmann not robust" finding.
""")
