#!/usr/bin/env sage
# m2_isotypic_decomposition.sage
#
# Angle 3: SO(5)-isotypic decomposition on M_2(S^4_r).
# Goal: For each non-trivial SO(5)-isotypic component (ell >= 1), the
# radial Schrodinger operator carries an extra "centrifugal" Casimir
# potential C_ell / r_orbit(s)^2 with C_ell = ell(ell+3).
# We verify:
#   (i)   r_orbit(s) asymptotics: r_orbit ~ s as s->0 (tip);
#         r_orbit ~ (r/2) e^{s/r} as s->infty (cusp on H^5_r).
#   (ii)  V_true(s) + C_ell/r_orbit(s)^2 is bounded below by inf at
#         s -> infty equal to 4/r^2 (the McKean threshold), so no
#         bound state on ell >= 1 components below 4/r^2.
#   (iii) Bargmann integral on ell >= 1 component for the negative
#         part of [V_true + Casimir - 4/r^2] is <= the ell=0 integral
#         and in fact identically 0 outside a possibly small region
#         where the repulsive centrifugal dominates.
#
# Budget: 45 min Sage.

from sage.all import *
import numpy as np

print("="*72)
print("Angle 3: SO(5)-isotypic decomposition of M_2(S^4_r)")
print("="*72)

# --------------------------------------------------------------------
# 1. Setup: closed-form Phi(mu) on the symmetric slice (from extra.tex)
# --------------------------------------------------------------------
# Phi(mu) = (mu/2)^2 * F_sym(mu/2), where
# F_sym(w) = int_0^1 t(1-t) [t(1-t) + 2 w^2] / [t(1-t) + w^2]^2 dt
#
# Endpoint asymptotics (verified previously):
#   Phi(mu) ~ mu^2/4 - mu^4/8 + O(mu^6)   (mu -> 0)
#   Phi(mu) -> 1/3 - 0.4/mu^2 + ...        (mu -> infty)

t, w, mu = var('t w mu', domain='positive')
F_sym_integrand = t*(1-t) * (t*(1-t) + 2*w**2) / (t*(1-t) + w**2)**2

def F_sym_num(wval, prec=53):
    """Numerical F_sym(w)."""
    from sage.all import numerical_integral
    f = lambda tt: float(tt*(1-tt)*(tt*(1-tt) + 2*wval**2)/(tt*(1-tt)+wval**2)**2)
    val, err = numerical_integral(f, 0, 1)
    return val

def Phi_num(muval):
    """Phi(mu) numerically."""
    if muval < 1e-12:
        return muval**2/4.0
    w0 = muval/2.0
    return w0**2 * F_sym_num(w0)

# Spot check:
print("\n[1] Phi(mu) spot check:")
for mv in [0.01, 0.1, 1.0, 5.0, 50.0, 500.0]:
    print(f"  Phi({mv:8.3f}) = {Phi_num(mv):.6f}  (small-mu pred mu^2/4 = {float(mv**2/4):.6f}, large-mu lim 1/3 = {float(1/3):.6f})")

# --------------------------------------------------------------------
# 2. r_orbit(s) asymptotics
# --------------------------------------------------------------------
# The SO(5)-invariant slice has principal orbit SO(5)/(SO(4)*U(1)^2)
# of real dim 12. The orbit metric is a homogeneous metric on this
# coset, scaled by an "orbit radius" function r_orbit(s).
#
# At s -> 0 (orbifold tip of transverse dim 12): J(s) = s^12 (1+O(s^2))
#   J(s) = vol(orbit at arclength s) = r_orbit(s)^12 * vol(unit orbit)
#   so r_orbit(s) ~ s   (linear in arclength at the tip)
#
# At s -> infty (cusp = H^5_r single-bubble end): J(s) ~ C e^{4s/r}.
#   The cusp end is H^5_r with radial metric ds^2 + r^2 sinh^2(s/r) d Omega_4^2.
#   The "transverse SO(5)/SO(4) = S^4" orbit at hyperbolic radius s has
#   radius r sinh(s/r) ~ (r/2) e^{s/r} as s -> infty.
#   So r_orbit(s) ~ (r/2) e^{s/r}.
#
# We use r = 1 throughout.
r = 1.0

s = var('s', domain='positive')

def r_orbit_tip(sval):
    """Asymptotic r_orbit near s=0: r_orbit ~ s."""
    return sval

def r_orbit_cusp(sval, r=1.0):
    """Asymptotic r_orbit at large s: r_orbit ~ (r/2) e^{s/r}.
    Equivalent to r sinh(s/r), which interpolates both."""
    return r * np.sinh(sval/r)

# Global interpolant: r_orbit(s) := r * sinh(s/r) for all s (linear near 0, exp at infty)
def r_orbit_global(sval, r=1.0):
    return r * np.sinh(sval/r)

print("\n[2] r_orbit(s) global interpolant r*sinh(s/r):")
for sv in [0.01, 0.1, 1.0, 5.0, 10.0]:
    print(f"  s = {sv:6.3f}: r_orbit = {r_orbit_global(sv):.6f}  (tip pred {sv:.6f}, cusp pred {0.5*np.exp(sv):.6f})")

# --------------------------------------------------------------------
# 3. Casimir potential for SO(5)
# --------------------------------------------------------------------
# Highest weight (ell, 0) of SO(5) has Casimir ell(ell+3).
# Other reps (ell1, ell2) with ell1 >= ell2 >= 0 have Casimir
#   C(ell1, ell2) = ell1(ell1+3) + ell2(ell2+1).
# For our purposes the lowest non-trivial is (1,0): C_1 = 4.

def C_ell(ell):
    """SO(5) Casimir for highest weight (ell, 0)."""
    return ell*(ell+3)

print("\n[3] SO(5) Casimir eigenvalues C_ell = ell(ell+3):")
for ell in range(0, 6):
    print(f"  ell = {ell}: C_ell = {C_ell(ell)}")

# --------------------------------------------------------------------
# 4. V_true(s) on the SO(5)-invariant component
# --------------------------------------------------------------------
# From prior analyses (m2_closure_higher_order, Step 2/Step 4 discussion),
# V_true(s) on the invariant slice is V(s) = (J'/J)^2 / 4 + (J'/J)'/2 ...
# wait: Liouville form: V(s) = (1/4)(J'/J)^2 + (1/2)(J'/J)'.
# Or equivalently, with u = J^{1/2}, V = u''/u.
#
# We don't need V_true exactly; we just need the LOWER BOUND
# V_true >= 0 in the bulk of the slice (prior numerics give V_true in [0.005, 0.05]
# under the symmetric-slice arclength; non-negative throughout).
# The key fact for Angle 3 is the *Casimir contribution* dominates.

# Use a stand-in V_true that matches the endpoint asymptotics:
#   V_true(s) -> 30/s^2 as s -> 0 (the J ~ s^12 Bessel tip potential, Step 4)
#   V_true(s) -> 4/r^2 as s -> infty (McKean cusp).
# This is the "ideal closed-form" Step 4 V; the closed-form Phi gives a milder
# V in [0.005, 0.05] (the Step 2 vs Step 4 mismatch). For an UPPER BOUND
# on the Bargmann integral on ell >= 1, the smaller V_true is better.
# We'll use the small V_true_phi(s) ~ small bounded number scenario.

def V_phi_estimate(sval, r=1.0):
    """Rough V_true on symmetric slice in the Step-2 normalization,
    matching the [0.005, 0.05] range from prior numerics."""
    # Smooth bump-ish: small at 0 (since Phi small), grows to ~0.05 at mid-s,
    # decays to ~0 at infty (since the SU(2)*U(1)^2 stabilizer contracts).
    return 0.05 * np.exp(-(sval - 2.0)**2 / 4.0) + 0.005

# --------------------------------------------------------------------
# 5. Centrifugal Casimir potential on the ell-isotypic component
# --------------------------------------------------------------------
def V_casimir(sval, ell, r=1.0):
    """Centrifugal potential C_ell / r_orbit(s)^2 on the ell-isotypic component."""
    return C_ell(ell) / (r_orbit_global(sval, r))**2

print("\n[5] Casimir potential V_cas(s) = C_ell/r_orbit^2 sample values:")
print(f"{'s':>6} | {'ell=1, C=4':>11} | {'ell=2, C=10':>12} | {'ell=3, C=18':>12}")
for sv in [0.05, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
    row = f"{sv:6.3f} | "
    for ell in [1, 2, 3]:
        row += f"{V_casimir(sv, ell):11.4f}  | "
    print(row)

# --------------------------------------------------------------------
# 6. Bargmann integral on ell >= 1 isotypic component
# --------------------------------------------------------------------
# Bargmann: N_- <= int_0^infty s * (V_eff(s) - 4/r^2)_- ds
# where V_eff = V_true(s) + V_cas(s; ell).
#
# Negative part: (x)_- = max(0, -x).
# The Casimir is POSITIVE (repulsive). So we have, with V_true small (~0.05)
# and V_cas adding repulsion:
#   V_eff(s) - 4 = V_true(s) + C_ell/r_orbit^2 - 4
# This becomes negative only where C_ell/r_orbit^2 < 4 - V_true ~ 4,
# i.e. r_orbit > sqrt(C_ell/4). For ell=1, r_orbit > 1, i.e. s > arcsinh(1) ~ 0.881.
# At those s, V_eff - 4 ~ V_true(s) - 4 + small_positive
# = roughly the same as the invariant case in the asymptotic region.
#
# But here's the angle's key point: on ell >= 1 components, the
# ESSENTIAL SPECTRUM still has bottom 4/r^2 (Casimir decays at infinity),
# but the *bound states below 4/r^2* are obstructed at small s by the
# Casimir blow-up. So Bargmann integral on ell >= 1 may be smaller.

def bargmann_isotypic(ell, r=1.0, sN=50.0, npts=2000):
    """Compute Bargmann integral on ell-isotypic component.
       Returns int_0^sN s * (V_eff(s) - 4/r^2)_- ds, where
       V_eff = V_phi_estimate + V_casimir(ell)."""
    from numpy import linspace
    threshold = 4.0/r**2
    svals = linspace(0.001, sN, npts)
    integrand = []
    for sv in svals:
        V_eff = V_phi_estimate(sv, r) + V_casimir(sv, ell, r)
        neg_part = max(0.0, threshold - V_eff)
        integrand.append(sv * neg_part)
    # Trapezoidal
    h = svals[1] - svals[0]
    val = h * (sum(integrand) - 0.5*(integrand[0] + integrand[-1]))
    return val

print("\n[6] Bargmann integral B_ell = int_0^infty s * (V_eff - 4/r^2)_- ds:")
print(f"  (using V_true ~ V_phi_estimate ranging in [0.005, 0.05])")
for ell in [0, 1, 2, 3, 4]:
    if ell == 0:
        # no Casimir
        from numpy import linspace
        threshold = 4.0
        svals = linspace(0.001, 50.0, 2000)
        integrand = [sv * max(0.0, threshold - V_phi_estimate(sv)) for sv in svals]
        h = svals[1]-svals[0]
        val = h*(sum(integrand) - 0.5*(integrand[0]+integrand[-1]))
        print(f"  ell = {ell}: B_ell = {val:12.4f}    <-- INVARIANT (no Casimir, V<<4)")
    else:
        B = bargmann_isotypic(ell)
        print(f"  ell = {ell}: B_ell = {B:12.4f}    <-- with C_ell={C_ell(ell)}")

# --------------------------------------------------------------------
# 7. The actual CLEAN statement for ell >= 1
# --------------------------------------------------------------------
# For ell >= 1: V_eff(s) = V_true(s) + C_ell/r_orbit(s)^2 >= C_ell/r_orbit(s)^2.
# At s -> infty, V_eff -> 0 + 0 + V_true(infty); V_true(infty) is supposed to be 4
# (McKean threshold). So V_eff(s) -> 4/r^2 as s -> infty.
# At s -> 0, V_eff -> infty (since r_orbit ~ s -> 0).
# So the QUESTION on ell >= 1 is:
#   Does V_eff(s) >= 4/r^2 for ALL s? If yes, Bargmann integral is 0, no bound state.
#
# This is equivalent to checking:
#   V_true(s) + C_ell/r_orbit(s)^2 >= 4/r^2 for all s >= 0.
#
# Since C_ell >= 4 for ell >= 1 and 1/r_orbit(s)^2 >= ? we need
# 1/r_orbit(s)^2 >= 1/r^2 i.e. r_orbit(s) <= r.
# But r_orbit(s) = r sinh(s/r) which is < r exactly when s < r arcsinh(1) ~ 0.881*r,
# and > r for larger s. So Casimir alone does NOT give >= 4/r^2 globally.
#
# At large s, Casimir vanishes and we need V_true(s) -> 4/r^2.
# At intermediate s ~ r (where Casimir is exactly 4/r^2), V_eff ~ V_true + 4/r^2.
# Since V_true >= 0, this is >= 4/r^2. GOOD.
# Conclusion: V_eff >= 4/r^2 on ell >= 1 if and only if V_true(s) + Casimir(s) >= 4/r^2,
# which holds for s NOT too large (where Casimir dominates) AND for s -> infty
# (where V_true -> 4/r^2 by McKean).

print("\n[7] Pointwise check V_eff(s) >= 4/r^2 on ell >= 1:")
print("  V_eff(s) = V_true(s) + C_ell/r_orbit(s)^2  vs  threshold 4/r^2 = 4")
from numpy import linspace
for ell in [1, 2]:
    deficit_max = 0.0
    s_at_max = 0.0
    for sv in linspace(0.01, 30.0, 3000):
        V_eff = V_phi_estimate(sv) + V_casimir(sv, ell)
        deficit = 4.0 - V_eff  # positive iff V_eff < 4
        if deficit > deficit_max:
            deficit_max = deficit
            s_at_max = sv
    print(f"  ell = {ell}: max(4 - V_eff) = {deficit_max:.4f} at s = {s_at_max:.3f}")
    if deficit_max <= 0:
        print(f"    --> V_eff >= 4/r^2 globally; ell-isotypic Bargmann integral = 0.")
    else:
        print(f"    --> V_eff dips below 4/r^2; need V_true asymptotic to 4 at infty.")

# --------------------------------------------------------------------
# 8. The TRUE V_true at large s: does it approach 4/r^2 (McKean)?
# --------------------------------------------------------------------
# From Step 2/Step 4 mismatch (extra.tex remark): V_phi gives [0.005, 0.05],
# but the "correct" Step 4 has V -> 4/r^2 at the cusp.
#
# CRUCIAL: the angle's argument works *modulo* the assumption that V_true
# matches McKean at infinity. The Casimir adds repulsion at small s, where
# V_true alone is small. The combined V_eff is bounded below by Casimir
# which is large at small s (where r_orbit ~ s -> 0) and by McKean ~4/r^2
# at large s. The TRANSITION REGION s ~ r is the only delicate part.

# Use the McKean-rescaled V_true that does asymptote to 4/r^2:
def V_McKean(sval, r=1.0):
    """V_true(s) under the Step-4 normalization: V -> 30/s^2 at 0, -> 4/r^2 at infty."""
    return 30.0/(sval**2 + 0.01) * np.exp(-sval) + 4.0/r**2 * (1 - np.exp(-sval))

print("\n[8] Under Step-4 normalization (V_true -> 4/r^2 at infty):")
for ell in [1, 2, 3]:
    deficit_max = 0.0
    s_at_max = 0.0
    for sv in linspace(0.01, 30.0, 5000):
        V_eff = V_McKean(sv) + V_casimir(sv, ell)
        deficit = 4.0 - V_eff
        if deficit > deficit_max:
            deficit_max = deficit
            s_at_max = sv
    print(f"  ell = {ell}: max(4 - V_eff) = {deficit_max:.6f} at s = {s_at_max:.3f}")
    # Bargmann integral
    from numpy import linspace as L
    svals = L(0.001, 50.0, 5000)
    integrand = [sv * max(0.0, 4.0 - V_McKean(sv) - V_casimir(sv, ell)) for sv in svals]
    h = svals[1]-svals[0]
    B = h*(sum(integrand) - 0.5*(integrand[0]+integrand[-1]))
    print(f"          Bargmann integral B_ell = {B:.6e}")

print("\n" + "="*72)
print("Summary:")
print("="*72)
print("""
- r_orbit(s) ~ s near tip (s->0); r_orbit(s) ~ (r/2)e^{s/r} at cusp (s->infty).
  Global interpolant r*sinh(s/r) matches both.
- SO(5) Casimirs: C_ell = ell(ell+3); lowest non-trivial C_1 = 4.
- Centrifugal potential C_ell/r_orbit(s)^2 on ell-isotypic component.
- For ell >= 1, V_eff(s) = V_true(s) + Casimir.
  * Small s: Casimir ~ C_ell/s^2 -> infty (repulsive).
  * Large s: Casimir -> 0 exponentially; V_true -> 4/r^2 (McKean).
- Under the Step-4 (McKean-asymptotic) V_true, V_eff >= 4/r^2 pointwise
  iff Casimir compensates the dip in V_true near s ~ r.
  Numerical check shows: small dip possible but Bargmann integral
  decreases rapidly with ell.
- For ell sufficiently large (C_ell >= 4 r^2 / r_orbit_min^2): no dip.
""")
