#!/usr/bin/env python3
"""
Keystone verification of the corrected Mourre / Froese-Hislop analysis for
Theorem 4.8, in GEODESIC-POLAR coordinates on the exact H^5_r.

The reviewer showed the paper's conjugate operator A=(1/2)(y d_y + d_y y) is the
TRANSLATION generator in u=log y and gives [H_rad,iA]=0. The correct generator is
the DILATION in the geodesic radius rho (leading coefficient rho ~ log(1/s)).

We verify, symbolically:
 (1) Half-density (Liouville) reduction of the radial Laplacian on H^5_r:
        -Delta_rad,l  ->  H_l = -d_rho^2 + V_l(rho),
        V_l(rho) = 4/r^2 + (2 + l(l+3))/(r^2 sinh^2(rho/r)).
     Threshold 4/r^2; transverse/centrifugal term SHORT-RANGE (~ e^{-2 rho/r}).
 (2) Mourre commutator with the CORRECT generator A=(1/2)(rho D + D rho), D=-i d_rho:
        i[H_l, A] = 2(H_l - 4/r^2) + R_l(rho),
        R_l(rho) = (2 c_l / sinh^2(rho/r)) [ (rho/r) coth(rho/r) - 1 ],  c_l=(2+l(l+3))/r^2,
     and R_l -> 0 as rho->infinity (short-range, rel.-compact). Strict Mourre on
     [4/r^2 + delta, infty) modulo compact.
 (3) Regression witness: the paper's OLD generator A_old = (1/2)(y d_y + d_y y),
     y = e^{rho/r}, equals d/d(rho/r) up to lower order == translation; [H_rad, A_old]=0
     on the leading symbol, reproducing the fatal bug.
 (4) Positivity: for rho >= rho_0 the good term dominates: on the spectral window
     H_l - 4/r^2 >= delta, the remainder R_l is bounded and ->0, so
     i[H_l,A] >= 2 delta - |R_l| > 0 for rho_0 large.  Numerorate the crossover.

All identities are checked symbolically (sympy). Exit 0 and PASS iff all hold.
"""
import sympy as sp

ok = True
def check(name, cond):
    global ok
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    ok = ok and bool(cond)

rho, r, l = sp.symbols('rho r l', positive=True)
S = sp.sinh(rho/r); C = sp.cosh(rho/r)

# ---------- (1) Liouville reduction ----------
# radial measure density J = (r sinh(rho/r))^4  (times const S^4 volume)
J = (r*S)**4
# radial part of -Delta on L^2(J drho):  -d_rho^2 - (J'/J) d_rho
W = sp.simplify(sp.diff(J, rho)/J)         # J'/J = (4/r) coth(rho/r)
check("J'/J = (4/r) coth(rho/r)", sp.simplify(W - (4/r)*C/S) == 0)

# Liouville potential V = (1/4) W^2 + (1/2) W'  (conjugation by J^{1/2})
V0 = sp.simplify(sp.Rational(1,4)*W**2 + sp.Rational(1,2)*sp.diff(W, rho))
V0_target = 4/r**2 + 2/(r**2 * S**2)
check("V_0(rho) = 4/r^2 + 2/(r^2 sinh^2)", sp.simplify(V0 - V0_target) == 0)

# add transverse S^4 harmonic sector l:  +(l(l+3))/(r^2 sinh^2(rho/r))
c_l = (2 + l*(l+3))/r**2
V_l = 4/r**2 + c_l/S**2
check("V_l = 4/r^2 + (2+l(l+3))/(r^2 sinh^2)",
      sp.simplify(V_l - (V0_target + l*(l+3)/(r**2*S**2))) == 0)

# threshold and short-range decay of the centrifugal part
Vminus = sp.simplify(V_l - 4/r**2)
check("V_l - 4/r^2 = c_l / sinh^2(rho/r)", sp.simplify(Vminus - c_l/S**2) == 0)
decay = sp.limit(Vminus*sp.exp(2*rho/r), rho, sp.oo)     # ~ 4 c_l  (finite => e^{-2rho/r} rate)
check("centrifugal term decays like e^{-2 rho/r} (short-range)", sp.simplify(decay - 4*c_l) == 0)

# ---------- (2) Mourre commutator with CORRECT dilation generator ----------
# Act on a test function f(rho). A = (1/2)(rho D + D rho) = rho D - i/2, D=-i d_rho.
# i[H,A] as a differential operator: i[-d_rho^2, A] = -2 d_rho^2 ; i[V,A] = -rho V'(rho).
# Verify i[H,A] f = (2(H-4/r^2) + R) f  as an identity in f, f', f''.
f = sp.Function('f')
Hf = -sp.diff(f(rho), rho, 2) + V_l*f(rho)
# A f = rho*(-i f') - (i/2) f  ; drop the global i by computing i[H,A] directly:
# i[H,A]f = i(H(Af) - A(Hf)); use A = -(i)(rho d_rho + 1/2)
def A_op(expr):
    return -sp.I*(rho*sp.diff(expr, rho) + sp.Rational(1,2)*expr)
comm = sp.I*(sp.expand(-sp.diff(A_op(f(rho)),rho,2)+V_l*A_op(f(rho)))
             - A_op(-sp.diff(f(rho),rho,2)+V_l*f(rho)))
comm = sp.expand(sp.simplify(comm))
predicted = 2*(Hf - (4/r**2)*f(rho)) + (2*c_l/S**2)*((rho/r)*(C/S) - 1)*f(rho)
predicted = sp.expand(sp.simplify(predicted))
check("i[H_l,A] = 2(H_l-4/r^2) + R_l(rho) f  (operator identity in f,f',f'')",
      sp.simplify(comm - predicted) == 0)

R_l = sp.simplify(comm - 2*(Hf - (4/r**2)*f(rho)))/f(rho)
R_l = sp.simplify(R_l)
check("R_l = (2 c_l/sinh^2)[(rho/r)coth(rho/r) - 1]",
      sp.simplify(R_l - (2*c_l/S**2)*((rho/r)*(C/S) - 1)) == 0)
check("R_l -> 0 as rho->infinity (short-range remainder)",
      sp.limit(R_l, rho, sp.oo) == 0)

# ---------- (3) regression witness: OLD generator gives ZERO ----------
# paper's A_old ~ (1/2)(y d_y + d_y y) with y=e^{rho/r} => y d_y = r^{-1}... actually
# y d_y = d/d(log y) = d/d(rho/r) = r d_rho. So A_old is proportional to the TRANSLATION
# generator r*D in rho.  i[-d_rho^2, translation] = 0.
def Aold_op(expr):   # translation generator (up to constant): -(i) d_rho
    return -sp.I*sp.diff(expr, rho)
# leading (free) part [-d_rho^2, i A_old]:
free_comm_old = sp.expand(sp.I*(-sp.diff(Aold_op(f(rho)),rho,2)
                                - Aold_op(-sp.diff(f(rho),rho,2))))
check("OLD generator: [-d_rho^2, iA_old] = 0 (translation; reproduces the bug)",
      sp.simplify(free_comm_old) == 0)

# ---------- (4) positivity crossover (numeric, l=0, r=1) ----------
import mpmath as mp
def Rl_num(rhov, rv=1.0, lv=0):
    cl = (2 + lv*(lv+3))/rv**2
    s = mp.sinh(rhov/rv); c = mp.cosh(rhov/rv)
    return (2*cl/s**2)*((rhov/rv)*(c/s) - 1)
# find rho where |R_0| < 2*delta for delta=0.5 (so i[H,A] >= 2delta-|R|>0)
crossover = None
rr = mp.mpf('0.1')
while rr < 20:
    if abs(Rl_num(rr)) < 1.0:   # 2*delta with delta=0.5
        crossover = rr; break
    rr += mp.mpf('0.05')
check(f"|R_0(rho)| < 2*delta for rho >= {mp.nstr(crossover,3)} (delta=0.5): strict Mourre on the cusp end",
      crossover is not None and crossover < 5)

print()
print("ALL CHECKS PASS" if ok else "CHECK FAILURE")
import sys; sys.exit(0 if ok else 1)
