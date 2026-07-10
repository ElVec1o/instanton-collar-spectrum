#!/usr/bin/env python3
"""Elementary / hyperbolic-distance closed form of the Lemma 4.1 cross-term.

Claim (Lemma 4.1, elementary form). With a = s1^2, b = s2^2, c = R^2,
Sigma = a + b + c, D = Sigma^2 - 4ab = (R^2+(s1+s2)^2)(R^2+(s1-s2)^2):

  I := int_0^1 t(1-t) (2X - t(1-t)c) / X^2 dt,   X = c t(1-t) + a(1-t) + b t
     = Sigma/D - (4ab / D^{3/2}) * arccosh( Sigma / (2 s1 s2) ).

Equivalently, with cosh d = Sigma/(2 s1 s2)  (d = hyperbolic distance between
(x1, s1) and (x2, s2) in the upper-half-space model of H^5):

  <dA1/ds1, dA2/ds2>_{L^2} = 48 pi^2 s1 s2 I = 24 pi^2 (cosh d sinh d - d)/sinh^3 d.

This script verifies:
  (0) the pointwise integrand identity  t(1-t)(2X - t(1-t)c)/X^2 = (1/c)(1 - L^2/X^2),
      L = a(1-t) + bt   (symbolic, exact);
  (1) the polynomial identity N = Sigma*D used in the derivation (symbolic, exact);
  (2) closed form vs 50-digit adaptive quadrature of the 1-D integral, 60-point grid;
  (3) closed form (x 48 pi^2 s1 s2) vs the Rust 4-D ground-truth values
      recorded in lemma_ct_rust/RESULTS.md (machine precision, ~2.5 ulp);
  (4) the limits: G(d) -> 2/3 as d -> 0   [24 pi^2 * 2/3 = 16 pi^2 = diagonal norm],
      G(d) ~ 2 e^{-d} as d -> infinity    [recovers 48 pi^2 s1 s2 / R^2],
      on-axis limit s2 -> 0: c*I -> 1/(1+u^2), u = s1/R;
  (5) the symmetric next-order expansion F(eps,eps) = 1 - 2 eps^2 + 8(1+log eps) eps^4 + ...
      (the paper's log term: log = -d + log-free part).

Exit code 0 and final PASS line iff all checks pass.
"""
import sympy as sp
import mpmath as mp
import sys

ok = True

def check(name, cond):
    global ok
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    ok = ok and bool(cond)

# ---------- (0) pointwise identity, symbolic ----------
t, a, b, c = sp.symbols('t a b c', positive=True)
X = c*t*(1-t) + a*(1-t) + b*t
L = a*(1-t) + b*t
P = t*(1-t)
res = sp.simplify(P*(2*X - P*c)/X**2 - (1 - L**2/X**2)/c)
print("(0) pointwise integrand identity")
check("t(1-t)(2X-t(1-t)c)/X^2 == (1/c)(1 - L^2/X^2)", res == 0)

# ---------- (1) polynomial identity N = Sigma*D ----------
Sig = a + b + c
D = sp.expand(Sig**2 - 4*a*b)
N = (a**3 - a**2*b + 3*a**2*c - a*b**2 + 2*a*b*c + 3*a*c**2
     + b**3 + 3*b**2*c + 3*b*c**2 + c**3)
print("(1) numerator identity")
check("N == Sigma * D", sp.expand(N - Sig*D) == 0)

# ---------- closed form ----------
mp.mp.dps = 50

def I_closed(s1, s2, R):
    av, bv, cv = mp.mpf(s1)**2, mp.mpf(s2)**2, mp.mpf(R)**2
    S = av + bv + cv
    Dv = S**2 - 4*av*bv
    return S/Dv - 4*av*bv/Dv**mp.mpf(1.5) * mp.acosh(S/(2*mp.mpf(s1)*mp.mpf(s2)))

def I_quad(s1, s2, R):
    av, bv, cv = mp.mpf(s1)**2, mp.mpf(s2)**2, mp.mpf(R)**2
    def g(tv):
        Xv = cv*tv*(1-tv) + av*(1-tv) + bv*tv
        Pv = tv*(1-tv)
        return Pv*(2*Xv - Pv*cv)/Xv**2
    return mp.quad(g, [0, 1])

def G(d):
    d = mp.mpf(d)
    return (mp.cosh(d)*mp.sinh(d) - d)/mp.sinh(d)**3

# ---------- (2) grid vs quadrature ----------
print("(2) closed form vs 50-digit quadrature, 60-point grid")
worst = mp.mpf(0)
for s1 in [0.01, 0.05, 0.3, 1.0, 2.0]:
    for s2 in [0.02, 0.5, 1.0, 3.0]:
        for R in [0.5, 1.0, 10.0]:
            e = abs(I_closed(s1, s2, R)/I_quad(s1, s2, R) - 1)
            worst = max(worst, e)
check(f"max rel err = {mp.nstr(worst, 3)} < 1e-45", worst < mp.mpf('1e-45'))

# ---------- (3) vs Rust 4-D ground truth (RESULTS.md rows) ----------
print("(3) 48 pi^2 s1 s2 I vs Rust 4-D integral values (RESULTS.md)")
rows = [
    (0.05, 0.05, 1,  '1.17831447946537e0'),
    (0.05, 0.1,  5,  '9.47007380292315e-2'),
    (0.05, 0.5,  10, '1.18136742582047e-1'),
    (0.05, 2,    10, '4.55499235638643e-1'),
    (0.1,  0.5,  1,  '1.85328596930434e1'),
    (0.5,  0.5,  50, '4.73646126502408e-2'),
    (0.5,  1,    1,  '8.92273557362241e1'),
]
worst3 = mp.mpf(0)
for s1, s2, R, val in rows:
    full = 48*mp.pi**2*mp.mpf(s1)*mp.mpf(s2)*I_closed(s1, s2, R)
    worst3 = max(worst3, abs(full/mp.mpf(val) - 1))
check(f"max rel err = {mp.nstr(worst3, 3)} < 5e-15 (table stored at ~16 digits)",
      worst3 < mp.mpf('5e-15'))

# ---------- (4) limits ----------
print("(4) limits and hyperbolic form")
# G(d) -> 2/3 (=> 24 pi^2 G -> 16 pi^2, the diagonal norm)
check("G(1e-8) = 2/3 to 15 digits", abs(G('1e-8') - mp.mpf(2)/3) < mp.mpf('1e-15'))
# hyperbolic form == closed form
s1, s2, R = mp.mpf('0.3'), mp.mpf('0.7'), mp.mpf('2.0')
S = s1**2 + s2**2 + R**2
d = mp.acosh(S/(2*s1*s2))
lhs = 48*mp.pi**2*s1*s2*I_closed(s1, s2, R)
rhs = 24*mp.pi**2*G(d)
check("48 pi^2 s1 s2 I == 24 pi^2 G(d)", abs(lhs/rhs - 1) < mp.mpf('1e-45'))
# large-d: G ~ 2 e^{-d}
check("G(30)/ (2 e^{-30}) = 1 + O(e^{-2d})", abs(G(30)/(2*mp.e**-30) - 1) < mp.mpf('1e-20'))
# on-axis: c*I -> 1/(1+u^2) as s2 -> 0
u = mp.mpf('0.3')
val = I_closed(u, mp.mpf('1e-12'), 1)
check("s2->0: I -> 1/(1+u^2) (R=1)", abs(val - 1/(1+u**2)) < mp.mpf('1e-20'))

# ---------- (5) symmetric next-order expansion ----------
print("(5) symmetric expansion F(eps,eps) = 1 - 2 eps^2 + 8(1+log eps) eps^4 + O(eps^6 log eps)")
worst5 = mp.mpf(0)
for ev in ['1e-3', '1e-4', '1e-5']:
    e_ = mp.mpf(ev)
    F = I_closed(e_, e_, 1)
    coeff = (F - (1 - 2*e_**2))/e_**4
    target = 8*(1 + mp.log(e_))
    worst5 = max(worst5, abs(coeff/target - 1))
check(f"coefficient matches 8(1+log eps), max rel dev = {mp.nstr(worst5, 3)} < 1e-4",
      worst5 < mp.mpf('1e-4'))

print()
print("ALL CHECKS PASS" if ok else "CHECK FAILURE")
sys.exit(0 if ok else 1)
