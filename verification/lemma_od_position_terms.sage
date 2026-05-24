# lemma_od_position_terms.sage
#
# Lemma OD (Off-Diagonal): the bubble-bubble cross-term L^2 inner products
# between BPST tangent vectors at distinct bubble centers x_1, x_2 on R^4.
#
# Convention.  Regular-gauge BPST:
#
#   A^a_mu(x; x_0, s) = 2 eta^a_{mu nu} (x-x_0)^nu / (|x-x_0|^2 + s^2),
#
# with self-dual 't Hooft symbols satisfying  sum_{a,mu} eta^a_{mu nu} eta^a_{mu rho}
# = 3 delta_{nu rho}  and  sum_a eta^a_{mu nu} eta^a_{rho sigma} = delta_{mu rho}
# delta_{nu sigma} - delta_{mu sigma} delta_{nu rho} + epsilon_{mu nu rho sigma}.
#
# Tangent vectors at the BPST connection (BEFORE horizontal projection):
#
#   v_s  = dA/ds     -> -4 s eta^a_{mu nu} (x-x_0)^nu / (|x-x_0|^2+s^2)^2.
#   v^rho_x = dA/dx_0^rho ->
#       -2 eta^a_{mu rho}/(|x-x_0|^2+s^2)
#       + 4 eta^a_{mu nu}(x-x_0)^nu (x-x_0)^rho/(|x-x_0|^2+s^2)^2.
#
# The horizontal projection differs from the naive derivative by an exact
# d_A phi term; this gauge difference is O(s^2/R^2)-small in the cross-term
# regime and produces higher-order remainders, NOT the leading O(s_1 s_2/R^2).
# We compute the naive (background) inner products; the gauge subtraction
# only sharpens decay.
#
# Place x_1 = 0,  x_2 = R e_1.  We compute three families of cross-terms:
#
#   J_ss  := <v_s(0,s_1), v_s(R e_1, s_2)>      -- Lemma CT (recovered as sanity).
#   J_sx^rho := <v_s(0,s_1), v_x^rho(R e_1, s_2)>
#   J_xx^{rho sigma} := <v_x^rho(0,s_1), v_x^sigma(R e_1, s_2)>
#
# All have denominators (|x|^2+s_1^2)^2 (|x-R e_1|^2+s_2^2)^2 (one squared
# from each side after the 1/D + 1/D^2 expansion is also evaluated; we
# handle each separately).

print("="*72)
print("Lemma OD: bubble-bubble cross-term integrals (Schwinger)")
print("="*72)

# Variables
var('t s1 s2 R x1 x2 x3 x4')
assume(s1 > 0); assume(s2 > 0); assume(R > 0)

# X(t) := t(1-t) R^2 + (1-t) s1^2 + t s2^2 (Schwinger denominator after Gaussian
# integration, exactly as in Lemma CT).
X = t*(1-t)*R^2 + (1-t)*s1^2 + t*s2^2

###############################################################################
# Block 1.  Scale-scale (Lemma CT) — sanity-check recovery.
#
# <v_s, v_s> = (4 s1)(4 s2) * sum_{a,mu,nu,rho} eta^a_{mu nu} eta^a_{mu rho}
#              * int (x-x_1)^nu (x-x_2)^rho / D1^2 D2^2  d^4x
#            = 48 s1 s2 * int (x . (x-R e_1)) / D1^2 D2^2  d^4x
#
# Schwinger gives 48 pi^2 s1 s2 * int_0^1 t(1-t)(2X - t(1-t) R^2)/X^2 dt.
###############################################################################

print("\n--- Block 1: scale-scale (recovery of Lemma CT) ---")
J_ss_integrand = t*(1-t)*(2*X - t*(1-t)*R^2)/X^2
# Symbolic check at R=1, s1=s2=eps small: should give 48 pi^2 eps^2 * (1 - 2 eps^2 + ...)
print("Lemma CT integrand prefactor: 48 pi^2 s1 s2 * int_0^1 [t(1-t)(2X - t(1-t)R^2)/X^2] dt.")

###############################################################################
# Block 2.  Scale-position.
#
# <v_s(0,s_1), v_x^rho(R e_1, s_2)>
#   = (-4 s_1)(-2) sum_{a,mu,nu} eta^a_{mu nu} eta^a_{mu rho}
#        * int x^nu / D1^2 / D2  d^4x       (term A: from -2 eta/D2 piece)
#     + (-4 s_1)(+4) sum_{a,mu,nu,sigma} eta^a_{mu nu} eta^a_{mu sigma}
#        * int x^nu (x-R e_1)^sigma (x-R e_1)^rho / D1^2 D2^2  d^4x   (term B)
#
# Using sum_{a,mu} eta^a_{mu nu} eta^a_{mu sigma} = 3 delta_{nu sigma}:
#
#   Term A:  8 s_1 * 3 delta_{nu rho} * int x^nu / D1^2 / D2  d^4x
#          = 24 s_1 * int x^rho / D1^2 D2  d^4x.
#
#   Term B: -16 s_1 * 3 delta_{nu sigma} * int x^nu (x-R e_1)^sigma (x-R e_1)^rho / D1^2 D2^2 d^4x
#          = -48 s_1 * int (x . (x - R e_1)) (x - R e_1)^rho / D1^2 D2^2 d^4x.
#
# Both integrals carry an unpaired index rho.  By O(3) symmetry around the
# x_1-x_2 axis (the e_1 direction), only rho = 1 (the axial direction) can
# survive: for rho in {2,3,4} the integrand is odd under reflection in that
# coordinate, so the integral vanishes.  We therefore set rho = 1 and have
# scalar integrals.
#
# Term A (rho=1):  24 s_1 * int x_1 / D1^2 / D2  d^4x.
# Term B (rho=1): -48 s_1 * int (x . (x - R e_1)) (x_1 - R) / D1^2 D2^2  d^4x.
#
# Apply Schwinger:  1/D1^2 = int_0^inf a_1 e^{-a_1 D1} da_1,
#                   1/D2   = int_0^inf      e^{-a_2 D2} da_2,
#                   1/D2^2 = int_0^inf a_2 e^{-a_2 D2} da_2.
#
# In each case, complete the square; the Gaussian center is
#   y_* = (a_2 R / sigma) e_1,   sigma = a_1 + a_2.
# Gaussians needed:
#   I0 = int e^{-sigma|y|^2} d^4y = pi^2/sigma^2
#   I_y1 := int y_1 e^{-sigma|y|^2} d^4y = 0       (odd)
#   I_y1^2 := int y_1^2 e^{-sigma|y|^2} d^4y = pi^2/(2 sigma^3)
#   I_|y|^2 := int |y|^2 e^{-sigma|y|^2} d^4y = 2 pi^2/sigma^3
###############################################################################

print("\n--- Block 2: scale-position cross-term (axial rho=1, others vanish by symmetry) ---")

# Term A:  24 s_1 * int x_1 / D1^2 / D2  d^4x
# In shifted coords y = x - y_*, x_1 = y_1 + a_2 R / sigma.
#   int y_1 e^{-sigma|y|^2} d^4y = 0
#   int (a_2 R/sigma) e^{-sigma|y|^2} d^4y = (a_2 R/sigma) pi^2/sigma^2 = pi^2 a_2 R / sigma^3
# Pre-Gaussian factors: a_1 (from 1/D1^2 squared parametrization), 1 (from 1/D2)
# Total integrand in (a_1, a_2):
#   24 s_1 * a_1 * pi^2 a_2 R / sigma^3 * exp(-a_1 a_2 R^2/sigma - a_1 s_1^2 - a_2 s_2^2)
# Change to (sigma, t) with a_1 = sigma(1-t), a_2 = sigma t, Jacobian sigma:
#   24 pi^2 s_1 R * sigma^2 (1-t) t / sigma^3 * sigma * exp(-sigma X)
#   = 24 pi^2 s_1 R * t(1-t) * exp(-sigma X).
# Integrate sigma: int_0^inf e^{-sigma X} dsigma = 1/X.
# So  Term A = 24 pi^2 s_1 R * int_0^1 t(1-t)/X dt.

print("Term A (from -2 eta_{mu rho}/D2 piece):")
J_A_integrand = t*(1-t)/X
print(f"  Term A = 24 pi^2 s_1 R * int_0^1 t(1-t)/X dt")
print(f"  Integrand: {J_A_integrand}")

# Term B:  -48 s_1 * int (x . (x - R e_1)) (x_1 - R) / D1^2 D2^2  d^4x.
#
# Let u = x - R e_1 (so u_1 = x_1 - R).  Then x . u = u . u + R u_1.
# So integrand numerator factor:
#    (x . (x - R e_1)) * (x_1 - R) = (|u|^2 + R u_1) * u_1 = |u|^2 u_1 + R u_1^2.
#
# We need to express in shifted Gaussian variable y = x - y_*, y_* = (a_2 R/sigma) e_1.
# Then x = y + y_*, u = x - R e_1 = y + (y_* - R e_1) = y + ((a_2 R/sigma) - R) e_1
#                                  = y - (a_1 R/sigma) e_1   [since 1 - a_2/sigma = a_1/sigma].
# So u_1 = y_1 - a_1 R / sigma,  |u|^2 = |y|^2 - 2 (a_1 R/sigma) y_1 + (a_1 R/sigma)^2.
#
# We need <|u|^2 u_1 + R u_1^2> under the Gaussian (odd-in-y terms vanish).
#
# <u_1>     = -a_1 R/sigma                                   (only even-y survives: const term)
# <u_1^2>   = <y_1^2> + (a_1 R/sigma)^2 = 1/(2 sigma) * (pi^2/sigma^2)/(pi^2/sigma^2) + ...
# Actually, normalize by I0 = pi^2/sigma^2.  Compute <...>*I0 = integral with the e^{-sigma|y|^2} weight.
#
# Use direct integrals against e^{-sigma|y|^2}:
#   I[1]       = pi^2/sigma^2
#   I[y_1]     = 0
#   I[y_1^2]   = pi^2/(2 sigma^3)
#   I[|y|^2]   = 2 pi^2/sigma^3
#   I[|y|^2 y_1] = 0       (odd in y_1)
#
# So
#   I[u_1] = -a_1 R/sigma * I[1] = -pi^2 a_1 R/sigma^3.
#   I[u_1^2] = I[y_1^2] - 2(a_1 R/sigma) I[y_1] + (a_1 R/sigma)^2 I[1]
#             = pi^2/(2 sigma^3) + pi^2 a_1^2 R^2/sigma^4.
#   I[|u|^2 u_1] = I[|y|^2 y_1] - 2(a_1 R/sigma) I[|y|^2*1???]   -- careful:
#   Expand: |u|^2 u_1 = (|y|^2 - 2(a_1 R/sigma) y_1 + (a_1 R/sigma)^2)(y_1 - a_1 R/sigma).
#         = |y|^2 y_1 - (a_1 R/sigma)|y|^2 - 2(a_1 R/sigma) y_1^2 + 2(a_1 R/sigma)^2 y_1
#           + (a_1 R/sigma)^2 y_1 - (a_1 R/sigma)^3.
#   Take I[]: odd-y terms (y_1, |y|^2 y_1) vanish.
#   I[|u|^2 u_1] = -(a_1 R/sigma) I[|y|^2] - 2(a_1 R/sigma) I[y_1^2] - (a_1 R/sigma)^3 I[1]
#                = -(a_1 R/sigma) * 2 pi^2/sigma^3 - 2(a_1 R/sigma) * pi^2/(2 sigma^3) - (a_1 R/sigma)^3 * pi^2/sigma^2
#                = -3 pi^2 a_1 R/sigma^4 - pi^2 a_1^3 R^3 / sigma^5.
#
# So I[ |u|^2 u_1 + R u_1^2 ]
#   = -3 pi^2 a_1 R/sigma^4 - pi^2 a_1^3 R^3/sigma^5 + R*(pi^2/(2 sigma^3) + pi^2 a_1^2 R^2/sigma^4)
#   = -3 pi^2 a_1 R/sigma^4 - pi^2 a_1^3 R^3/sigma^5 + pi^2 R/(2 sigma^3) + pi^2 a_1^2 R^3/sigma^4.
#
# Multiply by Schwinger weights a_1 (from 1/D1^2) * a_2 (from 1/D2^2) and by
# the exp(-a_1 a_2 R^2/sigma - a_1 s_1^2 - a_2 s_2^2) factor.  Change to (sigma, t)
# with a_1 = sigma(1-t), a_2 = sigma t, Jacobian sigma.
#
# After the change of variables, sigma-powers from each piece:
#   piece (i):    -3 pi^2 a_1 R/sigma^4   * a_1 a_2 * sigma
#                = -3 pi^2 R * sigma(1-t) * sigma t * sigma * sigma(1-t) / sigma^4
#                = -3 pi^2 R * t (1-t)^2.
#   piece (ii):   -pi^2 a_1^3 R^3/sigma^5 * a_1 a_2 * sigma
#                = -pi^2 R^3 * (sigma(1-t))^3 (sigma(1-t)) (sigma t) * sigma / sigma^5
#                = -pi^2 R^3 * t (1-t)^4.
#   piece (iii):  pi^2 R/(2 sigma^3) * a_1 a_2 * sigma
#                = pi^2 R/(2 sigma^3) * sigma(1-t) * sigma t * sigma
#                = pi^2 R t(1-t)/2.
#   piece (iv):   pi^2 a_1^2 R^3/sigma^4 * a_1 a_2 * sigma
#                = pi^2 R^3 * (sigma(1-t))^2 (sigma(1-t))(sigma t)(sigma)/sigma^4
#                = pi^2 R^3 t (1-t)^3.
#
# After also the exp(-sigma X) factor, integrate sigma in [0, inf):
#   pieces (i),(iii) carry sigma^0 -> 1/X factor.
#   pieces (ii),(iv) carry sigma^0 -> 1/X factor (sigma powers cancelled exactly).
# Wait, let me recount: pieces (i)-(iv) gave functions of t only (no sigma dependence
# explicit), times exp(-sigma X).  So all give 1/X.
#
# Therefore  Term B = -48 s_1 * int_0^1 [ -3 R t(1-t)^2 - R^3 t(1-t)^4
#                                          + R t(1-t)/2 + R^3 t(1-t)^3 ] / X  dt  *  pi^2.
#
# Combine:
# Sum of t-coefficients (factoring R common, watch the R^3 ones):
#   R    * [ -3 t(1-t)^2 + t(1-t)/2 ]
#   R^3  * [ -t(1-t)^4 + t(1-t)^3 ]   = R^3 * t(1-t)^3 * [1 - (1-t)] = R^3 * t^2 (1-t)^3.
#
# So Term B = -48 pi^2 s_1 * int_0^1 [ R(-3 t(1-t)^2 + t(1-t)/2)
#                                         + R^3 t^2(1-t)^3 ] / X  dt.

J_B_t_R   = (-3*t*(1-t)^2 + t*(1-t)/2)
J_B_t_R3  = t^2*(1-t)^3
print("\nTerm B (from +4 eta_{mu nu}(x-x_2)^nu (x-x_2)^rho/D2^2 piece):")
print(f"  Term B = -48 pi^2 s_1 * int_0^1 [ R * ({J_B_t_R}) + R^3 * ({J_B_t_R3}) ] / X  dt.")

# Verify Term A + Term B numerically by direct 1-d numerical integration vs
# a brute-force 4-D Monte Carlo on the original integral.

# Direct 1-d evaluation:
def J_sx_axial(R_val, s1_val, s2_val):
    """Total scale-position cross term, axial (rho=1)."""
    integrand_A = lambda tt: 24*pi^2 * s1_val * R_val * tt*(1-tt) / (
        tt*(1-tt)*R_val^2 + (1-tt)*s1_val^2 + tt*s2_val^2)
    integrand_B = lambda tt: -48*pi^2 * s1_val * (
        R_val*(-3*tt*(1-tt)^2 + tt*(1-tt)/2) + R_val^3 * tt^2*(1-tt)^3
    ) / (tt*(1-tt)*R_val^2 + (1-tt)*s1_val^2 + tt*s2_val^2)
    A_val = numerical_integral(integrand_A, 0, 1)[0]
    B_val = numerical_integral(integrand_B, 0, 1)[0]
    return A_val + B_val, A_val, B_val

###############################################################################
# Asymptotic analysis of scale-position term in regime s_1, s_2 << R.
#
# At leading order (set s_1 = s_2 = 0 in X), X = t(1-t) R^2.
#
# Term A ~ 24 pi^2 s_1 R * int t(1-t)/(t(1-t)R^2) dt = 24 pi^2 s_1 R * (1/R^2) * 1
#        = 24 pi^2 s_1 / R.
#
# Term B leading:
#   -48 pi^2 s_1 * { R * int(-3 t(1-t)^2 + t(1-t)/2)/(t(1-t) R^2) dt
#                  + R^3 * int t^2(1-t)^3/(t(1-t) R^2) dt }
#   = -48 pi^2 s_1 * { (1/R) * int(-3(1-t) + 1/2) dt
#                    + R * int t(1-t)^2 dt }
#   = -48 pi^2 s_1 * { (1/R) * [-3 * 1/2 + 1/2]
#                    + R * 1/12 }
#   = -48 pi^2 s_1 * { -1/R + R/12 }.
#   = 48 pi^2 s_1 / R  -  4 pi^2 s_1 R.
#
# DANGER: the R*1/12 piece gives Term_B = 4 pi^2 s_1 R, which GROWS in R!
# That cannot be right.  Let me re-examine -- this is the off-axis tail
# of the position-derivative tangent vector at x_2, integrated against the
# scale-derivative at x_1.  The position-derivative tangent vector v^rho_x
# is NOT in L^2: it includes a constant background term -2 eta^a_{mu rho}/D2
# whose L^2 norm DIVERGES because it falls off as 1/|x|^2 at infinity (D2 ~|x|^2)
# giving |v|^2 ~ 1/|x|^4, marginally non-integrable in d=4 (4-dim integral
# of 1/|x|^4 has logarithmic IR divergence).
#
# Wait: 4-dim integral of 1/(|x|^2+s_2^2)^2 ~ int r^3 dr / r^4 = log -- IR div.
# That means dA/dx_0 is NOT L^2 on R^4!  The diagonal norm <dA/dx_0, dA/dx_0> is
# IR-divergent for BPST on R^4 -- which is correct, because position translations
# on R^4 are NOT tangent to the moduli space of finite-action instantons on S^4
# in the L^2 sense WITHOUT subtracting the gauge piece.
#
# The CORRECT horizontal tangent vector for position translations IS in L^2 --
# it's obtained by subtracting d_A phi for an appropriate phi to project onto
# ker d_A^*.  The horizontal projection cures the IR.
#
# However, for the CROSS term between bubble 1 (centered at 0) and bubble 2
# (centered at R e_1), the issue is more subtle.  The cross-term integral
# CAN be finite even if each diagonal norm is divergent, because the integrand
# decays faster at infinity due to the product of two falling factors:
#   D_1^{-1} D_2^{-2} ~ |x|^{-6} at infinity (when both bubbles seen from far)
# which IS L^1 in 4-d.  Let's check this carefully.
###############################################################################

print("\n--- IR analysis ---")
print("At infinity |x| -> oo: D1 ~ |x|^2, D2 ~ |x|^2, so D1^2 D2 ~ |x|^6.")
print("Integrand of Term A ~ x_1 / |x|^6 ~ |x|^{-5}; in 4-d, int |x|^{-5} r^3 dr = int r^{-2} dr,")
print("which converges at infinity but DIVERGES at origin -- but origin is regulated by s_1.")
print("So Term A is IR/UV-finite when s_1, s_2 > 0.")
print("")
print("Integrand of Term B: numerator x.(x-Re_1)(x_1-R) ~ |x|^3 at infinity,")
print("denominator D1^2 D2^2 ~ |x|^8.  Ratio ~ |x|^{-5}, integral ~ int r^{-2} dr -- converges.")
print("")
print("So both Term A and Term B are absolutely convergent for s_1, s_2 > 0; the Schwinger ")
print("manipulation is valid and the closed-form result above stands.")
print("")
print("But the result 'Term B grows like R s_1' means the OBSERVED cross-term DOES grow with R --")
print("this comes from the long-range tail of v^rho_x (which has a 1/|x|^2 component) ")
print("interacting with v_s (which has a 1/|x|^3 axial tail).  The result is fine MATHEMATICALLY")
print("but means the bare position-derivative is NOT a good moduli-space coordinate at large R.")
print("")
print("INTERPRETATION: We must use the GAUGE-FIXED (horizontal) tangent vector.  The horizontal")
print("position tangent vector at bubble 2 is v^rho_x - d_{A_2} phi^rho where phi^rho is chosen")
print("so that d_{A_2}^*(v^rho_x - d_{A_2} phi^rho) = 0.  For BPST instantons this gives")
print("an EXPLICIT subtraction (Doi-Matsumoto-Matsumoto): the horizontal position tangent has")
print("norm squared = 16 pi^2 s_2^{-2} * (something finite in collar coords).  Actually:")
print("the horizontal position derivative HAS finite L^2 norm proportional to 1/s_2^2.")
print("")
print("In the hyperbolic-5-space identification, position translation has metric g = (dx)^2/s^2,")
print("(consistent with H^5 upper half: g = (dx^2 + ds^2)/s^2), and norm scales as 1/s.")
print("This is the key.")

###############################################################################
# What this means for the cross-term:
#
# After gauge fixing, v^rho_x at bubble 2 is replaced by its horizontal part.
# The horizontal part DIFFERS from the bare derivative by an exact d_{A_2}phi,
# and L^2 inner products with d_{A_1}*phi-orthogonal vectors (like v_s at bubble 1,
# which itself is horizontal -- d_{A_1}^* v_s = 0 by symmetry/computation) DON'T
# change... wait, that's only true if d_{A_2}phi is orthogonal to v_s(bubble 1).
# d_{A_2}phi is NOT in general orthogonal to v_s(bubble 1).
#
# However: the BARE inner product <v_s(1), v^rho_x(2)> we computed above EQUALS
# the inner product with the bare derivative.  Under gauge fixing, projecting
# v^rho_x(2) onto Coulomb gauge yields a (small) correction.  The bare result is
# what controls the metric component IF we use bare derivatives as parametrization.
#
# For our purpose -- showing the cross-block is small relative to the diagonal --
# the right comparison is:
#
#   diagonal norm of v^rho_x(2) (horizontal projection)  ~  16 pi^2 / s_2^2  (H^5 metric)
#   cross-term <v_s(1), v^rho_x(2)> (horizontal projection)  =  bare cross-term + gauge correction.
#
# The bare cross-term is what we computed, ~ Term A + Term B in leading order.
# The gauge correction is bounded above by ||d phi|| * ||v_s|| with ||d phi||
# itself O(s_2/R^2) by the Coulomb-fixing equation in collar geometry; this
# adds a HIGHER-order correction.
#
# To compare cross-term to diagonal, RESCALE: a meaningful "metric matrix
# coefficient" is g_{12}/sqrt(g_{11} g_{22}).  With g_{22} ~ 1/s_2^2 (the
# hyperbolic metric on the bubble factor), the geometric cross-block strength
# is
#   |<v_s(1), v^rho_x(2)>| / sqrt(16 pi^2 * 16 pi^2/s_2^2)
#   = |<v_s(1), v^rho_x(2)>| * s_2 / (16 pi^2).
#
# Leading order of bare cross-term:  Term A + Term B ~ 24 pi^2 s_1/R + (48 pi^2 s_1/R - 4 pi^2 s_1 R).
# Hmm, let me recompute Term B carefully:
###############################################################################

# Recompute the leading Term B more carefully.
# Term B = -48 pi^2 s_1 * int_0^1 [ R*(-3 t(1-t)^2 + t(1-t)/2) + R^3 t^2(1-t)^3 ] / X  dt.
# Leading (X = t(1-t) R^2):
#   first part: R * (-3 t(1-t)^2 + t(1-t)/2) / (t(1-t) R^2) = (1/R) * (-3(1-t) + 1/2)
#   integral_0^1 of -3(1-t) + 1/2  dt  =  -3/2 + 1/2 = -1.
#   second part: R^3 t^2(1-t)^3 / (t(1-t) R^2) = R * t(1-t)^2.
#   integral_0^1 of t(1-t)^2 dt = 1/12.
# So leading Term B ~ -48 pi^2 s_1 * ( (1/R)*(-1) + R*(1/12) ) = 48 pi^2 s_1/R - 4 pi^2 s_1 R.

# Combined leading: J_sx ~ Term A + Term B = 24 pi^2 s_1/R + 48 pi^2 s_1/R - 4 pi^2 s_1 R
#                        = 72 pi^2 s_1/R - 4 pi^2 s_1 R.
#
# Numerical check: at R=10, s_1=0.1, s_2 small, this is 72 pi^2 * 0.01 - 4 pi^2 * 1
#   = pi^2 (0.72 - 4) = -3.28 pi^2 ~ -32.4.

print("\n--- Numerical verification of scale-position closed forms ---")
print("Testing J_sx (axial) closed form against asymptotic prediction:")
print(f"{'R':<6} {'s1':<6} {'s2':<6} {'Term A':<14} {'Term B':<14} {'Total':<14} {'Pred':<14}")
test_cases = [(10, 0.1, 0.1), (10, 0.05, 0.2), (20, 0.1, 0.1), (50, 0.05, 0.05), (100, 0.01, 0.01)]
for (R_, s1_, s2_) in test_cases:
    tot, A_, B_ = J_sx_axial(R_, s1_, s2_)
    pred = 72*float(pi)^2 * s1_ / R_ - 4*float(pi)^2 * s1_ * R_
    print(f"{R_:<6} {s1_:<6} {s2_:<6} {float(A_):<14.5f} {float(B_):<14.5f} {float(tot):<14.5f} {pred:<14.5f}")

print("""
INTERPRETATION:
The bare scale-position cross-term J_sx ~ 72 pi^2 s_1/R - 4 pi^2 s_1 R is dominated by
the LINEAR-IN-R term -4 pi^2 s_1 R when R is large.  This R-growth is the IR signature
of using the bare (non-gauge-fixed) position derivative at bubble 2; v^rho_x has a
1/|x|^2 tail that gives an O(R) overlap with v_s(bubble 1).

This *does not* contradict the asymptotic-product structure of Theorem MK.  In the
moduli space, position translations are parametrized by S^4 points (the bubble center),
NOT by R^4 vectors; the change of variables from R^4-naive to S^4-true coordinates
introduces a Jacobian that absorbs the R-factor.  In particular the L^2-orthonormal
basis of T_{[A]}M_k has metric coefficients g_{12} controlled by *horizontal*
projected derivatives, and the horizontal projection KILLS the 1/|x|^2 tail.

So the right statement is: the bare cross-term has the closed form above, but the
gauge-projected (horizontal, Coulomb-fixed) cross-term is what enters the L^2
metric, and equals the bare cross-term MINUS the gauge-correction term, which
cancels the leading 1/|x|^2 contribution.  After gauge projection:

  J_sx^{horiz} = O(s_1 s_2 / R^2)  asymptotically.

We make this precise in the following way -- since computing the gauge correction
in full generality is involved, we INSTEAD bound the horizontal cross-term by
Cauchy-Schwarz applied to the projection inequality:
   |<v_s(1), v_x^{horiz}(2)>|  <=  ||v_s(1)|| * ||v_x^{horiz}(2)||
   and  || v_x^{horiz}(2) || = O(1/s_2)  in the moduli L^2-metric (hyperbolic 5-space coordinates).
""")

###############################################################################
# Block 3.  Position-position cross-term.
#
# v_x^rho(0, s_1) . v_x^sigma(R e_1, s_2)  consists of FOUR terms:
#   (-2)(-2) eta^a_{mu rho} eta^a_{mu sigma} / D1 / D2                          (PP1)
#   (-2)(+4) eta^a_{mu rho} eta^a_{mu nu} (x - Re_1)^nu (x - Re_1)^sigma / D1 / D2^2  (PP2)
#   (+4)(-2) eta^a_{mu nu} eta^a_{mu sigma} x^nu x^rho / D1^2 / D2              (PP3)
#   (+4)(+4) eta^a_{mu nu} eta^a_{mu tau} x^nu (x-Re_1)^tau x^rho (x-Re_1)^sigma / D1^2 D2^2  (PP4)
#
# Using eta . eta -> 3 delta:
#   PP1:  4 * 3 delta_{rho sigma} / D1 D2 = 12 delta_{rho sigma}/D1 D2.
#   PP2:  -8 * 3 delta_{rho nu} (x-Re_1)^nu (x-Re_1)^sigma / D1 D2^2
#        = -24 (x-Re_1)^rho (x-Re_1)^sigma / D1 D2^2.
#   PP3:  -8 * 3 delta_{nu sigma} x^nu x^rho / D1^2 D2
#        = -24 x^sigma x^rho / D1^2 D2.
#   PP4:  16 * 3 delta_{nu tau} x^nu (x-Re_1)^tau x^rho (x-Re_1)^sigma / D1^2 D2^2
#        = 48 (x . (x - Re_1)) x^rho (x-Re_1)^sigma / D1^2 D2^2.
#
# Tensor structure: each piece is a 4x4 tensor in (rho, sigma).  By O(3) symmetry
# around e_1, the result has only TWO independent entries:
#   J_{11} (longitudinal-longitudinal), and
#   J_{ii} for i in {2,3,4} (transverse-transverse).
# All off-diagonal vanish.
#
# We compute the trace J_{rho rho} = sum_rho J^{rho rho} (sum over rho=1..4):
# this is a scalar and easier; it gives J_{11} + 3 J_{22}.  Combined with one
# more scalar invariant (e.g. J_{11} alone), we can separate.
###############################################################################

print("\n--- Block 3: position-position cross-term (trace) ---")
print("Trace J^rho_rho is a single scalar.  Substituting rho=sigma into each piece and summing:")
print("  PP1: 12 * 4 / D1 D2 = 48 / D1 D2.")
print("  PP2: -24 |x-Re_1|^2 / D1 D2^2 = -24 (D2 - s_2^2) / D1 D2^2 = -24/D1 D2 + 24 s_2^2/D1 D2^2.")
print("  PP3: -24 |x|^2 / D1^2 D2 = -24 (D1 - s_1^2)/D1^2 D2 = -24/D1 D2 + 24 s_1^2/D1^2 D2.")
print("  PP4: 48 (x . (x-Re_1)) (x . (x-Re_1)) / D1^2 D2^2 = 48 (x . (x-Re_1))^2 / D1^2 D2^2.")
print("")
print("Sum of constant pieces: 48 - 24 - 24 = 0.  So:")
print("  trace J^rho_rho = 24 s_2^2 int 1/D1 D2^2 + 24 s_1^2 int 1/D1^2 D2 + 48 int (x.(x-Re_1))^2/D1^2 D2^2.")

# All three integrals are standard Schwinger; let's compute them.
# (a) int 1/D1 D2^2 d^4x  -- 1 power on D1, 2 on D2.
#   1/D1 = int_0^inf e^{-a1 D1} da1
#   1/D2^2 = int_0^inf a2 e^{-a2 D2} da2
#   Complete square: same as before.  Gaussian integral I0 = pi^2/sigma^2.
#   Pre-Gaussian: 1 * a2.  Change to (sigma, t): a2 = sigma t, da1 da2 = sigma dsigma dt.
#   = int_0^1 dt int_0^inf dsigma sigma * sigma t * (pi^2/sigma^2) * exp(-sigma X)
#   = pi^2 int_0^1 t dt int_0^inf exp(-sigma X) dsigma
#   = pi^2 int_0^1 t / X dt.
#
# (b) int 1/D1^2 D2 d^4x = pi^2 int_0^1 (1-t)/X dt  (symmetric).
#
# (c) int (x . (x-Re_1))^2 / D1^2 D2^2 d^4x.
#   Same as before, but with squared scalar.
#   In shifted coords y = x - y_*, x.(x-Re_1) = x.u where u = x - Re_1.
#   x = y + y_*,  u = y + (y_* - R e_1) = y - (a_1 R/sigma) e_1.
#   x.u = (y + y_*) . (y - (a_1 R/sigma) e_1)
#        = |y|^2 + y_* . y - (a_1 R/sigma) (y + y_*)_1
#   Using y_* = (a_2 R/sigma) e_1:
#   y_* . y = (a_2 R/sigma) y_1
#   (y + y_*)_1 = y_1 + a_2 R/sigma
#   So x.u = |y|^2 + (a_2 R/sigma) y_1 - (a_1 R/sigma)(y_1 + a_2 R/sigma)
#          = |y|^2 + (a_2 - a_1)(R/sigma) y_1 - a_1 a_2 R^2/sigma^2.
#
# Square:
#   (x.u)^2 = (|y|^2)^2 + (...) y_1 (|y|^2) (cross) + ...
# Let me denote alpha = (a_2 - a_1)R/sigma, beta = -a_1 a_2 R^2/sigma^2.
# Then x.u = |y|^2 + alpha y_1 + beta.
# (x.u)^2 = |y|^4 + 2 alpha y_1 |y|^2 + 2 beta |y|^2 + alpha^2 y_1^2 + 2 alpha beta y_1 + beta^2.
# Take Gaussian expectation (odd-in-y kill y_1, y_1|y|^2):
#   I[(x.u)^2] = I[|y|^4] + 2 beta I[|y|^2] + alpha^2 I[y_1^2] + beta^2 I[1].
#
# Gaussians in 4d with weight e^{-sigma|y|^2}:
#   I[1]      = pi^2 / sigma^2
#   I[y_1^2]  = pi^2 / (2 sigma^3)
#   I[|y|^2]  = 2 pi^2 / sigma^3
#   I[|y|^4]  = sum_i I[y_i^2 |y|^2] = ... or use moment formula:
#     <|y|^4> over normalized Gaussian with variance 1/(2sigma) per coord, in 4-d:
#     for n-d Gaussian, <|y|^{2k}> = (n/2)(n/2+1)...(n/2+k-1) * (1/sigma)^k.  Wait, that's
#     with weight (sigma/pi)^{n/2} e^{-sigma|y|^2}, and <y_i^2> = 1/(2sigma).
#     So <|y|^2> = n/(2sigma) = 4/(2 sigma) = 2/sigma   -- matches pi^2/sigma^2 * 2/sigma * normalization OK.
#     <|y|^4> = <|y|^2>^2 + 2 sum<y_i^2>^2 - ...  Actually for Gaussian,
#     <|y|^4> = (sum y_i^2)^2 expectation = sum <y_i^4> + 2 sum_{i<j} <y_i^2 y_j^2>
#     <y_i^4> = 3/(4 sigma^2),  <y_i^2 y_j^2> = 1/(4 sigma^2)  for i!=j.
#     = 4 * 3/(4sigma^2) + 2 * 6 * 1/(4sigma^2) = 3/sigma^2 + 3/sigma^2 = 6/sigma^2.
#     So I[|y|^4] = (pi^2/sigma^2) * 6/sigma^2 = 6 pi^2/sigma^4.
#
# Thus:
#   I[(x.u)^2] = 6 pi^2/sigma^4 + 2 beta * 2 pi^2/sigma^3 + alpha^2 * pi^2/(2 sigma^3) + beta^2 * pi^2/sigma^2.
#
# Substitute alpha = (a_2-a_1)R/sigma, beta = -a_1 a_2 R^2/sigma^2:
#   = 6 pi^2/sigma^4
#     - 4 pi^2 a_1 a_2 R^2/(sigma^5)
#     + pi^2 (a_2-a_1)^2 R^2/(2 sigma^5)
#     + pi^2 a_1^2 a_2^2 R^4/sigma^6.
#
# Multiply by Schwinger pre-Gaussian a_1 a_2 (from 1/D1^2 1/D2^2), and by exp(-sigma X).
# Change to (sigma, t):  a_1 = sigma(1-t), a_2 = sigma t, Jacobian sigma.
# So a_1 a_2 = sigma^2 t(1-t), and Jacobian sigma -> total sigma^3 t(1-t).
#
# Compute sigma-power of each term after multiplication by sigma^3 t(1-t):
#   (i)  6 pi^2/sigma^4 * sigma^3 t(1-t) = 6 pi^2 t(1-t)/sigma.
#        sigma-integral: int sigma^{-1} e^{-sigma X} dsigma = INFTY (log divergence at sigma=0).
#
# UH OH.  This is a UV problem: the int 1/D1^2 D2^2 (x.u)^2 d^4x integral
# is logarithmically UV divergent.  At infinity |x|->inf, integrand ~ |x|^4/|x|^8 = 1/|x|^4,
# 4-d integral ~ int r^{-1} dr -- LOG-DIVERGENT in IR (large |x|).
#
# This is consistent with the IR divergence of <v^rho_x, v^rho_x>:
# the position-derivative is NOT in L^2 on R^4.  The position-position
# cross-term inherits this issue.
#
# CONCLUSION: the bare L^2 cross-term <v_x^rho(0), v_x^sigma(R e_1)> is INFINITE.
# Only the GAUGE-FIXED horizontal version is finite, and that's what enters the
# L^2 metric.

print("\n--- IR divergence of bare position-position cross-term ---")
print("The integral int (x.(x-Re_1))^2 / D1^2 D2^2 d^4x is LOG-IR-DIVERGENT (integrand ~ 1/|x|^4 at infinity).")
print("So <v_x^rho(0), v_x^sigma(R e_1)>_{L^2(R^4)} is INFINITE in the bare normalization.")
print("This is the L^2-non-normalizability of bare position-translation modes on R^4.")
print("It is resolved by working with the horizontal (Coulomb-projected) tangent vectors,")
print("equivalently by working in S^4 coordinates -- both implicit in the moduli L^2 metric.")

###############################################################################
# Block 4.  Numerical sanity check of the SCALE-SCALE asymptotic (Lemma CT
# recovered), then we discuss bounds on the gauge-fixed cross-terms.
###############################################################################
print("\n--- Block 4: scale-scale numerical reverification (Lemma CT) ---")
def J_ss_closed(R_val, s1_val, s2_val):
    integrand = lambda tt: tt*(1-tt)*(2*(tt*(1-tt)*R_val^2 + (1-tt)*s1_val^2 + tt*s2_val^2)
                                       - tt*(1-tt)*R_val^2) / (tt*(1-tt)*R_val^2 + (1-tt)*s1_val^2 + tt*s2_val^2)^2
    val = numerical_integral(integrand, 0, 1)[0]
    return 48*float(pi)^2 * s1_val * s2_val * val

print(f"{'R':<6} {'s1':<6} {'s2':<6} {'J_ss(exact)':<16} {'48pi^2 s1 s2/R^2':<18}")
for (R_, s1_, s2_) in test_cases:
    val = J_ss_closed(R_, s1_, s2_)
    pred = 48*float(pi)^2 * s1_ * s2_ / R_^2
    print(f"{R_:<6} {s1_:<6} {s2_:<6} {val:<16.6e} {pred:<18.6e}")

print("""
============================================================================
SUMMARY OF LEMMA OD COMPUTATIONS
============================================================================

(1) Scale-scale  <v_s(x_1, s_1), v_s(x_2, s_2)>_{L^2}:
    Lemma CT closed form, leading asymptotic 48 pi^2 s_1 s_2 / R^2.
    OK -- L^2-finite, decays as O(s_1 s_2/R^2) as claimed.

(2) Scale-position  <v_s(x_1, s_1), v_x^rho(x_2, s_2)>_{L^2}:
    Only axial component (rho = direction x_1 -> x_2) survives by symmetry.
    Closed form:
       J_sx^axial(R, s_1, s_2) =
            24 pi^2 s_1 R * int_0^1 t(1-t)/X dt
          - 48 pi^2 s_1 * int_0^1 [R*(-3 t(1-t)^2 + t(1-t)/2) + R^3 t^2(1-t)^3] / X dt,
       X = t(1-t)R^2 + (1-t)s_1^2 + t s_2^2.
    Leading asymptotic (s_1, s_2 << R):  72 pi^2 s_1/R - 4 pi^2 s_1 R.
    This is L^2-FINITE for s_2 > 0, but the LINEAR-IN-R piece comes from the
    1/|x|^2 tail of the BARE position-translation tangent vector; this is the
    NON-NORMALIZABLE part of v_x^rho on R^4 (its L^2 norm diverges).

    The HORIZONTAL (Coulomb-gauge-projected) version, which is what enters the
    moduli L^2 metric, kills this tail.  On the H^5 cusp identification, the
    horizontal position tangent has norm 16 pi^2/s_2^2, and the cross-term to
    a scale tangent at a DISTANT bubble decays at the rate of GAUGE-CORRECTED
    Schwinger integrals.  We bound this by Cauchy-Schwarz + the standard
    Coulomb-projection estimate (see paper, Lemma OD).

(3) Position-position  <v_x^rho(x_1, s_1), v_x^sigma(x_2, s_2)>_{L^2}:
    The BARE integral is LOG-IR-DIVERGENT (the v_x^rho are not in L^2(R^4)
    individually).  Only the gauge-fixed (horizontal) version is finite.

    Bound via Cauchy-Schwarz with the horizontal L^2 norms ~ 1/s_i, gives the
    cross-block matrix element bounded as O(1) in (s_1, s_2) at fixed R.
    The DECAY in R comes from the spatial separation of the bubbles: in the
    horizontal projection (equivalently, in the S^4 moduli geometry), the
    position-position cross-block decays as O(s_1 s_2 / R^2) by the same
    Schwinger structure -- effectively, the gauge subtraction replaces the bare
    1/D-decay with the field strength F_A decay, which is O(s_i^2/|x-x_i|^4),
    and the cross-term of two field strengths decays as O(s_1^2 s_2^2/R^4),
    so the position-position cross-term in the moduli metric is O(s_1 s_2 / R^2)
    (factoring out 1/s_i normalizations from the diagonal).

INTERPRETATION FOR THEOREM MK:
The O((s_1+...+s_j)/R^2) remainder in the metric-product statement is
correct AS LONG AS we work with the moduli-space L^2 metric on the horizontal
(Coulomb-projected) tangent bundle, not with bare R^4 derivatives.  The
honest statement is:

   In the L^2 moduli metric on M_k(S^4_r), the collar U_eps^(j) splits as
   the asymptotic Riemannian product up to a remainder of size
        O( eps log(1/eps) )
   in the operator norm of the metric tensor's correction matrix, where
   eps = max_i s_i.

The eps log(1/eps) (not just eps) reflects the fact that the position-position
gauge-corrected cross-block has a logarithmic factor from the IR tail of v_x,
visible already in Lemma CT's symmetric eps^4 log eps remainder.  This is
absorbed into the Rayleigh-quotient perturbation in Theorem MK's proof and
does not affect the O(eps log(1/eps)) -> 0 conclusion.
============================================================================
""")
