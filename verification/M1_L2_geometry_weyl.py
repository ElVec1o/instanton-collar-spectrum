import mpmath as mp
mp.mp.dps=40
pi=mp.pi
def f(x): x=mp.mpf(x); return 4*pi**2/(x**2-1)**4*( x**4+10*x**2+1 - 12*x**2*(x**2+1)/(x**2-1)*mp.log(x) )
def g(x): x=mp.mpf(x); return pi**2/(2*(x**2-1)**2)*( x**4-8*x**2+1 + 24*x**4/(x**4-1)*mp.log(x) )
VolS4=8*pi**2/3                       # vol of unit round S^4
top=1-mp.mpf('1e-10'); bot=mp.mpf('1e-12')
# total 5-volume
I=mp.quad(lambda x: mp.sqrt(f(x))*g(x)**2, [bot,mp.mpf('0.3'),mp.mpf('0.7'),top])
Vol=VolS4*I
# boundary S^4 (x=0): induced metric g(0)*unit -> 4-vol = g(0)^2 * VolS4
g0=g(bot); Area=g0**2*VolS4
# diameter (radial length center x=1 to boundary x=0)
L=mp.quad(lambda x: mp.sqrt(f(x)), [bot,mp.mpf('0.3'),mp.mpf('0.7'),top])
print("Vol(M)/R^5      =",mp.nstr(Vol,14), " = ",mp.nstr(Vol/pi**4,12),"* pi^4")
print("g(0)=B^2 bdry   =",mp.nstr(g0,14)," (exact pi^2/2 =",mp.nstr(pi**2/2,14),")")
print("Area(dM)/R^4    =",mp.nstr(Area,14)," = ",mp.nstr(Area/pi**6,12),"* pi^6  (exact 2pi^6/3=",mp.nstr(2*pi**6/3,12),")")
print("radial length L =",mp.nstr(L,14))
# Weyl constants: N(lam) ~ C5 Vol lam^{5/2} - (1/4) C4 Area lam^2
w5=8*pi**2/15    # vol unit 5-ball
w4=pi**2/2       # vol unit 4-ball
C5=w5/(2*pi)**5; C4=w4/(2*pi)**4
print("\nWeyl: N(lam) = A*lam^{5/2} -/+ B*lam^2 + ...")
print("  A = w5/(2pi)^5 * Vol =",mp.nstr(C5*Vol,12))
print("  B = (1/4) w4/(2pi)^4 * Area =",mp.nstr(mp.mpf(1)/4*C4*Area,12))
# is Vol a clean multiple of pi^k? try identify
print("\nID Vol/pi^4:",mp.identify(Vol/pi**4, ['log(2)','pi']))
print("ID Vol:",mp.identify(Vol,['pi']))
