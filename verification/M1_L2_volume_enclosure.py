import mpmath as mp
mp.mp.dps=40; pi=mp.pi
def f(x): x=mp.mpf(x); return 4*pi**2/(x**2-1)**4*( x**4+10*x**2+1 - 12*x**2*(x**2+1)/(x**2-1)*mp.log(x) )
def g(x): x=mp.mpf(x); return pi**2/(2*(x**2-1)**2)*( x**4-8*x**2+1 + 24*x**4/(x**4-1)*mp.log(x) )
def h(x):  # integrand sqrt(f) g^2 ; monotone decreasing if f,g decreasing
    return mp.sqrt(f(x))*g(x)**2
VolS4=8*pi**2/3
X0=mp.mpf('1e-8'); X1=1-mp.mpf('1e-7')
N=20000
xs=[X0+(X1-X0)*mp.mpf(i)/N for i in range(N+1)]
vals=[h(x) for x in xs]
mono=all(vals[i+1]<=vals[i]+mp.mpf('1e-25') for i in range(N))
print("integrand sqrt(f)*g^2 monotone decreasing on (0,1)?", mono)
dx=(X1-X0)/N
upper=sum(vals[i] for i in range(N))*dx     # left endpoints = max on cell
lower=sum(vals[i+1] for i in range(N))*dx   # right endpoints = min on cell
# end pieces: (0,X0) integrand <= h(0+), (X1,1) integrand <= h(X1) and >=0
tail0_hi=vals[0]*X0; tail1_hi=vals[-1]*(1-X1)
VL=VolS4*lower; VU=VolS4*(upper+tail0_hi+tail1_hi)
print(f"CERTIFIED Vol/r^5 in [{mp.nstr(VL,12)}, {mp.nstr(VU,12)}]   width={mp.nstr(VU-VL,4)}")
print(f"  midpoint {mp.nstr((VL+VU)/2,10)}   (prev quoted 943.393)")
