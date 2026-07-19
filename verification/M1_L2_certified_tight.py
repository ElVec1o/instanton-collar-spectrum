import mpmath as mp, math, sys, time
mp.mp.dps=30; pi=mp.pi
def f(x): x=mp.mpf(x); return 4*pi**2/(x**2-1)**4*( x**4+10*x**2+1 - 12*x**2*(x**2+1)/(x**2-1)*mp.log(x) )
def g(x): x=mp.mpf(x); return pi**2/(2*(x**2-1)**2)*( x**4-8*x**2+1 + 24*x**4/(x**4-1)*mp.log(x) )
FLIM=mp.mpf('0.4')*pi**2
def fsafe(x):
    x=mp.mpf(x)
    if x>1-mp.mpf('1e-7'): return FLIM
    if x<mp.mpf('1e-9'):   return 4*pi**2
    return f(x)
H_=mp.mpf('1e-9')
def gs(x):  x=mp.mpf(x); return (g(x+H_)-g(x-H_))/(2*H_)/mp.sqrt(fsafe(x))
def gss(x): x=mp.mpf(x); return (gs(x+H_)-gs(x-H_))/(2*H_)/mp.sqrt(fsafe(x))
def Qk(x,k): x=mp.mpf(x); return float(mp.mpf(k*(k+3))/g(x) + gss(x)/g(x))
X0=mp.mpf('1e-6'); X1=1-mp.mpf('1e-6')
def build(N,k):
    t=time.time()
    xs=[X0+(X1-X0)*mp.mpf(i)/N for i in range(N+1)]
    sq=[float(mp.sqrt(fsafe(x))) for x in xs]
    H=[]; 
    for i in range(N):
        a=float(xs[i+1]-xs[i]); H.append(a*0.5*(sq[i]+sq[i+1]))   # trapezoid on sqrt(f); fine at this h
    QL=[Qk(xs[i],k) for i in range(N)]     # left endpoint = exact min on cell (Q increasing)
    h0=float(X0)*float(mp.sqrt(4*pi**2)); q0=min(Qk(X0/10,k),Qk(X0,k))
    h1=float(1-X1)*float(mp.sqrt(FLIM));  q1=Qk(X1,k)
    print(f"   build N={N} took {time.time()-t:.1f}s",flush=True)
    return [h0]+H+[h1], [q0]+QL+[q1]
def shoot(lam,H,QQ,bc0):
    v,vp=(0.0,1.0) if bc0=='D' else (1.0,0.0)
    for h,c in zip(H,QQ):
        k2=c-lam
        if k2>1e-14:
            kap=math.sqrt(k2); a=kap*h
            if a>500.0: a=500.0
            ch=math.cosh(a); sh=math.sinh(a); v,vp = ch*v+sh/kap*vp, kap*sh*v+ch*vp
        elif k2<-1e-14:
            mu=math.sqrt(-k2); ch=math.cos(mu*h); sh=math.sin(mu*h); v,vp = ch*v+sh/mu*vp, -mu*sh*v+ch*vp
        else: v,vp = v+h*vp, vp
        if abs(v)>1e150 or abs(vp)>1e150:
            n=max(abs(v),abs(vp)); v/=n; vp/=n
    return v
def root(H,QQ,bc0,lo,hi):
    fl=shoot(lo,H,QQ,bc0); fh=shoot(hi,H,QQ,bc0)
    assert (fl<0)!=(fh<0), f"no sign change in [{lo},{hi}]"
    for _ in range(80):
        m=0.5*(lo+hi); fm=shoot(m,H,QQ,bc0)
        if (fl<0)!=(fm<0): hi=m
        else: lo=m; fl=fm
    return 0.5*(lo+hi)
UP={'k=0 Dirichlet':0.8203177,'k=1 Neumann gap':0.9988386}
BR={'k=0 Dirichlet':(0,'D',0.815,0.8204),'k=1 Neumann gap':(1,'N',0.993,0.99884)}
N=int(sys.argv[1]) if len(sys.argv)>1 else 60000
print(f"=== certified lower bounds at N={N} ===")
for lbl,(k,bc,lo,hi) in BR.items():
    H,QL=build(N,k)
    r=root(H,QL,bc,lo,hi); u=UP[lbl]
    print(f"{lbl:18s} lower={r:.9f}  upper={u:.7f}  width={u-r:.3e}",flush=True)
