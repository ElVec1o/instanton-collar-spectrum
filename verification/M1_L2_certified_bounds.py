import mpmath as mp, math
mp.mp.dps=30; pi=mp.pi
def f(x): x=mp.mpf(x); return 4*pi**2/(x**2-1)**4*( x**4+10*x**2+1 - 12*x**2*(x**2+1)/(x**2-1)*mp.log(x) )
def g(x): x=mp.mpf(x); return pi**2/(2*(x**2-1)**2)*( x**4-8*x**2+1 + 24*x**4/(x**4-1)*mp.log(x) )
H_=mp.mpf('1e-9')
def gs(x):  x=mp.mpf(x); return (g(x+H_)-g(x-H_))/(2*H_)/mp.sqrt(f(x))
def gss(x): x=mp.mpf(x); return (gs(x+H_)-gs(x-H_))/(2*H_)/mp.sqrt(f(x))
def Qk(x,k): x=mp.mpf(x); return float(mp.mpf(k*(k+3))/g(x) + gss(x)/g(x))
X0=mp.mpf('1e-6'); X1=1-mp.mpf('1e-6')
# ---- is Q monotone increasing in x? (then min over a cell = left endpoint, exactly) ----
for k in [0,1,2]:
    grid=[X0+(X1-X0)*mp.mpf(i)/300 for i in range(301)]
    vals=[Qk(x,k) for x in grid]
    inc=all(vals[i+1]>=vals[i]-1e-9 for i in range(len(vals)-1))
    print(f"k={k}: Q monotone increasing on (0,1)? {inc}   Q(0+)={vals[0]:.6f}  Q(1-)={vals[-1]:.3e}")
FLIM=mp.mpf('0.4')*pi**2   # f(x) -> 2 pi^2/5 as x->1 (removable)
def fsafe(x):
    x=mp.mpf(x)
    if x>1-mp.mpf('1e-7'): return FLIM
    if x<mp.mpf('1e-9'):   return 4*pi**2
    return f(x)
def sw(a,b):  # s-width via Simpson on sqrt(f)
    m=(a+b)/2
    return float((b-a)/6*(mp.sqrt(fsafe(a))+4*mp.sqrt(fsafe(m))+mp.sqrt(fsafe(b))))
def build(N,k):
    xs=[X0+(X1-X0)*mp.mpf(i)/N for i in range(N+1)]
    H=[sw(xs[i],xs[i+1]) for i in range(N)]
    QL=[Qk(xs[i],k) for i in range(N)]          # left endpoint = min on cell (Q increasing)
    QM=[Qk((xs[i]+xs[i+1])/2,k) for i in range(N)]
    # end cells covering (0,X0) and (X1,1), using valid lower bounds of Q there
    h0=sw(mp.mpf(0),X0); q0=Qk(X0,k)            # Q increasing => Q>=Q(X0) is FALSE on (0,X0); use Q(0+)
    q0=min(Qk(X0/10,k),Qk(X0,k))
    h1=sw(X1,mp.mpf(1)); q1=Qk(X1,k)   # Q increasing => Q>=Q(X1) on (X1,1). valid lower bd
    return [h0]+H+[h1], [q0]+QL+[q1], [q0]+QM+[q1]
def shoot(lam,H,QQ,bc0):
    v,vp=(0.0,1.0) if bc0=='D' else (1.0,0.0)
    for h,c in zip(H,QQ):
        k2=c-lam
        if k2>1e-14:
            kap=math.sqrt(k2); a=min(kap*h,500.0)
            ch=math.cosh(a); sh=math.sinh(a); v,vp = ch*v+sh/kap*vp, kap*sh*v+ch*vp
        elif k2<-1e-14:
            mu=math.sqrt(-k2); ch=math.cos(mu*h); sh=math.sin(mu*h); v,vp = ch*v+sh/mu*vp, -mu*sh*v+ch*vp
        else: v,vp = v+h*vp, vp
        n=max(abs(v),abs(vp))
        if n>1e150: v/=n; vp/=n
    return v
def eigs(H,QQ,bc0,n=2,lo=-3.0,hi=25.0,step=0.004):
    out=[];lam=lo;prev=shoot(lam,H,QQ,bc0)
    while lam<hi and len(out)<n:
        lam+=step; cur=shoot(lam,H,QQ,bc0)
        if prev!=0 and (prev<0)!=(cur<0):
            a,b=lam-step,lam; fa=prev
            for _ in range(120):
                m=(a+b)/2; fm=shoot(m,H,QQ,bc0)
                if (fa<0)!=(fm<0): b=m
                else: a=m; fa=fm
            out.append((a+b)/2)
        prev=cur
    return out
print(f"\n{'mode':20s} {'CERTIFIED LOWER':>16} {'upper (R-Ritz)':>15} {'certified width':>16}")
UP={'k=0 Dirichlet':0.8203177,'k=1 Neumann gap':0.9988386,'k=1 Dirichlet':2.1806552,'k=0 Neumann n=1':2.0355386}
for k,bc,lbl in [(0,'D','k=0 Dirichlet'),(1,'N','k=1 Neumann gap'),(1,'D','k=1 Dirichlet'),(0,'N','k=0 Neumann n=1')]:
    H,QL,QM=build(6000,k)
    L=eigs(H,QL,bc,2)
    l0 = L[0] if not (bc=='N' and k==0) else L[1]
    u=UP[lbl]; print(f"{lbl:20s} {l0:16.7f} {u:15.7f} {u-l0:16.2e}")
