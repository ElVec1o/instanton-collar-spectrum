import numpy as np, mpmath as mp
from scipy.linalg import eigh_tridiagonal
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
mp.mp.dps=35; pi=mp.pi
def fmp(x): x=mp.mpf(x); return 4*pi**2/(x**2-1)**4*( x**4+10*x**2+1 - 12*x**2*(x**2+1)/(x**2-1)*mp.log(x) )
def gmp(x): x=mp.mpf(x); return pi**2/(2*(x**2-1)**2)*( x**4-8*x**2+1 + 24*x**4/(x**4-1)*mp.log(x) )
# 3-pt Gauss-Legendre on [0,1]
gx=np.array([0.5-0.5*np.sqrt(3/5),0.5,0.5+0.5*np.sqrt(3/5)]); gw=np.array([5/18,8/18,5/18])
def pqw_tab(nodes):
    # evaluate p=g^2/sqrt f, q0=sqrt f g (centrifugal per k(k+3)), w=sqrt f g^2 at all gauss pts of all elements
    P=[];Q=[];W=[]
    for e in range(len(nodes)-1):
        xa,xb=nodes[e],nodes[e+1]
        for xi in gx:
            x=xa+(xb-xa)*xi
            F=float(fmp(x)); G=float(gmp(x)); sF=np.sqrt(F)
            P.append(G**2/sF); Q.append(sF*G); W.append(sF*G**2)
    return np.array(P).reshape(-1,3),np.array(Q).reshape(-1,3),np.array(W).reshape(-1,3)
def solve(k,N,bc0,eps=1e-9):
    nodes=np.array([float(mp.mpf(eps)+(mp.mpf(1)-2*mp.mpf(eps))*i/N) for i in range(N+1)])
    Pg,Qg,Wg=pqw_tab(nodes)
    n=N+1
    Kd=np.zeros(n);Ko=np.zeros(n-1); Md=np.zeros(n);Mo=np.zeros(n-1)
    kc=k*(k+3)
    for e in range(N):
        h=nodes[e+1]-nodes[e]
        # shape N1=1-xi,N2=xi ; dNa/dx=-1/h,1/h
        pint=np.dot(gw,Pg[e])*h          # ∫p dx
        # stiffness (1/h^2)*∫p * [[1,-1],[-1,1]]
        s=pint/h**2
        # consistent mass/potential: ∫ (w or q) Na Nb dx  with Gauss
        N1=1-gx; N2=gx
        def mel(vals):
            m11=np.dot(gw, vals*N1*N1)*h; m22=np.dot(gw, vals*N2*N2)*h; m12=np.dot(gw, vals*N1*N2)*h
            return m11,m12,m22
        w11,w12,w22=mel(Wg[e]); q11,q12,q22=mel(Qg[e])
        # A = stiffness + kc*potential ;  B = mass(w)
        Kd[e]+=s+kc*q11; Kd[e+1]+=s+kc*q22; Ko[e]+=-s+kc*q12
        Md[e]+=w11; Md[e+1]+=w22; Mo[e]+=w12
    drop=[]
    if bc0=='D': drop.append(0)
    if k>=1: drop.append(n)  # marker; handle via center regularity: drop last node
    keep=np.array([i for i in range(n) if not (i==0 and bc0=='D') and not (i==n-1 and k>=1)])
    A=sp.diags([Ko,Kd,Ko],[-1,0,1]).tocsr()[keep][:,keep]
    B=sp.diags([Mo,Md,Mo],[-1,0,1]).tocsr()[keep][:,keep]
    vals=eigsh(A,k=6,M=B,sigma=-1e-6,which='LM',return_eigenvectors=False)
    return np.sort(vals.real)
def conv(k,bc0):
    Ns=[3000,6000,12000,24000]
    S=[solve(k,N,bc0) for N in Ns]
    print(f"  k={k} {bc0}: ", end="")
    for j in range(3):
        seq=[S[i][j] for i in range(len(Ns))]
        d=[seq[i+1]-seq[i] for i in range(len(seq)-1)]
        ratio=d[-1]/d[-2] if abs(d[-2])>1e-15 else 0
        rich=seq[-1]+(seq[-1]-seq[-2])/3   # O(h^2)->O(h^4) since consistent-mass is O(h^2)
        print(f"[{rich:.6f} r{ratio:.2f}]", end=" ")
    print()
print("Consistent-mass P1 (upper bounds, genuine O(h^2)); ratio->0.25 confirms O(h^2). Richardson value in []:")
for k in range(0,5):
    conv(k,'D')
for k in range(0,5):
    conv(k,'N')
