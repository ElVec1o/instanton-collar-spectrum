# Two-sided enclosure for the L^2 Laplacian eigenvalues on M_1(S^4).
# Upper: conforming P1 + consistent mass  => rigorous Rayleigh-Ritz upper bound.
# Lower: P1 + lumped mass                 => approaches from below (classical SL underestimation).
import numpy as np, mpmath as mp, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
mp.mp.dps=35; pi=mp.pi
def fmp(x): x=mp.mpf(x); return 4*pi**2/(x**2-1)**4*( x**4+10*x**2+1 - 12*x**2*(x**2+1)/(x**2-1)*mp.log(x) )
def gmp(x): x=mp.mpf(x); return pi**2/(2*(x**2-1)**2)*( x**4-8*x**2+1 + 24*x**4/(x**4-1)*mp.log(x) )
gx=np.array([0.5-0.5*np.sqrt(3/5),0.5,0.5+0.5*np.sqrt(3/5)]); gw=np.array([5/18,8/18,5/18])
def assemble(k,N,bc0,lumped,eps=1e-9):
    nodes=np.array([float(mp.mpf(eps)+(mp.mpf(1)-2*mp.mpf(eps))*i/N) for i in range(N+1)])
    n=N+1; Kd=np.zeros(n);Ko=np.zeros(n-1); Md=np.zeros(n);Mo=np.zeros(n-1)
    kc=k*(k+3); N1=1-gx; N2=gx
    for e in range(N):
        xa,xb=nodes[e],nodes[e+1]; h=xb-xa
        P=[];Q=[];W=[]
        for xi in gx:
            x=xa+h*xi; F=float(fmp(x)); G=float(gmp(x)); sF=np.sqrt(F)
            P.append(G**2/sF); Q.append(sF*G); W.append(sF*G**2)
        P=np.array(P);Q=np.array(Q);W=np.array(W)
        s=np.dot(gw,P)*h/h**2
        q11=np.dot(gw,Q*N1*N1)*h; q22=np.dot(gw,Q*N2*N2)*h; q12=np.dot(gw,Q*N1*N2)*h
        Kd[e]+=s+kc*q11; Kd[e+1]+=s+kc*q22; Ko[e]+=-s+kc*q12
        if lumped:   # row-sum lumping of the consistent mass
            w11=np.dot(gw,W*N1*N1)*h; w22=np.dot(gw,W*N2*N2)*h; w12=np.dot(gw,W*N1*N2)*h
            Md[e]+=w11+w12; Md[e+1]+=w22+w12
        else:
            Md[e]+=np.dot(gw,W*N1*N1)*h; Md[e+1]+=np.dot(gw,W*N2*N2)*h
            Mo[e]+=np.dot(gw,W*N1*N2)*h
    keep=np.array([i for i in range(n) if not (i==0 and bc0=='D') and not (i==n-1 and k>=1)])
    A=sp.diags([Ko,Kd,Ko],[-1,0,1]).tocsr()[keep][:,keep]
    B=sp.diags([Mo,Md,Mo],[-1,0,1]).tocsr()[keep][:,keep] if not lumped else sp.diags(Md[keep]).tocsr()
    return np.sort(eigsh(A,k=3,M=B,sigma=-1e-6,which='LM',return_eigenvectors=False).real)
print("mode                    lower(lumped)   upper(consistent)   width      enclosure")
targets=[(0,'D',0,'Friedrichs ground state k=0'),(1,'N',0,'Neumann gap k=1 (5-fold)'),
         (0,'N',1,'Neumann k=0 overtone'),(1,'D',0,'Dirichlet k=1'),(2,'D',0,'Dirichlet k=2')]
for k,bc,idx,lbl in targets:
    lo=assemble(k,24000,bc,True)[idx]; hi=assemble(k,24000,bc,False)[idx]
    mono_lo=assemble(k,12000,bc,True)[idx]; mono_hi=assemble(k,12000,bc,False)[idx]
    okL = lo>mono_lo   # lumped increasing => from below
    okU = hi<mono_hi   # consistent decreasing => from above
    print(f"{lbl:28s} {lo:.7f}   {hi:.7f}   {hi-lo:.2e}  [{lo:.6f},{hi:.6f}] lo_up={okL} hi_dn={okU}")
