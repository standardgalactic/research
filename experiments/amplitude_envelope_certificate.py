"""Time-dependent amplitude envelope (requires M0=max|z_i(0)|>0; the zero state stays zero).  At a node of maximal modulus the coupling term is nonnegative, so M(t)=max_i|z_i| satisfies
M' <= a M - M^3 (upper Dini derivative), hence M(t)^2 <= R(t)^2 = a/(1 - c e^{-2at}), c = 1 - a/M0^2.
Int_0^T R^2 dt = (1/2) ln((e^{2aT}-c)/(1-c)).  Replacing the constant R by R(t) in the forcing term gives
J(T) <= cA k^2 [ W0 + int(R^2)*sum_j w_j^2 /(4(beta_j-bs)) ]/bs   (mode resolved)  and the Young form with sum w_j^2 aggregated."""
import numpy as np, csv
from scipy.optimize import minimize_scalar
from mode_resolved_certificate import lap,traj,G,a,k,cA
def intR2(M0,T):
    c=1-a/M0**2
    return 0.5*np.log((np.exp(2*a*T)-c)/(1-c))
def run(L,l2,om,z0,T):
    w,v=np.linalg.eigh(L); lam=w[1:].copy(); lam[0]=l2
    omt=om-om.mean(); delta=np.abs(omt).max(); W0=0.5*np.real(np.vdot(z0,L@z0))
    beta=k-a/lam-delta/l2
    if beta.min()<=0: return None
    coef=v[:,1:].T@omt; b2=beta.min(); M0=abs(z0).max()
    R2c=max(M0,np.sqrt(a))**2*T; R2t=intR2(M0,T)
    def young(I):
        al=cA*k*k*W0/b2; kap=cA*k*k*np.sum(omt**2)*I/(4*b2*b2); return (np.sqrt(al+kap)+np.sqrt(kap))**2
    def mode(I):
        f=lambda bs: cA*k*k*(W0+I*np.sum(coef**2/(4*(beta-bs))))/bs
        return minimize_scalar(f,bounds=(1e-9*b2,b2*(1-1e-9)),method='bounded',options={'xatol':1e-14}).fun
    return young(R2c),young(R2t),mode(R2c),mode(R2t),R2t/R2c
def main():
    Ts=[1.0,5.0,20.0,100.0]; rng=np.random.default_rng(11)
    base=np.array([-0.03,-0.01,0.01,0.03])
    rows=[]
    for g,(e,l2) in G.items():
        L=lap(e)
        for scale in (1.0,3.0):
            om=scale*base
            for name,z0 in (('amp_0.45',np.array([.45,.21,.39,.27],complex)),
                            ('amp_0.20',0.2*np.exp(2j*np.pi*rng.uniform(size=4))),
                            ('amp_0.05',0.05*np.exp(2j*np.pi*rng.uniform(size=4))),
                            ('amp_0.01',0.01*np.exp(2j*np.pi*rng.uniform(size=4))),
                            ('amp_0.45_alt',np.array([.45,-.45,.45,-.45],complex))):
                J=traj(L,om,z0,Ts)
                for T,Ja in zip(Ts,J):
                    r=run(L,l2,om,z0,T)
                    if r is None: continue
                    yc,yt,mc,mt,ratio=r
                    rows.append(dict(graph=g,scale=scale,state=name,T=T,J_actual=Ja,young_const_R=yc,young_R_t=yt,mode_const_R=mc,mode_R_t=mt,
                       gain_young_pct=100*(1-yt/yc),gain_mode_pct=100*(1-mt/mc),gain_total_pct=100*(1-mt/yc),int_R2_ratio=ratio,valid=int(min(yc,yt,mc,mt)>=Ja)))
    with open('amplitude_envelope_results.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(len(rows),'all valid',all(r['valid'] for r in rows))
    import collections
    A=collections.defaultdict(list)
    for r in rows: A[(r['state'],r['T'])].append(r['gain_young_pct'])
    for kx in sorted(A): print(kx,'median %.2f max %.2f'%(np.median(A[kx]),max(A[kx])))
    print('max total',max(r['gain_total_pct'] for r in rows))

if __name__=='__main__':
    main()
