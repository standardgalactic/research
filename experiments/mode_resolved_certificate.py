"""Mode-resolved direct certificate.  For modes j>=2 with Laplacian eigenvalue lam_j and forcing coefficient w_j=<u_j,omega_tilde>:
beta_j = k - a_+/lam_j - delta/lam2.  For any 0<bs<beta_2:
  J(T) <= cA k^2 [ W0 + T*P(bs) ] / bs,   P(bs)=sum_j R^2 w_j^2 / (4 (beta_j-bs)).
Single-mode forcing recovers the Young form.  Compared with the global Young-optimized direct bound."""
import numpy as np, csv
from scipy.optimize import minimize_scalar
from scipy.integrate import solve_ivp
a,k,cA=0.19,0.7,0.1
G={'path_4':([(0,1),(1,2),(2,3)],0.585),'cycle_4':([(0,1),(1,2),(2,3),(3,0)],2.0),'complete_4':([(i,j) for i in range(4) for j in range(i+1,4)],4.0)}
def lap(e):
    L=np.zeros((4,4))
    for i,j in e: L[i,i]+=1;L[j,j]+=1;L[i,j]-=1;L[j,i]-=1
    return L
def traj(L,om,z0,Ts):
    def f(t,y):
        z=y[:4]+1j*y[4:8];Lz=L@z;dz=(a+1j*om)*z-abs(z)**2*z-k*Lz
        return np.concatenate([dz.real,dz.imag,[cA*k*k*np.sum(abs(Lz)**2)]])
    s=solve_ivp(f,[0,max(Ts)],np.concatenate([z0.real,z0.imag,[0]]),method='DOP853',rtol=1e-11,atol=1e-13,t_eval=Ts)
    return s.y[8]
def certs(L,l2,om,z0,T):
    w,v=np.linalg.eigh(L); lam=w[1:].copy(); lam[0]=min(lam[0],l2) if False else lam[0]
    lam_eff=lam.copy(); lam_eff[0]=l2           # valid lower bound for lambda_2 (as in earlier chapters)
    omt=om-om.mean(); delta=np.abs(omt).max()
    R=max(abs(z0).max(),np.sqrt(a)); W0=0.5*np.real(np.vdot(z0,L@z0))
    beta=k-a/lam_eff-delta/l2
    if beta.min()<=0: return None
    coef=v[:,1:].T@omt
    # group-free: lam_eff used per mode (modes sharing eigenvalue handled automatically)
    F2=R*R*np.sum(omt**2); b2=beta.min()
    al=cA*k*k*W0/b2; kap=cA*k*k*F2*T/(4*b2*b2)
    young=(np.sqrt(al+kap)+np.sqrt(kap))**2
    def obj(bs): return cA*k*k*(W0+T*np.sum(R*R*coef**2/(4*(beta-bs))))/bs
    r=minimize_scalar(obj,bounds=(1e-9*b2,b2*(1-1e-9)),method='bounded',options={'xatol':1e-14})
    return young,r.fun,r.x/b2
def main():
    base=np.array([-0.03,-0.01,0.01,0.03]);alt=np.array([-0.03,0.03,-0.03,0.03])*1.0
    Ts=[1.0,5.0,20.0,100.0]; rng=np.random.default_rng(7)
    rows=[]
    for g,(e,l2) in G.items():
        L=lap(e)
        for pname,om in (('base',base),('base3',3*base),('alternating',alt),('alternating_x3',3*alt)):
            for zn,z0 in (('aligned',np.array([.45,.21,.39,.27],complex)),('alt',np.array([.45,-.45,.45,-.45],complex)),('rand',(0.45*np.sqrt(rng.uniform(size=4))*np.exp(2j*np.pi*rng.uniform(size=4))))):
                J=traj(L,om,z0,Ts)
                for T,Ja in zip(Ts,J):
                    c=certs(L,l2,om,z0,T)
                    if c is None: continue
                    rows.append(dict(graph=g,pattern=pname,state=zn,T=T,J_actual=Ja,young=c[0],mode_resolved=c[1],gain_pct=100*(1-c[1]/c[0]),beta_star_frac=c[2],valid=int(c[1]>=Ja and c[0]>=Ja)))
    with open('mode_resolved_results.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    import collections
    print(len(rows),'rows; all valid:',all(r['valid'] for r in rows))
    agg=collections.defaultdict(list)
    for r in rows: agg[(r['graph'],r['pattern'])].append(r['gain_pct'])
    for kx,v in agg.items(): print(kx,'median gain %.2f max %.2f'%(np.median(v),max(v)))

if __name__=='__main__':
    main()
