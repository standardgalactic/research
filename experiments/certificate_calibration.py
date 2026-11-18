"""Calibration of network-level expenditure certificates against integrated nonlinear trajectories.
Model: zdot=(a+i*Omega)z-|z|^2 z-kLz, Jdot=cA k^2 ||Lz||^2.  Four-node path, cycle, complete graphs."""
import numpy as np, csv
from scipy.integrate import solve_ivp
a,k,cA=0.19,0.7,0.1
graphs={'path_4':([(0,1),(1,2),(2,3)],0.585,3.415),
        'cycle_4':([(0,1),(1,2),(2,3),(3,0)],2.0,4.0),
        'complete_4':([(i,j) for i in range(4) for j in range(i+1,4)],4.0,4.0)}
def lap(edges):
    L=np.zeros((4,4))
    for i,j in edges: L[i,i]+=1;L[j,j]+=1;L[i,j]-=1;L[j,i]-=1
    return L
def integrate(L,om,z0,Ts):
    def f(t,y):
        z=y[:4]+1j*y[4:8]
        dz=(a+1j*om)*z-np.abs(z)**2*z-k*(L@z)
        Lz=L@z
        return np.concatenate([dz.real,dz.imag,[cA*k*k*np.sum(np.abs(Lz)**2)]])
    y0=np.concatenate([z0.real,z0.imag,[0.0]])
    sol=solve_ivp(f,[0,max(Ts)],y0,method='DOP853',rtol=1e-11,atol=1e-13,t_eval=Ts)
    return sol.y[8]
def bounds(L,l2,lmax,om,z0,T):
    N=4; omt=om-om.mean(); delta=np.max(np.abs(omt))
    R=max(np.max(np.abs(z0)),np.sqrt(max(a,0)))
    F=R*np.linalg.norm(omt); beta=k-(a+delta)/l2
    W0=0.5*np.real(np.vdot(z0,L@z0))
    if beta<=0: return None
    al=cA*k*k*W0/beta; kap=cA*k*k*F*F*T/(4*beta*beta)
    half=2*al+4*kap
    opt=(np.sqrt(al+kap)+np.sqrt(kap))**2
    # earlier network bound (disagreement-norm route)
    gam=k*l2-a; Bf=delta*np.sqrt(N)*R; r=Bf/gam
    y0=np.linalg.norm(z0-z0.mean()); C=cA*k*k*lmax**2
    if gam<=0: old=None
    elif y0>=r:
        s=y0-r; old=C*(r*r*T+2*r*s/gam*(1-np.exp(-gam*T))+s*s/(2*gam)*(1-np.exp(-2*gam*T)))
    else: old=C*r*r*T
    return half,opt,old
rng=np.random.default_rng(20261118)
base=np.array([-0.03,-0.01,0.01,0.03])
def cases(L):
    w,v=np.linalg.eigh(L)
    cs={}
    cs['reference_real']=np.array([0.45,0.21,0.39,0.27],dtype=complex)
    cs['reference_complex']=np.array([0.45*np.exp(0.3j),0.21*np.exp(-1.1j),0.39*np.exp(2.0j),0.27*np.exp(0.7j)])
    for i in range(3):
        r=0.45*np.sqrt(rng.uniform(size=4)); th=rng.uniform(0,2*np.pi,4); cs[f'random_{i+1}']=r*np.exp(1j*th)
    cs['alternating_max']=np.array([0.45,-0.45,0.45,-0.45],dtype=complex)
    cs['fiedler_mode']=(0.45*v[:,1]/np.max(np.abs(v[:,1]))).astype(complex)
    cs['near_synchronized']=np.array([0.45,0.4505,0.4495,0.45],dtype=complex)
    cs['small_amplitude']=0.01*np.exp(1j*rng.uniform(0,2*np.pi,4))
    return cs
Ts=[1.0,5.0,20.0,100.0]
rows=[]
for gname,(edges,l2,lmax) in graphs.items():
    L=lap(edges)
    for scale in (1.0,3.0):
        om=scale*base
        for cname,z0 in cases(L).items():
            J=integrate(L,om,z0,Ts)
            for T,Ja in zip(Ts,J):
                b=bounds(L,l2,lmax,om,z0,T)
                if b is None: continue
                half,opt,old=b
                assert Ja<=opt*(1+1e-9) and Ja<=half*(1+1e-9) and (old is None or Ja<=old*(1+1e-9)),(gname,cname,T,Ja,opt,half,old)
                rows.append([gname,scale,cname,T,Ja,half,opt,old if old is not None else '',half/Ja if Ja>0 else '',opt/Ja if Ja>0 else '',(old/Ja) if (old is not None and Ja>0) else ''])
with open('certificate_calibration_results.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['graph','frequency_scale','initial_state','T','J_actual','bound_theta_half','bound_theta_optimal','bound_disagreement_route','ratio_half','ratio_optimal','ratio_disagreement_route']); w.writerows(rows)
print(len(rows),'rows; all bounds verified >= actual')
import collections
agg=collections.defaultdict(list)
for r in rows:
    if r[2] in('reference_real',): print(r[:7],r[8:10])
for gname in graphs:
  for T in Ts:
    sel=[r for r in rows if r[0]==gname and r[3]==T and r[1]==1.0 and r[4]>1e-6]
    if sel: print(gname,T,'median ratio half %.1f opt %.1f max opt %.1f'%(np.median([r[8] for r in sel]),np.median([r[9] for r in sel]),max(r[9] for r in sel)))
