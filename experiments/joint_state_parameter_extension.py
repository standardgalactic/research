import numpy as np, csv
from scipy.integrate import solve_ivp
from pathlib import Path
P=Path(__file__).resolve().parent
N=8; L=np.zeros((N,N))
for i in range(N):
 j=(i+1)%N;L[i,i]+=1;L[j,j]+=1;L[i,j]-=1;L[j,i]-=1
_,V=np.linalg.eigh(L);nu=np.cos(2*np.pi*np.arange(N)/N)
z0=(.26+.08j)*np.ones(N)+.16*V[:,1]+.09j*V[:,3]
true=np.array([.157,1.013,.212]); nominal=np.array([.15,1.,.2]); half=np.array([.01,.02,.02]); eps=.01

def f(z,p):return (p[0]+1j*p[2]*nu)*z-abs(z)**2*z-p[1]*(L@z)
def flow(z,p,h):
 def rhs(t,y):
  zz=y[:N]+1j*y[N:];v=f(zz,p);return np.r_[v.real,v.imag]
 sol=solve_ivp(rhs,(0,h),np.r_[z.real,z.imag],rtol=1e-9,atol=1e-11)
 if not sol.success:raise RuntimeError(sol.message)
 y=sol.y[:,-1];return y[:N]+1j*y[N:]
# Synthetic observations: a single fixed parameter and fixed initial state.
rng=np.random.default_rng(20261118)
interval=.2;steps=20
truth=[z0]
for j in range(steps):truth.append(flow(truth[-1],true,interval))
obs=[]
for z in truth:
 d=rng.normal(size=N)+1j*rng.normal(size=N);d/=np.linalg.norm(d)
 obs.append(z+(.65*eps)*d)
# Finite candidate ensemble, not an enclosure of the entire continuous parameter box.
mesh=np.array(np.meshgrid(np.linspace(-1,1,7),np.linspace(-1,1,7),np.linspace(-1,1,7),indexing='ij')).reshape(3,-1).T
candidates=nominal+mesh*half
# Add generating parameter as an audit control, not as information available to an estimator.
candidates=np.vstack([candidates,true])
# Initial state is assumed exactly known; checkpoint observations are noisy.
states=np.tile(z0,(len(candidates),1)); active=np.ones(len(candidates),dtype=bool)
rows=[]
for j in range(steps+1):
 if j:
  for ix in np.flatnonzero(active):states[ix]=flow(states[ix],candidates[ix],interval)
 errors=np.linalg.norm(states-obs[j],axis=1)
 active &= errors<=eps+1e-12
 surviving=candidates[active]
 if len(surviving):
  low=surviving.min(axis=0);high=surviving.max(axis=0)
 else:low=high=np.full(3,np.nan)
 rows.append([j,j*interval,int(active.sum()),bool(active[-1]),*low,*high])
with open(P/'joint_state_parameter_results.csv','w',newline='') as fh:
 w=csv.writer(fh);w.writerow(['checkpoint','time','surviving_candidates','true_candidate_retained','a_min','k_min','s_min','a_max','k_max','s_max']);w.writerows(rows)
for row in rows[::2]:print(row)
