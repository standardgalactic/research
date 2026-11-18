import numpy as np,csv
from scipy.integrate import solve_ivp
from scipy.special import exprel
from numpy.polynomial.legendre import leggauss
from pathlib import Path
p=Path(__file__).resolve().parent
N=8;L=np.zeros((N,N))
for i in range(N):
 j=(i+1)%N;L[i,i]+=1;L[j,j]+=1;L[i,j]-=1;L[j,i]-=1
lam,V=np.linalg.eigh(L);a=.15;k=1.;s=.2
z0=(.26+.08j)*np.ones(N)+.16*V[:,1]+.09j*V[:,3]
omega=s*np.cos(2*np.pi*np.arange(N)/N);Om=np.diag(omega)
K=.5j*(L@Om-Om@L);kap=np.linalg.norm(K,2)
B=max(np.linalg.norm(z0),np.sqrt(N*a));eta=kap*B*B
Lip=abs(a)+np.max(abs(omega))+k*lam[-1]+3*B*B

def field(z):return (a+1j*omega)*z-abs(z)**2*z-k*L@z

def rhs(t,y):
 z=y[:N]+1j*y[N:2*N]; dz=field(z)
 return np.r_[dz.real,dz.imag,np.vdot(z,K@z).real]
sol=solve_ivp(rhs,(0,20),np.r_[z0.real,z0.imag,0.],rtol=2e-11,atol=2e-13,dense_output=True)
assert sol.success
nodes,weights=leggauss(24)
rows=[]
for dt in [.05,.1,.25,.5,1.]:
 for eps in [0.,.001,.005,.01,.025,.05]:
  total=0.; t=0.
  while t<20-1e-12:
   h=min(dt,20-t); y=sol.sol(t); z=y[:N]+1j*y[N:2*N]
   # The true state is at the center of the uncertainty ball, but the certificate only uses the supplied radius.
   # The true trajectory is globally bounded by B, and the center is inside B.
   v=np.linalg.norm(field(z)); F=np.vdot(z,K@z).real
   tau=(nodes+1)*h/2
   rho=np.exp(Lip*tau)*eps+v*tau*np.exp(Lip*tau)*exprel(-Lip*tau)
   env=np.minimum(eta,F+kap*(2*np.linalg.norm(z)*rho+rho*rho))
   total+=h/2*np.dot(weights,env);t+=h
  actual=float(sol.sol(20)[2*N]);rows.append([dt,eps,actual,total,eta*20,total-actual,total/(actual if actual>0 else np.nan)])
  assert total>=actual-1e-8
with open(p/'uncertain_checkpoint_results.csv','w',newline='') as f:
 w=csv.writer(f);w.writerow(['checkpoint_spacing','certified_error_radius','actual_signed_work','local_upper','static_upper','absolute_slack','upper_to_actual_ratio']);w.writerows(rows)
for r in rows:print(' '.join(f'{v:.6f}' for v in r))
