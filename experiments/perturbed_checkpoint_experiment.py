import numpy as np,csv
from scipy.integrate import solve_ivp
from scipy.special import exprel
from numpy.polynomial.legendre import leggauss
from pathlib import Path
P=Path(__file__).resolve().parent
N=8;L=np.zeros((N,N))
for i in range(N):
 j=(i+1)%N;L[i,i]+=1;L[j,j]+=1;L[i,j]-=1;L[j,i]-=1
lam,V=np.linalg.eigh(L); lm=lam[-1]; a0=.15;k0=1.;s=.2
z0=(.26+.08j)*np.ones(N)+.16*V[:,1]+.09j*V[:,3]
w0=s*np.cos(2*np.pi*np.arange(N)/N)
B=max(np.linalg.norm(z0),np.sqrt(N*(a0+.01)))
node,weight=leggauss(32)
def field(z,a,k,w):return (a+1j*w)*z-abs(z)**2*z-k*(L@z)
def kmatrix(w):return .5j*(L@np.diag(w)-np.diag(w)@L)
def run(dt,eps,da,dw,dk):
 # Actual parameters differ from the nominal parameters, with certified componentwise deviations.
 a=a0+da*.7;k=k0+dk*.65
 w=w0+dw*np.sin(2*np.pi*np.arange(N)/N+.37)
 K=kmatrix(w); K0=kmatrix(w0); kap=np.linalg.norm(K0,2)
 mismatch=(da+dw+dk*lm)*B
 Lip=abs(a0)+np.max(abs(w0))+k0*lm+3*(B+eps)**2
 def rhs(t,y):
  z=y[:N]+1j*y[N:2*N];dz=field(z,a,k,w)
  return np.r_[dz.real,dz.imag,np.vdot(z,K@z).real]
 sol=solve_ivp(rhs,(0,20),np.r_[z0.real,z0.imag,0.],rtol=2e-11,atol=2e-13,dense_output=True)
 assert sol.success
 # ||K-K0|| <= ||L|| ||Delta Omega|| <= lm*dw
 eta=(kap+lm*dw)*B**2
 total=0.;max_observation_error=0.;t=0.;j=0
 while t<20-1e-12:
  h=min(dt,20-t);y=sol.sol(t);z=y[:N]+1j*y[N:2*N]
  direction=np.exp(1j*(.71*np.arange(N)+1.37*j))
  direction/=np.linalg.norm(direction)
  zhat=z+eps*.85*direction
  max_observation_error=max(max_observation_error,np.linalg.norm(zhat-z))
  v=np.linalg.norm(field(zhat,a0,k0,w0))
  F0=np.vdot(zhat,K0@zhat).real
  tau=(node+1)*h/2
  rho=np.exp(Lip*tau)*eps+(v+mismatch)*tau*np.exp(Lip*tau)*exprel(-Lip*tau)
  # True state norm bounded by B, and true frequency matrix differs by at most dw.
  env=np.minimum(eta,F0+kap*(2*np.linalg.norm(zhat)*rho+rho*rho)+lm*dw*B**2)
  total+=h/2*np.dot(weight,env)
  t+=h;j+=1
 actual=float(sol.sol(20)[2*N]);assert total>=actual-1e-8,(dt,eps,da,dw,dk,total,actual)
 return [dt,eps,da,dw,dk,actual,total,eta*20,total/actual,max_observation_error]
rows=[]
for dt in [.05,.1,.25,.5,1.]:
 for eps in [0,.005,.01,.025]:
  for scale in [0,.5,1]:rows.append(run(dt,eps,.01*scale,.02*scale,.02*scale))
with open(P/'perturbed_checkpoint_results.csv','w',newline='') as f:
 w=csv.writer(f);w.writerow(['checkpoint_spacing','observation_radius','growth_uncertainty','frequency_uncertainty','coupling_uncertainty','actual_signed_work','certified_upper','static_upper','upper_to_actual','maximum_observed_error']);w.writerows(rows)
for row in rows:
 if row[0] in [.1,.5,1.] and row[1] in [.005,.01] and row[3] in [0,.02]:print(' '.join(f'{v:.6f}' for v in row))
