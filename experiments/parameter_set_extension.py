import numpy as np,csv
from scipy.integrate import solve_ivp
from scipy.optimize import linprog
from pathlib import Path
P=Path(__file__).resolve().parent
N=8; L=np.zeros((N,N))
for i in range(N):
 j=(i+1)%N;L[i,i]+=1;L[j,j]+=1;L[i,j]-=1;L[j,i]-=1
lam,V=np.linalg.eigh(L); lm=lam[-1]; template=np.cos(2*np.pi*np.arange(N)/N)
a0,k0,s0=.15,1.,.2; da,dk,ds=.01,.02,.02
actual=np.array([a0+.007,k0+.013,s0+.012]); z0=(.26+.08j)*np.ones(N)+.16*V[:,1]+.09j*V[:,3]
B=max(np.linalg.norm(z0),np.sqrt(N*(a0+da)))
def f(z,p):return (p[0]+1j*p[2]*template)*z-abs(z)**2*z-p[1]*(L@z)
def rhs(t,y):
 z=y[:N]+1j*y[N:];v=f(z,actual);return np.r_[v.real,v.imag]
sol=solve_ivp(rhs,(0,20),np.r_[z0.real,z0.imag],rtol=2e-11,atol=2e-13,dense_output=True)
assert sol.success
bounds=[(a0-da,a0+da),(k0-dk,k0+dk),(s0-ds,s0+ds)]
# Global bounds valid for all admissible parameters and all true states in the B ball.
Lip=(a0+da)+(s0+ds)+ (k0+dk)*lm+3*B**2
M=((a0+da)+(s0+ds)+(k0+dk)*lm)*B+B**3
# Error in a forward finite-difference velocity using noisy endpoints, and model evaluated at noisy center.
def observe(t,eps,seed):
 z=sol.sol(t);z=z[:N]+1j*z[N:]
 direction=np.exp(1j*(.53*np.arange(N)+.91*seed));direction/=np.linalg.norm(direction)
 return z+.7*eps*direction

def fit(eps,h,spacing):
 A=[]; b=[]; ts=np.arange(0, min(4.,20-h)+1e-9,spacing)
 for j,t in enumerate(ts):
  x=observe(t,eps,2*j); y=observe(t+h,eps,2*j+1)
  v=(y-x)/h+abs(x)**2*x
  features=np.stack([x,-L@x,1j*template*x],axis=1)
  # Norm bound implies each real coordinate is bounded by R.
  R=2*eps/h+Lip*eps+.5*Lip*M*h
  for comp in [np.real,np.imag]:
   A.extend(comp(features)); b.extend(comp(v)+R)
   A.extend(-comp(features));b.extend(-comp(v)+R)
 A=np.array(A);b=np.array(b)
 out=[]
 for i in range(3):
  c=np.eye(3)[i]; lo=linprog(c,A_ub=A,b_ub=b,bounds=bounds,method='highs');hi=linprog(-c,A_ub=A,b_ub=b,bounds=bounds,method='highs')
  if not (lo.success and hi.success):raise RuntimeError((eps,h,spacing,lo.message,hi.message))
  out.append((lo.fun,-hi.fun))
 assert all(l-1e-8<=p<=u+1e-8 for (l,u),p in zip(out,actual))
 return out,len(ts),R
rows=[]
for eps in [0,1e-6,1e-5,1e-4]:
 for h in [.001,.005,.02]:
  out,n,R=fit(eps,h,.1)
  rows.append([eps,h,n,R,*[v for interval in out for v in interval],(out[2][1]-out[2][0])/(2*ds)])
with open(P/'parameter_set_results.csv','w',newline='') as fh:
 w=csv.writer(fh);w.writerow(['observation_radius','difference_step','number_pairs','certified_residual_radius','a_min','a_max','k_min','k_max','s_min','s_max','frequency_interval_width_fraction']);w.writerows(rows)
for r in rows:print(' '.join(f'{v:.7g}' for v in r))
# Independent velocity measurements: explicitly different information from state-only checkpoints.
with open(P/'parameter_set_velocity_results.csv','w',newline='') as fh:
 w=csv.writer(fh);w.writerow(['velocity_error_radius','a_min','a_max','k_min','k_max','s_min','s_max','frequency_width_fraction'])
 for ev in [1e-5,1e-4,1e-3,.01]:
  A=[];b=[]
  for j,t in enumerate(np.arange(0,4.001,.1)):
   v=sol.sol(t);z=v[:N]+1j*v[N:]
   # Velocity data are independent observations with a certified error ball.
   direction=np.exp(1j*(.37*np.arange(N)+.7*j));direction/=np.linalg.norm(direction)
   measured=f(z,actual)+.7*ev*direction
   target=measured+abs(z)**2*z
   features=np.stack([z,-L@z,1j*template*z],axis=1)
   for comp in [np.real,np.imag]:
    A.extend(comp(features));b.extend(comp(target)+ev)
    A.extend(-comp(features));b.extend(-comp(target)+ev)
  A=np.array(A);b=np.array(b);intervals=[]
  for i in range(3):
   c=np.eye(3)[i];lo=linprog(c,A_ub=A,b_ub=b,bounds=bounds,method='highs');hi=linprog(-c,A_ub=A,b_ub=b,bounds=bounds,method='highs')
   assert lo.success and hi.success
   intervals.append((lo.fun,-hi.fun))
  assert all(l-1e-8<=p<=u+1e-8 for (l,u),p in zip(intervals,actual))
  w.writerow([ev,*[v for interval in intervals for v in interval],(intervals[2][1]-intervals[2][0])/(2*ds)])
  print('velocity',ev,intervals)
