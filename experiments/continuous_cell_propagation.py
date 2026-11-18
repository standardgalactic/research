import numpy as np,csv
from scipy.integrate import solve_ivp
from pathlib import Path
P=Path(__file__).resolve().parent
N=3; L=np.array([[1.,-1,0],[-1,2,-1],[0,-1,1]])
nu=np.array([-1.,0.,1.]); k=.7;s=.16
z0=np.array([.45+.10j,.21-.06j,.33+.04j]); a_true=.173
prior=(.12,.22); eps=.003; h=.002; checkpoints=[0,.2,.4,.6,.8]

def f(z,a):return (a+1j*s*nu)*z-abs(z)**2*z-k*(L@z)
def rhs(t,y):
 z=y[:N]+1j*y[N:]; dz=f(z,a_true);return np.r_[dz.real,dz.imag]
sol=solve_ivp(rhs,(0,.8),np.r_[z0.real,z0.imag],t_eval=checkpoints,rtol=1e-12,atol=1e-13)
truth=[v[:N]+1j*v[N:] for v in sol.y.T]
rng=np.random.default_rng(22); obs=[]
for z in truth:
 noise=rng.normal(size=N)+1j*rng.normal(size=N);noise*=.35*eps/np.linalg.norm(noise);obs.append(z+noise)
# Constant-in-time parameter cells. Global Euclidean state/derivative bounds.
a_max=prior[1];B=max(np.linalg.norm(z0),np.sqrt(N*a_max));lam=np.linalg.eigvalsh(L)[-1]
Lip=a_max+3*B*B+k*lam+abs(s);F=(a_max+abs(s)+k*lam+B*B)*B
M=Lip*F
rows=[]; survivors=[]
for ncell in [8,16,32,64,128]:
 width=(prior[1]-prior[0])/ncell; retained=0; min_a=np.inf;max_a=-np.inf; min_r=np.inf; max_r=0
 for j in range(ncell):
  lo=prior[0]+j*width;hi=lo+width; mid=(lo+hi)/2
  z=z0.copy();radius=0.; keep=True
  for ix,t in enumerate(checkpoints):
   if np.linalg.norm(z-obs[ix])>radius+eps:
    keep=False;break
   if ix==len(checkpoints)-1:break
   for step in range(round((checkpoints[ix+1]-t)/h)):
    # Conservative Euclidean ball enclosure of the exact flow about Euler center.
    # Floating point arithmetic is not outward-rounded; extra padding is illustrative.
    radius=(1+h*Lip)*radius+h*(width/2)*B+.5*h*h*M+1e-14
    z=z+h*f(z,mid)
  if keep:
   retained+=1;min_a=min(min_a,lo);max_a=max(max_a,hi);min_r=min(min_r,radius);max_r=max(max_r,radius)
  if lo<=a_true<=hi: truth_retained=keep
 rows.append([ncell,retained,min_a if retained else '',max_a if retained else '',min_r if retained else '',max_r if retained else '',truth_retained])
with open(P/'continuous_cell_results.csv','w',newline='') as fh:
 w=csv.writer(fh);w.writerow(['cells','retained_cells','retained_a_min','retained_a_max','minimum_final_radius','maximum_final_radius','true_parameter_cell_retained']);w.writerows(rows)
for r in rows:print(r)
print('B,Lip,F,M=',B,Lip,F,M)
