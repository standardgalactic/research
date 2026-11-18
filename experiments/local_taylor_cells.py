import numpy as np,csv
from scipy.integrate import solve_ivp
N=3;L=np.array([[1,-1,0],[-1,2,-1],[0,-1,1.]],float);nu=np.array([-1.,0.,1.]);k=.7;s=.16
z0=np.array([.45+.10j,.21-.06j,.33+.04j]); astar=.173; lo0=.12;hi0=.22; obs=.003
B=max(np.linalg.norm(z0),np.sqrt(N*hi0));Lam=hi0+3*B**2+s+k*3

def rhs(t,y,a):
 z=y[:N]+1j*y[N:2*N];S=y[2*N:3*N]+1j*y[3*N:4*N]
 f=(a+1j*s*nu)*z-np.abs(z)**2*z-k*L@z
 JS=(a+1j*s*nu)*S-(2*np.abs(z)**2*S+z*z*np.conj(S))-k*L@S
 return np.r_[f.real,f.imag,(JS+z).real,(JS+z).imag]
y0=np.r_[z0.real,z0.imag,np.zeros(2*N)]
times=[.2,.4,.6,.8];true=solve_ivp(lambda t,y:rhs(t,y,astar),(0,.8),y0,rtol=1e-12,atol=1e-14,dense_output=True)
def major(t,y):
 U,V=y;return [Lam*U+B,Lam*V+6*B*U*U+2*U]
maj=solve_ivp(major,(0,.8),[0.,0.],rtol=1e-12,atol=1e-14,dense_output=True)
rows=[]
for n in [8,16,32,64,128]:
 persistent=np.ones(n,dtype=bool)
 for t in times:
  yt=true.sol(t);zt=yt[:N]+1j*yt[N:2*N];V=maj.sol(t)[1]
  retained=np.zeros(n,dtype=bool)
  for i in range(n):
   low=lo0+(hi0-lo0)*i/n;high=lo0+(hi0-lo0)*(i+1)/n;m=(low+high)/2;d=(high-low)/2
   sol=solve_ivp(lambda tt,y:rhs(tt,y,m),(0,t),y0,rtol=2e-11,atol=2e-13)
   y=sol.y[:,-1];z=y[:N]+1j*y[N:2*N];S=y[2*N:3*N]+1j*y[3*N:4*N]
   rad=np.linalg.norm(S)*d+.5*V*d*d
   retained[i]=np.linalg.norm(z-zt)<=rad+obs
  persistent &= retained
  true_i=min(n-1,int((astar-lo0)/(hi0-lo0)*n))
  rows.append([n,t,int(retained.sum()),int(persistent.sum()),bool(retained[true_i]),bool(persistent[true_i]),V])
  print(rows[-1],flush=True)
with open('/mnt/data/local_taylor_cells_results.csv','w',newline='') as f:
 w=csv.writer(f);w.writerow(['cells','time','checkpoint_retained','persistent_retained','true_cell_checkpoint_retained','true_cell_persistent_retained','global_second_sensitivity_bound']);w.writerows(rows)
