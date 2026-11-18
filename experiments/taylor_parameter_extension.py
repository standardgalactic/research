import numpy as np, csv
from scipy.integrate import solve_ivp
N=3; L=np.array([[1,-1,0],[-1,2,-1],[0,-1,1]],float); nu=np.array([-1.,0.,1.]); k=.7;s=.16
z0=np.array([.45+.10j,.21-.06j,.33+.04j]); a0=.17; astar=.173

def rhs(t,y,a):
 z=y[:N]+1j*y[N:2*N]; S=y[2*N:3*N]+1j*y[3*N:4*N]; T=y[4*N:5*N]+1j*y[5*N:6*N]
 f=(a+1j*s*nu)*z-np.abs(z)**2*z-k*L@z
 J=lambda v:(a+1j*s*nu)*v-(2*np.abs(z)**2*v+z*z*np.conj(v))-k*L@v
 fs=J(S)+z
 ft=J(T)- (2*np.conj(z)*S*S+4*z*np.abs(S)**2)+2*S
 return np.r_[f.real,f.imag,fs.real,fs.imag,ft.real,ft.imag]
y0=np.r_[z0.real,z0.imag,np.zeros(4*N)]
sol=solve_ivp(lambda t,y:rhs(t,y,a0),(0,.8),y0,rtol=2e-12,atol=2e-14,dense_output=True)
true=solve_ivp(lambda t,y:rhs(t,y,astar),(0,.8),y0,rtol=2e-12,atol=2e-14,dense_output=True)
B=max(np.linalg.norm(z0),np.sqrt(N*.22)); Lam=.22+3*B**2+s+ k*3
# Uniform second-sensitivity bound (real-arithmetic theorem)
def bound_rhs(t,y):
 S,T=y
 return [Lam*S+B,Lam*T+6*B*S*S+2*S]
bound=solve_ivp(bound_rhs,(0,.8),[0.,0.],rtol=1e-12,atol=1e-14,dense_output=True)
rows=[]
for t in [.2,.4,.6,.8]:
 yy=sol.sol(t); z=yy[:N]+1j*yy[N:2*N]; S=yy[2*N:3*N]+1j*yy[3*N:4*N]
 ztrue=true.sol(t)[:N]+1j*true.sol(t)[N:2*N]
 curvature=bound.sol(t)[1]; obs=.003
 for cells in [8,16,32,64,128]:
  kept=0; kept_true=False
  for i in range(cells):
   lo=.12+.10*i/cells;hi=.12+.10*(i+1)/cells; mid=(lo+hi)/2;d=(hi-lo)/2
   # Taylor expansion about a0, with remainder for full distance to a0
   center=z+S*(mid-a0)
   rad=np.linalg.norm(S)*d+.5*curvature*(abs(mid-a0)+d)**2
   admissible=np.linalg.norm(center-ztrue)<=rad+obs
   kept+=admissible
   if lo<=astar<=hi:kept_true=admissible
  rows.append([t,cells,kept,kept_true,curvature])
with open('/mnt/data/taylor_parameter_results.csv','w',newline='') as f:
 w=csv.writer(f);w.writerow(['time','cells','retained','true_cell_retained','second_sensitivity_bound']);w.writerows(rows)
print('B',B,'Lambda',Lam)
for r in rows:print(r)
