import numpy as np,csv
from scipy.integrate import solve_ivp,quad
from pathlib import Path
p=Path(__file__).resolve().parent
N=8; L=np.zeros((N,N))
for i in range(N):
 j=(i+1)%N; L[i,i]+=1;L[j,j]+=1;L[i,j]-=1;L[j,i]-=1
lam,V=np.linalg.eigh(L); a=.15;k=1.;cA=.02;H=20.
z0=(.26+.08j)*np.ones(N)+.16*V[:,1]+.09j*V[:,3]
Q0=float(np.real(np.vdot(z0,L@z0))); B=max(np.linalg.norm(z0),np.sqrt(N*a)); gamma=k*lam[1]-a
rows=[]
for spread in [0,.02,.05,.1,.2]:
 omega=spread*np.cos(2*np.pi*np.arange(N)/N); delta=max(abs(omega-omega.mean()))
 def rhs(t,y):
  z=y[:N]+1j*y[N:2*N]; dz=(a+1j*omega-np.abs(z)**2)*z-k*L@z
  return np.r_[dz.real,dz.imag,cA*k*k*np.linalg.norm(L@z)**2]
 sol=solve_ivp(rhs,[0,H],np.r_[z0.real,z0.imag,0.],rtol=1e-10,atol=1e-12)
 assert sol.success
 observed=sol.y[-1,-1]
 h0=np.sqrt(Q0); r=delta*B*np.sqrt(lam[-1])/gamma
 # h(t) <= r+(h0-r) exp(-gamma t)
 I=r*r*H+2*r*(h0-r)*(-np.expm1(-gamma*H))/gamma+(h0-r)**2*(-np.expm1(-2*gamma*H))/(2*gamma)
 M=Q0/2+a*I
 D=delta*B*np.sqrt(H)
 bound=cA/4*(D+np.sqrt(D*D+4*k*M))**2
 assert bound+1e-9>=observed,(spread,bound,observed)
 rows.append((spread,delta,observed,bound,bound/observed,Q0,B,gamma))
with (p/'heterogeneous_edge_results.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['frequency_spread','delta','observed_cost','edge_certificate','overestimate_factor','initial_edge_energy','amplitude_bound','contraction_margin']);w.writerows(rows)
for row in rows:print('spread %.3f observed %.9f bound %.9f ratio %.2f'%(row[0],row[2],row[3],row[4]))
