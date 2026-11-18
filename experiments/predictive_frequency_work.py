import numpy as np,csv
from scipy.integrate import solve_ivp
from pathlib import Path
p=Path(__file__).resolve().parent
N=8; L=np.zeros((N,N))
for i in range(N):
 j=(i+1)%N; L[i,i]+=1;L[j,j]+=1;L[i,j]-=1;L[j,i]-=1
lam,V=np.linalg.eigh(L); a=.15;k=1.;cA=.02
z0=(.26+.08j)*np.ones(N)+.16*V[:,1]+.09j*V[:,3]
B=max(np.linalg.norm(z0),np.sqrt(N*a)); rows=[]
for s in [0.02,0.05,0.10,0.20]:
 omega=s*np.cos(2*np.pi*np.arange(N)/N); Om=np.diag(omega)
 K=0.5j*(L@Om-Om@L); hnorm=np.linalg.norm(K,2); eta=hnorm*B**2
 # uniform speed bound in rotating frame, omega centered
 M=(a+np.max(np.abs(omega))+k*lam[-1])*B+B**3
 def rhs(t,y):
  z=y[:N]+1j*y[N:2*N];dz=(a+1j*omega)*z-np.abs(z)**2*z-k*L@z
  F=np.vdot(z,K@z).real; Y2=np.vdot(L@z,L@z).real
  return np.r_[dz.real,dz.imag,F,Y2]
 sol=solve_ivp(rhs,(0,20),np.r_[z0.real,z0.imag,0.,0.],rtol=1e-11,atol=1e-13,dense_output=True)
 assert sol.success
 for H in [.1,.5,2.,20.]:
  y=sol.sol(H); J=y[2*N]; actual=cA*k*k*y[2*N+1]
  F0=np.vdot(z0,K@z0).real; slope=2*hnorm*B*M
  # integrate min(eta, F0+slope*t) exactly; always a certified upper envelope
  if slope<1e-14: U=H*min(eta,F0)
  else:
   tclip=np.clip((eta-F0)/slope,0,H)
   U=F0*tclip+.5*slope*tclip*tclip+eta*(H-tclip)
  comm=eta*H
  # exact accumulated signed frequency work should be below U
  assert J <= U+1e-9,(s,H,J,U)
  rows.append([s,H,F0,hnorm,B,M,J,U,comm,U/comm if comm else 0,actual])
with open(p/'predictive_frequency_work_results.csv','w',newline='') as f:
 w=csv.writer(f);w.writerow(['frequency_spread','horizon','initial_signed_frequency_work_rate','commutator_operator_norm','amplitude_bound','uniform_speed_bound','actual_integrated_frequency_work','predictive_initial_state_work_upper','static_commutator_work_upper','predictive_to_static_ratio','actual_actuator_cost']);w.writerows(rows)
for r in rows:print('s %.2f H %4.1f J % .6g U %.6g static %.6g ratio %.3f'%(r[0],r[1],r[6],r[7],r[8],r[9]))
