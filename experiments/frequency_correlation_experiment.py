import csv
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
p=Path(__file__).resolve().parent
N=8; L=np.zeros((N,N))
for i in range(N):
 j=(i+1)%N;L[i,i]+=1;L[j,j]+=1;L[i,j]-=1;L[j,i]-=1
lam,V=np.linalg.eigh(L);a=.15;k=1.;cA=.02;H=20.
z0=(.26+.08j)*np.ones(N)+.16*V[:,1]+.09j*V[:,3]
Q0=np.vdot(z0,L@z0).real;B=max(np.linalg.norm(z0),np.sqrt(N*a));gamma=k*lam[1]-a
rows=[]
for s in [0,.02,.05,.1,.2]:
 omega=s*np.cos(2*np.pi*np.arange(N)/N);Om=np.diag(omega);delta=max(abs(omega-omega.mean()))
 # State augmented by integrated quantities: actuation Y², Q, F, positive F, cubic edge dissipation
 def rhs(t,y):
  z=y[:N]+1j*y[N:2*N];Lz=L@z;g=np.abs(z)**2*z
  dz=(a+1j*omega)*z-g-k*Lz
  Y2=np.vdot(Lz,Lz).real;Q=np.vdot(z,Lz).real
  F=np.vdot(Lz,1j*omega*z).real
  T=np.vdot(Lz,g).real
  return np.r_[dz.real,dz.imag,Y2,Q,F,max(F,0.),T]
 sol=solve_ivp(rhs,(0,H),np.r_[z0.real,z0.imag,np.zeros(5)],rtol=1e-11,atol=1e-13,dense_output=True)
 assert sol.success
 X2,IQ,JF,JFpos,IT=sol.y[2*N:,-1]
 zH=sol.y[:N,-1]+1j*sol.y[N:2*N,-1];QH=np.vdot(zH,L@zH).real
 cost=cA*k*k*X2
 h0=np.sqrt(Q0);r=delta*B*np.sqrt(lam[-1])/gamma
 I=r*r*H+2*r*(h0-r)*(-np.expm1(-gamma*H))/gamma+(h0-r)**2*(-np.expm1(-2*gamma*H))/(2*gamma)
 M=Q0/2+a*I;D=delta*B*np.sqrt(H)
 old=cA/4*(D+np.sqrt(D*D+4*k*M))**2
 # conditional a-posteriori signed-correlation bound with actual measured J_F and upper Q envelope I
 conditional=cA*k*max(0.,M+JF)
 # conditional positive forcing budget only (also trajectory dependent)
 positive=cA*k*(M+JFpos)
 # a priori commutator estimate F<= ||i[L,Omega]|| ||z||²/2
 comm=L@Om-Om@L;eta=.5*np.linalg.norm(comm,2)*B*B
 comm_bound=cA*k*(M+eta*H)
 # identity: k integral Y²=Q0/2-QH/2+a IQ - IT + JF
 identity=Q0/2-QH/2+a*IQ-IT+JF
 assert np.isclose(k*X2,identity,atol=1e-8,rtol=1e-8),(s,X2,identity)
 assert min(old,conditional,positive,comm_bound)+1e-9>=cost
 rows.append([s,cost,old,conditional,positive,comm_bound,JF,JFpos,IT,IQ,QH,eta,old/cost,conditional/cost])
with (p/'frequency_correlation_results.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['frequency_spread','observed_cost','prior_heterogeneous_certificate','signed_correlation_conditional_certificate','positive_forcing_conditional_certificate','commutator_apriori_certificate','integrated_signed_frequency_work','integrated_positive_frequency_work','integrated_cubic_dissipation','integrated_edge_energy','final_edge_energy','commutator_source_bound','prior_overestimate','conditional_overestimate']);w.writerows(rows)
for x in rows:print('s=%.2f cost=%.8f prior=%.8f signed=%.8f positive=%.8f comm=%.8f J=%.7f ratio=%.2f'%(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[13]))
