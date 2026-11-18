import numpy as np,csv
from scipy.integrate import solve_ivp
from scipy.special import exprel
from pathlib import Path
p=Path(__file__).resolve().parent
N=8;L=np.zeros((N,N))
for i in range(N):
 j=(i+1)%N;L[i,i]+=1;L[j,j]+=1;L[i,j]-=1;L[j,i]-=1
lam,V=np.linalg.eigh(L);a=.15;k=1.;cA=.02
z0=(.26+.08j)*np.ones(N)+.16*V[:,1]+.09j*V[:,3]
B=max(np.linalg.norm(z0),np.sqrt(N*a)); rows=[]
for s in [.02,.05,.10,.20]:
 omega=s*np.cos(2*np.pi*np.arange(N)/N); Om=np.diag(omega)
 K=.5j*(L@Om-Om@L);kap=np.linalg.norm(K,2);eta=kap*B**2
 M=(a+np.max(np.abs(omega))+k*lam[-1])*B+B**3
 Lip=a+np.max(np.abs(omega))+k*lam[-1]+3*B**2
 def field(z):return (a+1j*omega)*z-np.abs(z)**2*z-k*L@z
 def rhs(t,y):
  z=y[:N]+1j*y[N:2*N];dz=field(z)
  return np.r_[dz.real,dz.imag,np.vdot(z,K@z).real,cA*k*k*np.vdot(L@z,L@z).real]
 sol=solve_ivp(rhs,(0,20),np.r_[z0.real,z0.imag,0.,0.],rtol=2e-11,atol=2e-13,dense_output=True)
 assert sol.success
 for H in [.1,.5,2.,20.]:
  def integrate_envelope(dt,method):
   total=0.; t0=0.
   while t0<H-1e-12:
    h=min(dt,H-t0);y=sol.sol(t0);z=y[:N]+1j*y[N:2*N];F0=np.vdot(z,K@z).real
    v=np.linalg.norm(field(z));zz=np.linalg.norm(z)
    nodes,weights=np.polynomial.legendre.leggauss(24)
    ts=(nodes+1)*h/2
    if method=='local':
     rho=np.minimum(2*B,v*ts*np.exp(Lip*ts)*exprel(-Lip*ts)) # v*(exp(Lip*t)-1)/Lip stable
     env=np.minimum(eta,F0+kap*(2*zz*rho+rho*rho))
    elif method=='reset_global':env=np.minimum(eta,F0+2*kap*B*M*ts)
    total+=h/2*np.dot(weights,env)
    t0+=h
   return total
  F0=np.vdot(z0,K@z0).real;slope=2*kap*B*M
  tc=np.clip((eta-F0)/slope,0,H) if slope else H
  initial=F0*tc+.5*slope*tc**2+eta*(H-tc)
  actual=sol.sol(H)[2*N]; actcost=sol.sol(H)[2*N+1]
  for step in [.1,.25,.5]:
   loc=integrate_envelope(step,'local');glob=integrate_envelope(step,'reset_global')
   assert actual<=loc+1e-8,(s,H,step,actual,loc)
   rows.append([s,H,step,actual,eta*H,initial,glob,loc,actcost])
with open(p/'local_frequency_work_results.csv','w',newline='') as f:
 w=csv.writer(f);w.writerow(['spread','horizon','checkpoint_spacing','actual_signed_work','static_upper','initial_state_upper','reset_global_speed_upper','local_state_upper','actual_actuator_cost']);w.writerows(rows)
for r in rows:
 if r[0]==.2 and r[2]==.1:print('H',r[1],'actual',r[3],'static',r[4],'initial',r[5],'reset',r[6],'local',r[7])
