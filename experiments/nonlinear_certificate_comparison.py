import numpy as np, csv
from scipy.integrate import solve_ivp
from pathlib import Path
N=8; L=np.zeros((N,N))
for i in range(N):
 j=(i+1)%N; L[i,i]+=1;L[j,j]+=1;L[i,j]-=1;L[j,i]-=1
vals,V=np.linalg.eigh(L)
a=.15;k=1.;cA=.02;H=20.;times=np.linspace(0,H,4001)
# Real and imaginary modal content, with nonzero collective mean
z0=(.26+.08j)*np.ones(N)+.16*V[:,1]+.09j*V[:,3]
max0=np.max(np.abs(z0)); maxamp=max(max0,np.sqrt(a)); B=np.sqrt(N)*maxamp
w0=np.linalg.norm(z0-z0.mean()); gamma=k*vals[1]-a; C=a+3*B**2
assert gamma>0

def rhs(t,y):
 z=y[:N]+1j*y[N:]; dz=(a-np.abs(z)**2)*z-k*L@z
 return np.r_[dz.real,dz.imag]
sol=solve_ivp(rhs,[0,H],np.r_[z0.real,z0.imag],t_eval=times,rtol=2e-10,atol=1e-12)
assert sol.success
z=sol.y[:N].T+1j*sol.y[N:].T
w=z-z.mean(axis=1,keepdims=True)
act=np.sum(np.abs(z@L.T)**2,axis=1)*cA*k*k
actual=np.trapezoid(act,times)
W=w0*np.exp(-gamma*times)
global_bound=cA*k*k*vals[-1]**2*np.trapezoid(W**2,times)
# Mode envelopes: constant forcing bound M_j=C*B; time-decay bound C*W(t).
# Because common frequency and delta=0, only the nonlinear and growth terms remain.
mode0=np.abs(V.T@z0)
constant=0.;decay=0.
for j in range(1,N):
 r=k*vals[j]
 Uconst=mode0[j]*np.exp(-r*times)+(C*B)/r*(1-np.exp(-r*times))
 if abs(r-gamma)<1e-9: convolution=times*np.exp(-r*times)
 else: convolution=(np.exp(-gamma*times)-np.exp(-r*times))/(r-gamma)
 Udecay=mode0[j]*np.exp(-r*times)+C*w0*convolution
 constant+=cA*k*k*vals[j]**2*np.trapezoid(Uconst**2,times)
 decay+=cA*k*k*vals[j]**2*np.trapezoid(Udecay**2,times)
rows=[('Actual nonlinear expenditure',actual),('Global transverse bound',global_bound),('Constant-forcing spectral bound',constant),('Decay-sensitive spectral bound',decay)]
assert all(v+1e-8>=actual for _,v in rows[1:]),rows
with open('/mnt/data/nonlinear_certificate_comparison.csv','w',newline='') as f:
 wr=csv.writer(f);wr.writerow(['method','cumulative_cost','factor_over_actual']);wr.writerows((name,float(v),float(v/actual)) for name,v in rows)
print('gamma',gamma,'B',B,'max observed',np.max(np.abs(z)),'w0',w0)
for row in rows:print(row[0],row[1],row[1]/actual)
