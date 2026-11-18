"""Reproduce the edge-dissipative certificate comparison from the earlier experiment."""
import csv
import numpy as np
from scipy.integrate import solve_ivp
from pathlib import Path
out=Path(__file__).resolve().parent
N=8; L=np.zeros((N,N))
for i in range(N):
 j=(i+1)%N
 L[i,i]+=1; L[j,j]+=1; L[i,j]-=1; L[j,i]-=1
vals,V=np.linalg.eigh(L)
a=.15; k=1.; cA=.02; H=20.; times=np.linspace(0,H,4001)
z0=(.26+.08j)*np.ones(N)+.16*V[:,1]+.09j*V[:,3]
def rhs(t,y):
 z=y[:N]+1j*y[N:]
 dz=(a-np.abs(z)**2)*z-k*L@z
 return np.r_[dz.real,dz.imag]
sol=solve_ivp(rhs,[0,H],np.r_[z0.real,z0.imag],t_eval=times,rtol=2e-10,atol=1e-12)
assert sol.success
z=sol.y[:N].T+1j*sol.y[N:].T
actual=float(np.trapezoid(cA*k*k*np.sum(np.abs(z@L.T)**2,axis=1),times))
Q0=float(np.real(np.vdot(z0,L@z0)))
gamma=k*vals[1]-a
edge_bound=float(cA*k*Q0/2*(1+a/gamma*(-np.expm1(-2*gamma*H))))
rows=[('Observed nonlinear expenditure',actual),('Edge-dissipative certificate',edge_bound)]
with (out/'edge_dissipation_results.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['method','cumulative_cost','factor_over_actual']);w.writerows((name,v,v/actual) for name,v in rows)
assert edge_bound>=actual
print('Q0=',Q0,'gamma=',gamma)
for name,v in rows:print(name,v,'ratio=',v/actual)
