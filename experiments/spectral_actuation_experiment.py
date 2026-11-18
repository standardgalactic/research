"""Exact spectral vs global contraction bound for linear network consensus."""
import csv
import numpy as np
from pathlib import Path
out = Path(__file__).resolve().parent
n=8
L=np.zeros((n,n))
for i in range(n):
    j=(i+1)%n
    L[i,i]+=1; L[j,j]+=1; L[i,j]-=1; L[j,i]-=1
lam,V=np.linalg.eigh(L)
lam[np.abs(lam)<1e-12]=0
# Initial disagreement only in the first and third nonzero modes.
a=np.zeros(n); a[1]=1.; a[3]=0.35
k=1.; cA=0.01; H=20.
exact=cA*k/2*sum(lam[j]*a[j]**2*(-np.expm1(-2*k*lam[j]*H)) for j in range(1,n))
global_bound=cA*k*k*lam[-1]**2*np.dot(a,a)*(-np.expm1(-2*k*lam[1]*H))/(2*k*lam[1])
assert exact<=global_bound+1e-12
rows=[]
for t in [1,2,5,10,20]:
    exact_t=cA*k/2*sum(lam[j]*a[j]**2*(-np.expm1(-2*k*lam[j]*t)) for j in range(1,n))
    bound_t=cA*k*k*lam[-1]**2*np.dot(a,a)*(-np.expm1(-2*k*lam[1]*t))/(2*k*lam[1])
    rows.append(dict(horizon=t,exact_spectral_cost=exact_t,global_decay_bound=bound_t,overestimate_factor=bound_t/exact_t))
with (out/'spectral_actuation_results.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
print('lambda_2',lam[1],'lambda_N',lam[-1]);print('exact',exact,'global',global_bound,'factor',global_bound/exact)
