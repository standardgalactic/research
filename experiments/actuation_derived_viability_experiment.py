import numpy as np, csv
from pathlib import Path
out=Path(__file__).resolve().parent
n=8
L=np.zeros((n,n))
for i in range(n):
    j=(i+1)%n
    L[i,i]+=1;L[j,j]+=1;L[i,j]-=1;L[j,i]-=1
lam,V=np.linalg.eigh(L)
l=float(lam[1])
a=1.; H=40.; tau=5.; eps=.1; E0=1.; Estar=.1; q=d=.1; c=2.
ksync=np.log(a/eps)/(l*tau)
def cost(k): return .5*c*a*a*k*l*(-np.expm1(-2*k*l*H))
def energy(k): return E0+(q-d)*H-cost(k)
rows=[]
for k in [0.2,0.5,1.,1.5,2.,3.,4.,ksync]:
    rows.append(dict(k=k,disagreement_at_deadline=a*np.exp(-k*l*tau),
        energy_at_horizon=energy(k), synchronized=a*np.exp(-k*l*tau)<=eps+1e-12,
        viable=energy(k)>=Estar-1e-12))
with (out/'actuation_derived_viability_results.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
print('lambda_2:',l,'k_sync:',ksync,'E(k_sync):',energy(ksync),'feasible:',energy(ksync)>=Estar)
for r in rows: print(r)
