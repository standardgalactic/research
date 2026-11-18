#!/usr/bin/env python3
"""Exact rational upper certificate for piecewise Riccati disagreement."""
from fractions import Fraction as F
import csv
q=lambda s:F(str(s))
C=8*q('.1')*q('.7')**2
M=q('.121')**2
D0=2*q('.7')-q('.19')+3*q('.325')**2
D1=2*q('.7')-q('.19')+3*q('.329')**2

def exp_neg_bounds(x,n=60):
    assert 0<=x<=4 and n%2==0
    term=F(1); total=term
    for i in range(1,n+1):
        term*= -x/F(i);total+=term
    lower=total-term*x/F(n+1)
    assert 0<lower<=total<1
    return lower,total

def log1p_upper(x,n=31):
    assert 0<=x<1 and n%2==1
    return sum(((-1)**(j+1))*x**j/F(j) for j in range(1,n+1))

def log1p_lower(x,n=32):
    assert 0<=x<1 and n%2==0
    return sum(((-1)**(j+1))*x**j/F(j) for j in range(1,n+1))

elo,ehi=exp_neg_bounds(2*D0)
def v1(e):return D0*M*e/(D0+M*(1-e))
Vlo,Vhi=v1(elo),v1(ehi)
# integral first stage = (1/2) log((D0+M)/(D0+V1));
# second stage = (1/2) log(1+V1/D1).
# In the sum the correction decreases with V1 (D1>D0),
# so Vlo maximizes the total bound.
Jpiece=C/F(2)*(log1p_upper(M/D0) + log1p_upper(Vlo/D1)-log1p_lower(Vlo/D0))
Jric=C/F(2)*log1p_upper(M/D0)
Jconst=C*M/(2*D0)
assert D1>D0 and Jpiece<Jric<Jconst
# This certificate uses p>=.329 after t=1, as established in prior report.
rows=[]
for name,J in [('constant linear dissipation',Jconst),('constant Riccati',Jric),('piecewise Riccati',Jpiece)]:
 rows.append({'method':name,'expenditure_upper':f'{float(J):.15f}', 'energy_lower_all_t':f'{float(q(".1")-J):.15f}', 'improvement_over_constant_riccati_percent':f'{float(100*(Jric-J)/Jric):.8f}'})
with open('/mnt/data/piecewise_riccati_results.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
for r in rows:print(r)
print('D0',D0,'D1',D1,'V1 interval',float(Vlo),float(Vhi),'gain',float(Jric-Jpiece))
