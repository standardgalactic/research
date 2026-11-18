from fractions import Fraction as F
from math import isqrt
import csv
from pathlib import Path
import sympy as S
out=Path('/mnt/data');n=4;k=F(7,10);a=F(19,100);delta=F(3,100);R=F(45,100);c=F(1,10);supply=F(5,1000);E0=F(1,10)
x=[F(45,100),F(21,100),F(39,100),F(27,100)]
w=[F(-3,100),F(-1,100),F(1,100),F(3,100)]
def sqrt_up(v,bits=48):
 q=isqrt((v.numerator << (2*bits))//v.denominator)
 if F(q,1<<bits)**2<v:q+=1
 assert F(q,1<<bits)**2>=v
 return F(q,1<<bits)
def lap(edges):
 L=S.zeros(n)
 for i,j in edges:L[i,i]+=1;L[j,j]+=1;L[i,j]-=1;L[j,i]-=1
 return L
rows=[]
for name,edges,l2,lmax in [('path_4',[(0,1),(1,2),(2,3)],F(585,1000),F(3415,1000)),('cycle_4',[(0,1),(1,2),(2,3),(3,0)],F(2),F(4)),('complete_4',[(i,j) for i in range(n) for j in range(i+1,n)],F(4),F(4))]:
 L=lap(edges); xx=S.Matrix([S.Rational(v.numerator,v.denominator) for v in x]); W0s=(xx.T*L*xx)[0]/2; W0=F(int(S.numer(W0s)),int(S.denom(W0s)))
 beta=k-(a+delta)/l2
 assert beta>0
 Fsq=R*R*sum(t*t for t in w)
 # Young's inequality F sqrt(Q) <= beta*Q/2 + F^2/(2beta)
 transient=2*c*k*k*W0/beta
 rate=c*k*k*Fsq/(beta*beta)
 floor=E0-transient if supply>=rate else None
 # Prior graph-aware ceiling: disagreement y^2 <= R^2 h / gamma; Q <= lmax^2 y^2
 Pj=S.eye(n)-S.ones(n)/n;A=k*L-(a+delta)*Pj;inv=(A+S.ones(n)/n).inv()-S.ones(n)/n
 ww=S.Matrix([S.Rational(v.numerator,v.denominator) for v in w]);h=S.simplify((ww.T*inv*ww)[0]);h=F(int(S.numer(h)),int(S.denom(h)))
 gamma=k*l2-a-delta
 oldrate=c*k*k*lmax*lmax*R*R*h/gamma
 rows.append(dict(graph=name,lambda2_lower=float(l2),beta=float(beta),initial_edge_energy=float(W0),direct_transient_expenditure_buffer=float(transient),direct_persistent_power_ceiling=float(rate),previous_graph_aware_power_ceiling=float(oldrate),persistent_power_reduction_pct=float(100*(1-rate/oldrate)),net_supply=float(supply),resource_floor_certified=float(floor) if floor is not None else '',floor_above_0_08=bool(floor is not None and floor>F(8,100))))
with (out/'direct_spectral_expenditure_results.csv').open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
for r in rows:print(r)
