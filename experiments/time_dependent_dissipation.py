#!/usr/bin/env python3
"""Exact-rational, outward-bounded staged dissipation certificate."""
from fractions import Fraction as F
import csv
Q=lambda x:F(str(x))
k=Q('.7'); ca=Q('.1'); C=8*ca*k*k
P=Q('.325'); P1=Q('.329'); M=Q('.121')**2
D0=2*k-Q('.19')+3*P*P
D1=2*k-Q('.19')+3*P1*P1
# Alternating Taylor series, even truncation upper and odd lower for exp(-x), x>=0.
def expminus(x,n=40):
 assert x>=0 and n%2==0
 term=F(1);s=term
 for j in range(1,n+1):
  term=term*(-x)/j;s+=term
 lower=s+term*(-x)/(n+1)
 assert lower<=s
 return lower,s
# Starting at t=.5, the lower bound for p' in the strip P <= p <= P1 is
# P*(amin-P1^2-3*M*exp(-2D0*.5)). This is strictly positive.
_,e_half_up=expminus(D0)
v_strip=P*(Q('.15')-P1*P1-3*M*e_half_up)
assert v_strip>0 and v_strip*Q('.5')>P1-P
# Therefore p(t)>=P1 for t>=1, by first crossing and inward vector field.
# Exact piecewise exponential comparison for V=m^2, using D0 before 1 and D1 after.
e1_lo,e1_up=expminus(2*D0)
# J<= C*M/(2D0)*(1-exp(-2D0)) + C*M*exp(-2D0)/(2D1)
# The expression decreases in exp(-2D0) because D1>D0.
Jstage=C*M/F(2)*( (1-e1_lo)/D0 + e1_lo/D1 )
Jconstant=C*M/(2*D0)
# Riccati sharpening: V' <= -2(D0+V)V. Then J_inf <= C/2 log(1+M/D0).
# Alternating log series odd truncation is an upper bound for 0<x<1.
x=M/D0
log_upper=sum(((-1)**(j+1))*x**j/j for j in range(1,16,2) for _ in [0]) if False else sum(((-1)**(j+1))*x**j/F(j) for j in range(1,16))
assert x<1
Jriccati=C*log_upper/2
# Strengthen staged estimate using V' <= -2(D_i+V)V on each segment.
# Let w=1/V; w' >= 2D_i w +2. The maximal V(t) is logistic,
# and the exact integral of V over a segment is (1/2)log((D+Vstart)/(D+Vend)).
# Conservative bound via constant-D Riccati globally suffices as final certificate.
assert Jstage<Jconstant and Jriccati<Jstage
rows=[]
for name,j in [('constant dissipation',Jconstant),('time-dependent two-stage',Jstage),('quadratic Riccati',Jriccati)]:
 rows.append(dict(method=name,expenditure_upper=f'{float(j):.15f}',uniform_energy_lower=f'{float(Q(".1")-j):.15f}',improvement_percent=f'{float(100*(Jconstant-j)/Jconstant):.8f}',floor_pass=Q('.1')-j>Q('.08')))
with open('/mnt/data/time_dependent_dissipation_results.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
for r in rows:print(r)
print('D0,D1',float(D0),float(D1),'strip derivative lower',float(v_strip),'strip displacement lower',float(v_strip/2))
print('exp(-2D0) enclosure',float(e1_lo),float(e1_up))
