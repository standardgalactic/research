from fractions import Fraction as F
from math import isqrt
import csv
from pathlib import Path
import sympy as S
P=Path('/mnt/data'); n=4;k=F(7,10);a=F(19,100);delta=F(3,100);R=F(45,100); c=F(1,10)
w=[F(-3,100),F(-1,100),F(1,100),F(3,100)]
def sqhi(v):
 t=(v.numerator<<80)//v.denominator; q=isqrt(t);q+=F(q,1<<40)**2<v;return F(q,1<<40)
def lap(edges):
 L=S.zeros(n)
 for i,j in edges:L[i,i]+=1;L[j,j]+=1;L[i,j]-=1;L[j,i]-=1
 return L
rows=[]
for name,edges,l2,lmax in [('path_4',[(0,1),(1,2),(2,3)],F(585,1000),F(3415,1000)),('cycle_4',[(0,1),(1,2),(2,3),(3,0)],F(2),F(4)),('complete_4',[(i,j) for i in range(n) for j in range(i+1,n)],F(4),F(4))]:
 L=lap(edges);Pj=S.eye(n)-S.ones(n)/n
 gamma=k*l2-a-delta
 assert gamma>0
 # A=kL -(a+delta)P, invert restricted to mean-zero using projector shift
 A=k*L-(a+delta)*Pj
 inv=(A+S.ones(n)/n).inv()-S.ones(n)/n
 omega=S.Matrix([S.Rational(x.numerator,x.denominator) for x in w]);omega=Pj*omega
 h=S.simplify((omega.T*inv*omega)[0]); h=F(int(S.numer(h)),int(S.denom(h)))
 bound=R*sqhi(h/gamma)
 old=delta*sqhi(F(n))*R/(k*l2-a)
 # refined norm forcing without graph inverse, using delta y term
 simple=R*sqhi(sum(x*x for x in w))/(k*l2-a-delta)
 rows.append(dict(graph=name,lambda2_lower=float(l2),lambda_max_upper=float(lmax),graph_quadratic_form=float(h),graph_aware_disagreement_ceiling=float(bound),original_disagreement_ceiling=float(old),simple_centered_disagreement_ceiling=float(simple),graph_vs_original_pct=float(100*(1-bound/old)),graph_vs_simple_pct=float(100*(1-bound/simple))))
with (P/'graph_aware_frequency_results.csv').open('w',newline='') as f:
 d=csv.DictWriter(f,fieldnames=rows[0].keys());d.writeheader();d.writerows(rows)
for row in rows:print(row)
