from fractions import Fraction as F
import csv
from pathlib import Path
P=Path('/mnt/data')
k=F(7,10); a=F(19,100); c=F(1,10)
x=[F(45,100),F(21,100),F(39,100),F(27,100)]
mu=sum(x)/len(x); V0=sum((v-mu)**2 for v in x)/2
# Rational conservative eigenvalue bounds, justified in accompanying report.
graphs=[('path_4',F(585,1000),F(3415,1000)),('cycle_4',F(2),F(4)),('complete_4',F(4),F(4))]
rows=[]
for name,l2,lmax in graphs:
    margin=k*l2-a
    assert margin>0
    # two independent inequalities: direct exponential envelope and edge-energy budget
    direct=c*k*k*lmax*lmax*V0/margin
    edge=c*k*lmax*V0*(k*l2)/margin
    bound=min(direct,edge)
    rows.append(dict(graph=name,lambda2_lower=str(l2),lambda_max_upper=str(lmax),V0=str(V0),margin=str(margin),direct_bound=float(direct),edge_bound=float(edge),best_J_bound=float(bound),resource_floor=float(F(1,10)-bound)))
with open(P/'network_dissipation_results.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
print('V0',V0)
for r in rows: print(r)
