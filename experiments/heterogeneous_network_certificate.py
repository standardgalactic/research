from fractions import Fraction as F
from math import sqrt
import csv
from pathlib import Path
P=Path('/mnt/data')
k=F(7,10); a=F(19,100); c=F(1,10); net=F(5,1000)
x=[F(45,100),F(21,100),F(39,100),F(27,100)]
w=[F(-3,100),F(-1,100),F(1,100),F(3,100)]
N=len(x); mean=sum(x)/N; y0sq=sum((v-mean)**2 for v in x)
# certified rational square-root upper bounds by ceiling integer square root
from math import isqrt
def sqrt_upper(q,bits=40):
    n=(q.numerator << (2*bits))//q.denominator
    v=isqrt(n)
    if F(v,1<<bits)**2<q: v+=1
    return F(v,1<<bits)
y0=sqrt_upper(y0sq); R=F(45,100); delta=max(abs(t) for t in w)
B=delta*sqrt_upper(F(N))*R
rows=[]
for name,l2,lmax in [('path_4',F(585,1000),F(3415,1000)),('cycle_4',F(2),F(4)),('complete_4',F(4),F(4))]:
    gamma=k*l2-a; assert gamma>0
    r=B/gamma; C=c*k*k*lmax*lmax
    rate=C*r*r
    if y0>=r:
        s=y0-r
        buffer=C*(2*r*s/gamma+s*s/(2*gamma))
        floor=F(1,10)-buffer if net>=rate else None
    else:
        s=F(0);buffer=F(0);floor=F(1,10) if net>=rate else None
    rows.append(dict(graph=name,gamma=float(gamma),initial_disagreement_upper=float(y0),forcing_upper=float(B),asymptotic_disagreement_upper=float(r),steady_expenditure_rate_upper=float(rate),net_supply=float(net),supply_covers_rate=net>=rate,transient_expenditure_buffer=float(buffer),all_time_energy_floor=(float(floor) if floor is not None else ''),certifies_floor_0_08=(floor>=F(8,100) if floor is not None else False)))
with (P/'heterogeneous_network_results.csv').open('w',newline='') as f:
    dw=csv.DictWriter(f,fieldnames=list(rows[0]));dw.writeheader();dw.writerows(rows)
for row in rows: print(row)
