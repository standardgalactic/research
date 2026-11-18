"""Exact rational rotating-map enclosure benchmark, with fixed parameter and expenditure."""
from fractions import Fraction as F
import csv
from pathlib import Path
out=Path(__file__).with_name('affine_wrapping_results.csv')
# Rotation matrix R=[[3/5,-4/5],[4/5,3/5]]. Initial (x,y)=(p,0), p in [0.9,1.1].
# Box method forgets shared p at every step; affine method retains (x,y)=v_n*p.
# J_n=sum_{j=0}^{n-1}(x_j^2+y_j^2)=n*p^2 exactly.
a,b=F(3,5),F(4,5);d=F(1,10);p=(1-d,1+d)
def add(u,v):return (u[0]+v[0],u[1]+v[1])
def mul(c,u):return (c*u[0],c*u[1]) if c>=0 else (c*u[1],c*u[0])
rows=[];x=p;y=(F(0),F(0));vx=F(1);vy=F(0)
for n in range(1,25):
 x,y=add(mul(a,x),mul(-b,y)),add(mul(b,x),mul(a,y))
 vx,vy=a*vx-b*vy,b*vx+a*vy
 affine_x=mul(vx,p);affine_y=mul(vy,p)
 box_width=(x[1]-x[0])+(y[1]-y[0]);affine_width=(affine_x[1]-affine_x[0])+(affine_y[1]-affine_y[0])
 rows.append(dict(steps=n,box_total_width=float(box_width),affine_total_width=float(affine_width),ratio=float(box_width/affine_width),J_lower=float(n*(1-d)**2),J_upper=float(n*(1+d)**2),J_true_p1=float(n)))
with out.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
for r in rows:
 if r['steps'] in (4,8,12,16,20,24):print(r)
print('Wrote',out)
