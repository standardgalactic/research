from fractions import Fraction as F
import csv
Q=lambda s:F(str(s))
h=Q('.01'); a=(Q('.15'),Q('.19')); k=Q('.7'); ca=Q('.1'); bits=40
P=(Q('.329'),Q('.331')); M=(Q('.119'),Q('.121')); Z=(Q(0),Q(0))
def add(x,y):return (x[0]+y[0],x[1]+y[1])
def neg(x):return (-x[1],-x[0])
def sub(x,y):return add(x,neg(y))
def mul(x,y):
 p=[v*w for v in x for w in y];return min(p),max(p)
def sc(x,c):return mul(x,(c,c))
def sq(x):return (F(0) if x[0]<=0<=x[1] else min(x[0]**2,x[1]**2),max(x[0]**2,x[1]**2))
def cube(x):return (x[0]**3,x[1]**3)
def width(x):return x[1]-x[0]
def dyad(x):
 d=1<<bits;return F((x[0]*d).__floor__(),d),F((x[1]*d).__ceil__(),d)
def contains(big,small):return big[0]<=small[0] and big[1]>=small[1]
def inflate(x,r):return (x[0]-r,x[1]+r)
def rhs_ph(X):
 u,v,j=X;delta=sub(u,v)
 return [sub(sub(mul(a,u),cube(u)),sc(delta,k)),add(sub(mul(a,v),cube(v)),sc(delta,k)),sc(sq(delta),2*ca*k*k)]
def rhs_modal(X):
 p,m,j=X
 return [sub(sub(mul(a,p),cube(p)),sc(mul(p,sq(m)),3)),sub(sub(mul(sub(a,(2*k,2*k)),m),sc(mul(sq(p),m),3)),cube(m)),sc(sq(m),8*ca*k*k)]
def step(X,rhs):
 # Picard existence and enclosure over the whole step.
 f=rhs(X);radius=[h*max(abs(v[0]),abs(v[1]))+Q('.00001') for v in f]
 for _ in range(32):
  Y=[inflate(x,r) for x,r in zip(X,radius)];g=rhs(Y)
  if all(contains(y,add(x,sc(fi,h))) for x,y,fi in zip(X,Y,g)):break
  radius=[2*r for r in radius]
 else:raise RuntimeError('Picard inclusion failed')
 return [dyad(add(x,sc(fi,h))) for x,fi in zip(X,g)]
def physical(X):
 p,m,j=X;return [add(p,m),sub(p,m),j]
def report(t,physical_box,modal_box):
 mod_phys=physical(modal_box)
 return dict(horizon=float(t),direct_state_width=float(sum(map(width,physical_box[:2]))),modal_frame_state_width=float(sum(map(width,mod_phys[:2]))),direct_J_upper=float(physical_box[2][1]),modal_J_upper=float(mod_phys[2][1]),direct_floor_lower=float(Q('.1')-Q('.005')*t-physical_box[2][1]),modal_floor_lower=float(Q('.1')-Q('.005')*t-mod_phys[2][1]))
X=[dyad(v) for v in physical([P,M,Z])];Y=[dyad(v) for v in [P,M,Z]]
rows=[report(F(0),X,Y)]
for i in range(1,81):
 X=step(X,rhs_ph);Y=step(Y,rhs_modal)
 if i in (10,20,40,60,80):rows.append(report(i*h,X,Y))
for row in rows:print(row)
with open('/mnt/data/modal_frame_results.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
