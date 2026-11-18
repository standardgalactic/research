from fractions import Fraction as F
import csv
Q=lambda s:F(str(s))
h=Q('.01'); a=(Q('.15'),Q('.19')); k=Q('.7'); ca=Q('.1'); bits=40
P=(Q('.329'),Q('.331')); M=(Q('.119'),Q('.121')); Z=(F(0),F(0))
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[1],-x[0]
def sub(x,y):return add(x,neg(y))
def mul(x,y):
 p=[v*w for v in x for w in y];return min(p),max(p)
def sc(x,c):return mul(x,(c,c))
def sq(x):return (F(0) if x[0]<=0<=x[1] else min(x[0]**2,x[1]**2),max(x[0]**2,x[1]**2))
def cube(x):return x[0]**3,x[1]**3
def width(x):return x[1]-x[0]
def dyad(x):
 d=1<<bits;return F((x[0]*d).__floor__(),d),F((x[1]*d).__ceil__(),d)
def contains(big,small):return big[0]<=small[0] and big[1]>=small[1]
def inflate(x,r):return x[0]-r,x[1]+r
def physical(Y):p,m,j=Y;return [add(p,m),sub(p,m),j]
def rhs_ph(X):
 u,v,j=X;delta=sub(u,v)
 return [sub(sub(mul(a,u),cube(u)),sc(delta,k)),add(sub(mul(a,v),cube(v)),sc(delta,k)),sc(sq(delta),2*ca*k*k)]
def rhs_modal(Y):
 p,m,j=Y
 return [sub(sub(mul(a,p),cube(p)),sc(mul(p,sq(m)),3)),sub(sub(mul(sub(a,(2*k,2*k)),m),sc(mul(sq(p),m),3)),cube(m)),sc(sq(m),8*ca*k*k)]
def step(X,rhs):
 f=rhs(X);r=[h*max(abs(v[0]),abs(v[1]))+Q('.00001') for v in f]
 for _ in range(32):
  Y=[inflate(x,v) for x,v in zip(X,r)];g=rhs(Y)
  if all(contains(y,add(x,sc(fi,h))) for x,y,fi in zip(X,Y,g)):break
  r=[2*v for v in r]
 else:raise RuntimeError('Picard inclusion failed')
 return [dyad(add(x,sc(fi,h))) for x,fi in zip(X,g)]
# Correlated physical-coordinate affine frame is u=p+m, v=p-m with shared p,m generators.
# This is precisely the same enclosure as modal intervals under the fixed linear map.
X=[dyad(v) for v in physical([P,M,Z])];Y=[dyad(v) for v in (P,M,Z)]
rows=[]
for i in range(81):
 if i in (0,10,20,40,60,80):
  C=physical(Y)
  rows.append(dict(horizon=float(i*h),coordinate_hull_width=float(sum(map(width,X[:2]))),physical_affine_width=float(sum(map(width,C[:2]))),modal_width=float(sum(map(width,C[:2]))),coordinate_hull_J_upper=float(X[2][1]),physical_affine_J_upper=float(C[2][1]),modal_J_upper=float(Y[2][1]),physical_affine_floor_lower=float(Q('.1')-Q('.005')*i*h-C[2][1])))
 if i<80:X=step(X,rhs_ph);Y=step(Y,rhs_modal)
with open('/mnt/data/fair_affine_comparison.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
for r in rows:print(r)
