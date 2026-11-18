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
def modal(Y):
 p,m,j=Y
 return [sub(sub(mul(a,p),cube(p)),sc(mul(p,sq(m)),3)),sub(sub(mul(sub(a,(2*k,2*k)),m),sc(mul(sq(p),m),3)),cube(m)),sc(sq(m),8*ca*k*k)]
def rhs_shear(Y,s):
 v,m,j=Y;p=add(v,sc(m,s)); fp,fm,fj=modal([p,m,j]);return [sub(fp,sc(fm,s)),fm,fj]
def step(X,rhs):
 f=rhs(X);r=[h*max(abs(v[0]),abs(v[1]))+Q('.00001') for v in f]
 for _ in range(32):
  Y=[inflate(x,v) for x,v in zip(X,r)];g=rhs(Y)
  if all(contains(y,add(x,sc(fi,h))) for x,y,fi in zip(X,Y,g)):break
  r=[2*v for v in r]
 else:raise RuntimeError('Picard inclusion failed')
 return [dyad(add(x,sc(fi,h))) for x,fi in zip(X,g)]
# Fixed modal frame vs time-varying shear frame. Both start with the SAME physical set.
fixed=[P,M,Z];s=F(0);moving=[P,M,Z]
rows=[]
for i in range(81):
 if i in (0,10,20,40,60,80):
  pp=add(moving[0],sc(moving[1],s))
  rows.append(dict(horizon=float(i*h),fixed_width=float(2*(width(fixed[0])+width(fixed[1]))),shear_width=float(2*(width(pp)+width(moving[1]))),fixed_J_upper=float(fixed[2][1]),shear_J_upper=float(moving[2][1]),shear=float(s)))
 if i==80:break
 fixed=step(fixed,modal)
 moving=step(moving,lambda Y:rhs_shear(Y,s))
 # Shear follows local Jacobian: dp/dm = -6 p m, dm/dm = a-2k-3p^2-3m^2.
 # Conservative rational midpoint selection is an algorithmic heuristic, not a proof of optimality.
 p_mid=(moving[0][0]+moving[0][1])/2+s*(moving[1][0]+moving[1][1])/2
 m_mid=(moving[1][0]+moving[1][1])/2
 alpha=Q('.17')-3*p_mid*p_mid-3*m_mid*m_mid
 beta=Q('.17')-2*k-3*p_mid*p_mid-3*m_mid*m_mid
 coupling=-6*p_mid*m_mid
 s_next=s+h*(coupling+(alpha-beta)*s)
 # Change coordinates exactly: v_new=v_old+(s_old-s_new)*m.
 moving[0]=dyad(add(moving[0],sc(moving[1],s-s_next)))
 s=s_next
with open('/mnt/data/matrix_frame_results.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
for r in rows:print(r)
