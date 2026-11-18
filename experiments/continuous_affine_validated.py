"""Exact-rational, Picard-validated continuous-time enclosure for a real
invariant Stuart--Landau subsystem and accumulated actuator cost.
All arithmetic and inclusion tests use fractions.Fraction.
"""
from fractions import Fraction as Q
import csv
h=Q('0.025'); a0=Q('0.17'); d=Q('0.02'); k=Q('0.7'); ca=Q('0.1'); steps=4
Aparam=(a0-d,a0+d)
def add(x,y):return (x[0]+y[0],x[1]+y[1])
def neg(x):return (-x[1],-x[0])
def sub(x,y):return add(x,neg(y))
def mul(x,y):
 p=[u*v for u in x for v in y];return min(p),max(p)
def scale(x,s):return (min(s*x[0],s*x[1]),max(s*x[0],s*x[1]))
def point(x):return x,x
def sym(r):return -r,r
def cube(x):return x[0]**3,x[1]**3
def sq(x):return (Q(0) if x[0]<=0<=x[1] else min(x[0]**2,x[1]**2),max(x[0]**2,x[1]**2))
def absmax(x):return max(abs(x[0]),abs(x[1]))
def contains(x,y):return x[0]<=y[0] and y[1]<=x[1]
def fb(X):
 u,v,_=X; dif=sub(u,v)
 return [sub(sub(mul(Aparam,u),cube(u)),scale(dif,k)),sub(sub(mul(Aparam,v),cube(v)),scale(neg(dif),k)),scale(sq(dif),2*ca*k*k)]
def jacb(X):
 u,v,_=X; dif=sub(u,v)
 return [[sub(sub(Aparam,scale(sq(u),3)),point(k)),point(k),point(Q(0))],
         [point(k),sub(sub(Aparam,scale(sq(v),3)),point(k)),point(Q(0))],
         [scale(dif,4*ca*k*k),scale(dif,-4*ca*k*k),point(Q(0))]]
def fpoint(c):
 u,v,j=c
 return [a0*u-u**3-k*(u-v),a0*v-v**3-k*(v-u),2*ca*k*k*(u-v)**2]
def jacpoint(c):
 u,v,_=c
 return [[a0-3*u*u-k,k,Q(0)],[k,a0-3*v*v-k,Q(0)],[4*ca*k*k*(u-v),-4*ca*k*k*(u-v),Q(0)]]
def local_defect(X):
 # Certified tube, inflated uniformly, verified via Picard inclusion.
 # All values exact rationals. Inflate to accommodate the entire trajectory.
 F=fb(X); radii=[h*(absmax(v)+Q('0.2')) for v in F]
 for _ in range(20):
  Y=[add(x,sym(r)) for x,r in zip(X,radii)]
  FY=fb(Y)
  if all(contains(y,add(x,scale(f,h))) for x,y,f in zip(X,Y,FY)):
   break
  radii=[r*2 for r in radii]
 else:raise RuntimeError('Picard tube could not be verified')
 # For every point in tube, y''=Df(y)f(y).
 D=jacb(Y)
 acceleration=[point(Q(0)) for _ in range(3)]
 for i in range(3):
  for j in range(3):acceleration[i]=add(acceleration[i],mul(D[i][j],FY[j]))
 return [h*h*absmax(v)/2 for v in acceleration],Y
c=[Q('0.45'),Q('0.21'),Q(0)]; sens=[Q(0)]*3; r=[Q(0)]*3
box=[point(v) for v in c]
rows=[]
for n in range(steps+1):
 affine=[add(point(c[i]),sym(abs(sens[i])+r[i])) for i in range(3)]
 rows.append([n,float(n*h),float(sum(x[1]-x[0] for x in affine[:2])),float(sum(x[1]-x[0] for x in box[:2])),float(affine[2][1]),float(box[2][1]),float(max(r))])
 if n==steps:break
 defect,Y=local_defect(affine)
 f=fpoint(c); J=jacpoint(c)
 U=[abs(sens[i])+r[i] for i in range(3)]
 rem=[d*U[i]+3*abs(c[i])*U[i]**2+U[i]**3 for i in range(2)]
 rem.append(2*ca*k*k*(U[0]+U[1])**2)
 newr=[r[i]+h*(sum(abs(J[i][j])*r[j] for j in range(3))+rem[i])+defect[i] for i in range(3)]
 newsens=[sens[i]+h*(sum(J[i][j]*sens[j] for j in range(3))+(d*c[i] if i<2 else Q(0))) for i in range(3)]
 newc=[c[i]+h*f[i] for i in range(3)]
 boxdef,_=local_defect(box)
 F=fb(box)
 box=[add(add(box[i],scale(F[i],h)),sym(boxdef[i])) for i in range(3)]
 c,sens,r=newc,newsens,newr
with open('/mnt/data/continuous_affine_results.csv','w',newline='') as o:
 w=csv.writer(o);w.writerow(['step','time','affine_state_total_width','box_state_total_width','affine_J_upper','box_J_upper','affine_remainder_max']);w.writerows(rows)
for i in (0,1,2,3,4):print(rows[i])
