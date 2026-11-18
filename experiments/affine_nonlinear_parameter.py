"""Validated parameter-polynomial (degree 1 + interval quadratic remainder) propagation.
mpmath.iv directed interval operations; local Picard tubes and Taylor-1 time steps.
"""
from mpmath import iv
from decimal import Decimal
import csv
iv.dps=35
L=((1,-1,0),(-1,2,-1),(0,-1,1)); nu=(-1,0,1)
k=iv.mpf('0.7'); s=iv.mpf('0.16'); ca=iv.mpf('0.1')
z0=['0.45','0.21','0.33','0.10','-0.06','0.04','0']
def ends(v):return float(v.a),float(v.b)
def wid(v):return ends(v)[1]-ends(v)[0]
def magnitude(v):return max(abs(q) for q in ends(v))
def expand(v,r):return v+iv.mpf([str(-r),str(r)])
def inside(v,w):return ends(v)[0]>ends(w)[0] and ends(v)[1]<ends(w)[1]
class P:
 def __init__(self,c):self.c=c
 def __add__(self,b):
  if not isinstance(b,P):b=P([b,0,0])
  return P([self.c[i]+b.c[i] for i in range(3)])
 __radd__=__add__
 def __neg__(self):return P([-x for x in self.c])
 def __sub__(self,b):return self+-b
 def __rsub__(self,b):return -self+b
 def __mul__(self,b):
  if not isinstance(b,P):b=P([b,0,0])
  return P([sum(self.c[j]*b.c[i-j] for j in range(i+1)) for i in range(3)])
 __rmul__=__mul__
def polynomial_rhs(z,a):
 out=[]
 for off in (0,3):
  for i in range(3):
   r=z[i];v=z[i+3];q=r*r+v*v
   out.append(a*z[i+off]+(r if off else -v)*(s*nu[i])-q*z[i+off]-sum((L[i][j]*z[j+off] for j in range(3)),P([0,0,0]))*k)
 lap=[sum((L[i][j]*z[j+off] for j in range(3)),P([0,0,0])) for off in (0,3) for i in range(3)]
 out.append(sum((v*v for v in lap),P([0,0,0]))*(ca*k*k))
 return out
def augmented_rhs(y,a):
 z=[P([y[i],y[i+7],y[i+14]]) for i in range(7)]
 p=polynomial_rhs(z,P([a,1,0]))
 return [p[i].c[j] for j in range(3) for i in range(7)]
def picard_step(x,a,h):
 f=augmented_rhs(x,a); size=max(1e-8,2*h*max(magnitude(v) for v in f))
 for _ in range(32):
  y=[expand(v,size) for v in x]
  fy=augmented_rhs(y,a)
  tube=[v+iv.mpf(['0',str(h)])*w for v,w in zip(x,fy)]
  if all(inside(v,w) for v,w in zip(tube,y)):break
  size*=1.7
 else:raise RuntimeError('Picard failure')
 # Integral Taylor 1: x+h*f(x)+h^2/2*f'(Y)f(Y); bound directional derivative
 # by symmetric interval differentiation of polynomial jet in time.
 fy=augmented_rhs(y,a)
 # Bound f'(Y)f(Y) using interval automatic differentiation in a dual number.
 class D:
  def __init__(self,v,d=0):self.v=v;self.d=d
  def __add__(self,b):
   if not isinstance(b,D):b=D(b)
   return D(self.v+b.v,self.d+b.d)
  __radd__=__add__
  def __neg__(self):return D(-self.v,-self.d)
  def __sub__(self,b):return self+-b
  def __rsub__(self,b):return -self+b
  def __mul__(self,b):
   if not isinstance(b,D):b=D(b)
   return D(self.v*b.v,self.d*b.v+self.v*b.d)
  __rmul__=__mul__
 dy=[D(v,d) for v,d in zip(y,fy)]
 # P coefficients themselves may be duals; polynomial_rhs handles both
 zp=[P([dy[i],dy[i+7],dy[i+14]]) for i in range(7)]
 fp=polynomial_rhs(zp,P([D(a),D(1),D(0)]))
 second=[fp[i].c[j].d for j in range(3) for i in range(7)]
 hh=iv.mpf(str(h))
 return [v+hh*w+hh*hh*iv.mpf('0.5')*d for v,w,d in zip(x,f,second)]
def direct_projection(y,d,remainder):
 # Taylor in parameter: z(m+u) = z(m)+u*z_a(m)+u^2*z_aa(m)/2
 # augmented components 14:21 store coefficient of u^2, so remainder d^2*T
 u=iv.mpf([str(-d),str(d)]);u2=iv.mpf([0,str(d*d)])
 return [y[i]+u*y[i+7]+u2*remainder[i+14] for i in range(7)]
def direct_rhs(x,a):
 return augmented_rhs(x+[iv.mpf(0)]*14,a)[:7]
def direct_step(x,a,h):
 f=direct_rhs(x,a);size=max(1e-8,2*h*max(magnitude(v) for v in f))
 for _ in range(32):
  y=[expand(v,size) for v in x];fy=direct_rhs(y,a)
  test=[v+iv.mpf(['0',str(h)])*w for v,w in zip(x,fy)]
  if all(inside(v,w) for v,w in zip(test,y)):break
  size*=1.7
 else:raise RuntimeError('direct Picard failure')
 hh=iv.mpf(str(h))
 return [v+hh*w+hh*hh*iv.mpf('0.5')*d for v,w,d in zip(x,f,second_direct(y,fy,a))]
def second_direct(y,fy,a):
 class D:
  def __init__(self,v,d=0):self.v=v;self.d=d
  def __add__(self,b):
   if not isinstance(b,D):b=D(b)
   return D(self.v+b.v,self.d+b.d)
  __radd__=__add__
  def __neg__(self):return D(-self.v,-self.d)
  def __sub__(self,b):return self+-b
  def __rsub__(self,b):return -self+b
  def __mul__(self,b):
   if not isinstance(b,D):b=D(b)
   return D(self.v*b.v,self.d*b.v+self.v*b.d)
  __rmul__=__mul__
 z=[P([D(y[i],fy[i]),0,0]) for i in range(7)]
 return [v.c[0].d for v in polynomial_rhs(z,P([D(a),0,0]))]
rows=[]
for n in (4,8):
 for horizon in ((0.05,0.1,0.2) if n==4 else (0.05,0.1)):
  max_affine=0;max_direct=0;max_ja=0;max_jd=0
  for j in range(n):
   lo=Decimal('0.12')+Decimal('0.10')*j/n;hi=Decimal('0.12')+Decimal('0.10')*(j+1)/n
   m=(lo+hi)/2;d=(hi-lo)/2
   x=[iv.mpf(v) for v in z0];y=x+[iv.mpf(0)]*14;yw=x+[iv.mpf(0)]*14
   for i in range(7):y[i+7]=iv.mpf(0)
   for i in range(3):y[i+7]=iv.mpf(z0[i])
   for i in range(3,6):y[i+7]=iv.mpf(z0[i])
   # Initial state independent of a: first and second sensitivities zero
   y=x+[iv.mpf(0)]*14
   a=iv.mpf([str(lo),str(hi)])
   for _ in range(round(horizon/.005)):
    x=direct_step(x,a,.005);y=picard_step(y,iv.mpf(str(m)),.005);yw=picard_step(yw,a,.005)
   p=direct_projection(y,d,yw)
   max_affine=max(max_affine,max(wid(v) for v in p[:6]));max_direct=max(max_direct,max(wid(v) for v in x[:6]))
   max_ja=max(max_ja,ends(p[6])[1]);max_jd=max(max_jd,ends(x[6])[1])
  rows.append((n,horizon,max_direct,max_affine,max_jd,max_ja,max_direct/max_affine))
  print(rows[-1],flush=True)
with open('/mnt/data/affine_nonlinear_results.csv','w',newline='') as f:
 w=csv.writer(f);w.writerow(['cells','horizon','direct_max_width','parameter_affine_max_width','direct_max_J','parameter_affine_max_J','width_improvement']);w.writerows(rows)
