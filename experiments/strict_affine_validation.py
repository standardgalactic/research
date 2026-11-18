"""Directed-interval Picard/Taylor verification with no float endpoint tests.
Reuses the earlier polynomial vector field, but all verification decisions are
performed on interval endpoints, and time steps are exact decimal inputs.
"""
import importlib.util, pathlib, ast, csv
from decimal import Decimal
from mpmath import iv
iv.dps=45
source=pathlib.Path('/mnt/data/affine_nonlinear_parameter.py').read_text()
# Load definitions without running its experiment loop.
source=source.split('rows=[]')[0]
ns={};exec(compile(source,'prior_core','exec'),ns)
P=ns['P']; rhs=ns['augmented_rhs']; polynomial_rhs=ns['polynomial_rhs']; z0=ns['z0']

def strict_inside(v,w): return bool(v.a>w.a and v.b<w.b)
def width(v): return float(v.b-v.a)
def upper(v): return float(v.b)
def mag(v): return max(abs(float(v.a)),abs(float(v.b)))
def expand(x,r):return x+iv.mpf([str(-r),str(r)])
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

def step(x,a,h,augmented):
 f=rhs(x,a) if augmented else rhs(x+[iv.mpf(0)]*14,a)[:7]
 hh=iv.mpf(h);size=Decimal('0.00000001')+Decimal(h)*Decimal('2')*Decimal(str(max(mag(v) for v in f)))
 for trial in range(40):
  y=[expand(v,size) for v in x]
  fy=rhs(y,a) if augmented else rhs(y+[iv.mpf(0)]*14,a)[:7]
  tube=[v+iv.mpf(['0',h])*w for v,w in zip(x,fy)]
  if all(strict_inside(v,w) for v,w in zip(tube,y)):break
  size*=Decimal('1.7')
 else:raise RuntimeError('Picard inclusion failed')
 if augmented:
  dy=[D(v,d) for v,d in zip(y,fy)]
  zp=[P([dy[i],dy[i+7],dy[i+14]]) for i in range(7)]
  fp=polynomial_rhs(zp,P([D(a),D(1),D(0)]))
  second=[fp[i].c[j].d for j in range(3) for i in range(7)]
 else:
  zp=[P([D(y[i],fy[i]),0,0]) for i in range(7)]
  second=[v.c[0].d for v in polynomial_rhs(zp,P([D(a),0,0]))]
 return [v+hh*w+(hh*hh/2)*d for v,w,d in zip(x,f,second)]

def project(y,d,wide):
 u=iv.mpf([str(-d),str(d)]);u2=iv.mpf(['0',str(d*d)])
 return [y[i]+u*y[i+7]+u2*wide[i+14] for i in range(7)]
rows=[]
for n in (4,8):
 for horizon in ('0.05','0.1','0.2'):
  aw=dw=aj=dj=0.0
  for j in range(n):
   lo=Decimal('0.12')+Decimal('0.10')*j/n;hi=Decimal('0.12')+Decimal('0.10')*(j+1)/n
   m=(lo+hi)/2;d=(hi-lo)/2
   x=[iv.mpf(v) for v in z0];y=x+[iv.mpf(0)]*14;wide=list(y)
   a=iv.mpf([str(lo),str(hi)])
   for _ in range(int(Decimal(horizon)/Decimal('0.005'))):
    x=step(x,a,'0.005',False)
    y=step(y,iv.mpf(str(m)),'0.005',True)
    wide=step(wide,a,'0.005',True)
   p=project(y,d,wide)
   aw=max(aw,max(width(v) for v in p[:6]));dw=max(dw,max(width(v) for v in x[:6]))
   aj=max(aj,upper(p[6]));dj=max(dj,upper(x[6]))
  rows.append([n,horizon,dw,aw,dj,aj,100*(1-aw/dw)])
  print(rows[-1],flush=True)
with open('/mnt/data/strict_affine_results.csv','w',newline='') as f:
 w=csv.writer(f);w.writerow(['cells','horizon','direct_state_width','affine_state_width','direct_J_upper','affine_J_upper','state_width_reduction_pct']);w.writerows(rows)
