"""Interval Taylor-2 with validated Picard tube, fixed parameter cells.
Requires mpmath; all arithmetic uses mpmath.iv with directed interval rounding.
"""
from mpmath import iv
from decimal import Decimal
import csv
iv.dps=35
L=((1,-1,0),(-1,2,-1),(0,-1,1)); nu=(-1,0,1)
k=iv.mpf('0.7'); s=iv.mpf('0.16'); ca=iv.mpf('0.1')
z0=['0.45','0.21','0.33','0.10','-0.06','0.04','0']
def bounds(x):return float(x.a),float(x.b)
def rad(x):return max(abs(bounds(x)[0]),abs(bounds(x)[1]))
def hull(x,r):return iv.mpf([bounds(x)[0]-r,bounds(x)[1]+r])
def subset(x,y):return bounds(x)[0]>bounds(y)[0] and bounds(x)[1]<bounds(y)[1]
def width(x):return bounds(x)[1]-bounds(x)[0]
def rhs(x,a):
    z=x[:6];o=[]
    for i in range(3):
        r=z[i];v=z[i+3]; m=r*r+v*v
        o.append(a*r-s*nu[i]*v-m*r-k*sum(L[i][j]*z[j] for j in range(3)))
    for i in range(3):
        r=z[i];v=z[i+3];m=r*r+v*v
        o.append(a*v+s*nu[i]*r-m*v-k*sum(L[i][j]*z[j+3] for j in range(3)))
    lap=[sum(z[j+off]*L[i][j] for j in range(3)) for off in (0,3) for i in range(3)]
    o.append(ca*k*k*sum(v*v for v in lap))
    return o
# interval truncated jets, coefficients x^(j)/j!; evaluates polynomial RHS exactly via jet arithmetic
class Jet:
    def __init__(self,c):self.c=c
    def __add__(self,b):
        if not isinstance(b,Jet):b=Jet([b,0,0])
        return Jet([self.c[i]+b.c[i] for i in range(3)])
    __radd__=__add__
    def __neg__(self):return Jet([-x for x in self.c])
    def __sub__(self,b):return self+-b if isinstance(b,Jet) else self+(-b)
    def __rsub__(self,b):return -self+b
    def __mul__(self,b):
        if not isinstance(b,Jet):b=Jet([b,0,0])
        return Jet([sum(self.c[j]*b.c[i-j] for j in range(i+1)) for i in range(3)])
    __rmul__=__mul__
def rhsjet(x,a):
    z=x[:6];o=[]
    for i in range(3):
        r=z[i];v=z[i+3];m=r*r+v*v
        o.append(r*a-v*(s*nu[i])-m*r-sum((z[j]*L[i][j] for j in range(3)),Jet([0,0,0]))*k)
    for i in range(3):
        r=z[i];v=z[i+3];m=r*r+v*v
        o.append(v*a+r*(s*nu[i])-m*v-sum((z[j+3]*L[i][j] for j in range(3)),Jet([0,0,0]))*k)
    lap=[sum((z[j+off]*L[i][j] for j in range(3)),Jet([0,0,0])) for off in (0,3) for i in range(3)]
    o.append(sum((v*v for v in lap),Jet([0,0,0]))*(ca*k*k))
    return o
def coefficients(x,a):
    jets=[Jet([v,iv.mpf(0),iv.mpf(0)]) for v in x]
    for order in (0,1):
        fj=rhsjet(jets,a)
        for i in range(7):jets[i].c[order+1]=fj[i].c[order]/(order+1)
    fj=rhsjet(jets,a)
    return [[jets[i].c[0],jets[i].c[1],jets[i].c[2],fj[i].c[2]/3] for i in range(7)]
def step(x,a,h,order):
    # Construct a Picard tube Y enclosing all trajectories for 0<=t<=h.
    f0=rhs(x,a)
    size=max(0.0001, 2*h*max(rad(v) for v in f0))
    for _ in range(35):
        y=[hull(v,size) for v in x]
        fy=rhs(y,a)
        test=[v+iv.mpf([0,str(h)])*q for v,q in zip(x,fy)]
        if all(subset(v,w) for v,w in zip(test,y)):break
        size*=1.6
    else:raise RuntimeError('Picard tube not verified')
    # Taylor remainder bounded by interval coefficients evaluated over the whole tube.
    coef0=coefficients(x,a);coefY=coefficients(y,a)
    hh=iv.mpf(str(h));out=[]
    for c,d in zip(coef0,coefY):
        if order==1:
            # second derivative coefficient from the full tube
            out.append(c[0]+hh*c[1]+hh**2*d[2])
        else:
            out.append(c[0]+hh*c[1]+hh**2*c[2]+hh**3*d[3])
    return out
rows=[]
for n in (4,8):
 for T in (0.05,0.1,0.2):
    for order in (1,2):
     maxw=0;maxj=0;successful=True
     try:
      for i in range(n):
       lo=Decimal('0.12')+Decimal('0.10')*Decimal(i)/Decimal(n)
       hi=Decimal('0.12')+Decimal('0.10')*Decimal(i+1)/Decimal(n)
       a=iv.mpf([str(lo),str(hi)])
       x=[iv.mpf(v) for v in z0]
       h=.005
       for j in range(round(T/h)):x=step(x,a,h,order)
       maxw=max(maxw,max(width(v) for v in x[:6]));maxj=max(maxj,bounds(x[6])[1])
     except (RuntimeError,OverflowError,ValueError) as e:successful=False;print('fail',n,T,order,e)
     rows.append((n,T,order,maxw,maxj,successful));print(rows[-1],flush=True)
with open('/mnt/data/nonlinear_taylor_results.csv','w',newline='') as f:
 w=csv.writer(f);w.writerow(['cells','horizon','order','max_state_width','max_expenditure_upper','all_picard_tubes_verified']);w.writerows(rows)
