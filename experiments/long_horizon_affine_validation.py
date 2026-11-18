from fractions import Fraction as F
import csv
Q=lambda x:F(str(x))
h=Q('.01');k=Q('.7');ca=Q('.1');bits=40; a= (Q('.15'),Q('.19'))
def I(x,y=None):return (x,x if y is None else y)
def add(x,y):return x[0]+y[0],x[1]+y[1]
def sub(x,y):return x[0]-y[1],x[1]-y[0]
def mul(x,y):
 v=[t*s for t in x for s in y];return min(v),max(v)
def sc(x,v):return mul(x,I(v))
def sq(x):return (F(0) if x[0]<=0<=x[1] else min(x[0]**2,x[1]**2),max(x[0]**2,x[1]**2))
def cube(x):return x[0]**3,x[1]**3
def dyad(x):
 d=1<<bits;return F((x[0]*d).__floor__(),d),F((x[1]*d).__ceil__(),d)
def absmax(x):return max(abs(x[0]),abs(x[1]))
def hull(c,gs,r):return [I(c[i]-sum(abs(g[i]) for g in gs)-r[i],c[i]+sum(abs(g[i]) for g in gs)+r[i]) for i in range(3)]
def rhs(x,aa):
 p,m,j=x
 return [sub(sub(mul(aa,p),cube(p)),sc(mul(p,sq(m)),3)),sub(sub(mul(sub(aa,I(2*k)),m),sc(mul(sq(p),m),3)),cube(m)),sc(sq(m),8*ca*k*k)]
def jac(c,aa):
 p,m,j=c
 return [[aa-3*p*p-3*m*m,-6*p*m,F(0)],[-6*p*m,aa-2*k-3*p*p-3*m*m,F(0)],[F(0),16*ca*k*k*m,F(0)]]
def jac_interval(X,aa):
 p,m,j=X
 return [[sub(sub(aa,sc(sq(p),3)),sc(sq(m),3)),sc(mul(p,m),-6),I(0)], [sc(mul(p,m),-6),sub(sub(sub(aa,I(2*k)),sc(sq(p),3)),sc(sq(m),3)),I(0)],[I(0),sc(m,16*ca*k*k),I(0)]]
def tube(X):
 # A priori interval T with X+h*f(T) strictly contained in its interior.
 T=[I(x[0]-Q('.01'),x[1]+Q('.01')) for x in X]
 for attempt in range(30):
  f=rhs(T,a); Y=[add(x,sc(v,h)) for x,v in zip(X,f)]
  if all(T[i][0]<Y[i][0] and Y[i][1]<T[i][1] for i in range(3)):return T
  T=[I(min(t[0],y[0])-Q('.001'),max(t[1],y[1])+Q('.001')) for t,y in zip(T,Y)]
 raise RuntimeError('Picard tube inclusion failed')
def defect(T):
 # x''=D_x f(x,a) f(x,a), a is fixed in time
 f=rhs(T,a);G=jac_interval(T,a)
 return [h*h*sum(absmax(mul(G[i][j],f[j])) for j in range(3))/2 for i in range(3)]
def step_aff(c,gs,r):
 aa=Q('.17');G=jac(c,aa);fc=[v[0] for v in rhs([I(v) for v in c],I(aa))]
 newc=[c[i]+h*fc[i] for i in range(3)]
 newgs=[]
 for idx,g in enumerate(gs):
  dg=[sum(G[i][j]*g[j] for j in range(3)) for i in range(3)]
  if idx==2:dg=[dg[i]+(c[i]*Q('.02') if i<2 else 0) for i in range(3)]
  newgs.append([g[i]+h*dg[i] for i in range(3)])
 X=hull(c,gs,r);T=tube(X);D=defect(T);JX=jac_interval(X,a)
 rad=[sum(abs(g[i]) for g in gs)+r[i] for i in range(3)]
 rr=[]
 for i in range(3):
  diff=sum(max(abs(JX[i][j][0]-G[i][j]),abs(JX[i][j][1]-G[i][j]))*rad[j] for j in range(3))
  if i<2:diff+=Q('.02')*rad[i]
  transported=sum(abs(G[i][j])*r[j] for j in range(3))
  rr.append(r[i]+h*(transported+diff)+D[i])
 for i in range(3):
  cc=dyad(I(newc[i]));rr[i]+=max(abs(cc[0]-newc[i]),abs(cc[1]-newc[i]));newc[i]=(cc[0]+cc[1])/2
  for g in newgs:
   gi=dyad(I(g[i]));rr[i]+=max(abs(gi[0]-g[i]),abs(gi[1]-g[i]));g[i]=(gi[0]+gi[1])/2
  rr[i]=dyad(I(rr[i]))[1]
 return newc,newgs,rr,T

def step_box(X):
 T=tube(X);D=defect(T)
 f=rhs(X,a)
 return [dyad(add(add(x,sc(f[i],h)),I(-D[i],D[i]))) for i,x in enumerate(X)],T
c=[Q('.33'),Q('.12'),F(0)];gs=[[Q('.001'),F(0),F(0)],[F(0),Q('.001'),F(0)],[F(0),F(0),F(0)]];r=[F(0)]*3
box=hull(c,gs,r);rows=[];floor=Q('.08');E0=Q('.1');q=Q('.01');b=Q('.005')
for n in range(201):
 if n in (0,20,40,80,120,160,200):
  A=hull(c,gs,r);t=n*h
  rows.append(dict(time=float(t),box_state_width=float(sum(x[1]-x[0] for x in box[:2])),affine_state_width=float(sum(x[1]-x[0] for x in A[:2])),box_J_upper=float(box[2][1]),affine_J_upper=float(A[2][1]),affine_energy_lower=float(E0+(q-b)*t-A[2][1]),floor_pass=(E0+(q-b)*t-A[2][1]>=floor)))
 if n<200:
  box,_=step_box(box)
  c,gs,r,_=step_aff(c,gs,r)
with open('/mnt/data/long_horizon_affine_results.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
for row in rows:print(row)
