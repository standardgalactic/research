from fractions import Fraction as F
import csv
Q=lambda s:F(str(s))
h=Q('.01'); k=Q('.7'); ca=Q('.1'); bits=40

def I(x,y=None):return (x,x if y is None else y)
def add(x,y):return x[0]+y[0],x[1]+y[1]
def sub(x,y):return x[0]-y[1],x[1]-y[0]
def mul(x,y):
 v=[a*b for a in x for b in y];return min(v),max(v)
def sc(x,v):return mul(x,I(v))
def sq(x):return (F(0) if x[0]<=0<=x[1] else min(x[0]**2,x[1]**2),max(x[0]**2,x[1]**2))
def cube(x):return x[0]**3,x[1]**3
def dyad(x):
 d=1<<bits;return F((x[0]*d).__floor__(),d),F((x[1]*d).__ceil__(),d)
def hull(c,gs,r):
 return [add(I(c[i]),I(-sum(abs(g[i]) for g in gs)-r[i],sum(abs(g[i]) for g in gs)+r[i])) for i in range(3)]
def rhs(x,a):
 p,m,j=x
 return [sub(sub(mul(a,p),cube(p)),sc(mul(p,sq(m)),3)),sub(sub(mul(sub(a,I(2*k)),m),sc(mul(sq(p),m),3)),cube(m)),sc(sq(m),8*ca*k*k)]
def jac(c,a):
 p,m,j=c
 return [[a-3*p*p-3*m*m,-6*p*m,0],[-6*p*m,a-2*k-3*p*p-3*m*m,0],[0,16*ca*k*k*m,0]]
def matvec(A,v):return [sum(A[i][j]*v[j] for j in range(3)) for i in range(3)]
def step_aff(c,gs,r):
 # xi0,xi1 are independent initial-state generators, xi2 fixed uncertain growth parameter.
 a0=Q('.17'); G=jac(c,a0); fc=[v[0] for v in rhs([I(v) for v in c],I(a0))]
 newc=[c[i]+h*fc[i] for i in range(3)]
 newgs=[]
 for idx,g in enumerate(gs):
  dg=matvec(G,g)
  if idx==2:dg=[dg[i]+(c[i] if i<2 else 0)*Q('.02') for i in range(3)]
  newgs.append([g[i]+h*dg[i] for i in range(3)])
 X=hull(c,gs,r); actual=rhs(X,I(Q('.15'),Q('.19')))
 # residual evaluated as correlated polynomial minus affine model, using a safe interval extension.
 # subtract linear terms as intervals; enclosure is sound but can overestimate from dependency loss.
 linear=[]
 for i in range(3):
  z=I(fc[i])
  for g in gs:z=add(z,I(-abs(sum(G[i][j]*g[j] for j in range(3))),abs(sum(G[i][j]*g[j] for j in range(3)))))
  z=add(z,I(-sum(abs(G[i][j])*r[j] for j in range(3)),sum(abs(G[i][j])*r[j] for j in range(3))))
  if i<2:z=add(z,I(-Q('.02')*max(abs(X[i][0]),abs(X[i][1])),Q('.02')*max(abs(X[i][0]),abs(X[i][1]))))
  linear.append(z)
 # Instead of naive interval subtraction (which doubles the linear term), use an analytic Hessian remainder.
 # f polynomial in (p,m,a); bound second-order terms by Jacobian variation over the full enclosure.
 rad=[sum(abs(g[i]) for g in gs)+r[i] for i in range(3)]
 pmax=max(abs(X[0][0]),abs(X[0][1]));mmax=max(abs(X[1][0]),abs(X[1][1]));
 # Jacobian deviation from center bounded by interval evaluation of exact Jacobian over X and a interval.
 JX=[[sub(sub(I(Q('.15'),Q('.19')),sc(sq(X[0]),3)),sc(sq(X[1]),3)),sc(mul(X[0],X[1]),-6),I(0)],
     [sc(mul(X[0],X[1]),-6),sub(sub(I(Q('.15')-2*k,Q('.19')-2*k),sc(sq(X[0]),3)),sc(sq(X[1]),3)),I(0)],
     [I(0),sc(X[1],16*ca*k*k),I(0)]]
 # Variational remainder due to state deviation, plus parameter/state product.
 # Uses integral mean-value theorem, and is guaranteed for the discrete map.
 rr=[]
 for i in range(3):
  diff=sum(max(abs(JX[i][j][0]-G[i][j]),abs(JX[i][j][1]-G[i][j]))*rad[j] for j in range(3))
  if i<2:diff+=Q('.02')*rad[i]
  # Old remainder transport by G plus new nonlinear remainder
  transported=sum(abs(G[i][j])*r[j] for j in range(3))
  rr.append(r[i]+h*(transported+diff))
 # round center, generator, and remainder outward by bounding rounding errors in r
 for i in range(3):
  cc=dyad(I(newc[i])); err=max(abs(cc[0]-newc[i]),abs(cc[1]-newc[i]));newc[i]=(cc[0]+cc[1])/2;rr[i]+=err
  for g in newgs:
   gi=dyad(I(g[i])); err=max(abs(gi[0]-g[i]),abs(gi[1]-g[i]));g[i]=(gi[0]+gi[1])/2;rr[i]+=err
  rr[i]=dyad(I(rr[i]))[1]
 return newc,newgs,rr

def step_box(X):return [dyad(add(x,sc(f,h))) for x,f in zip(X,rhs(X,I(Q('.15'),Q('.19'))))]
c=[Q('.33'),Q('.12'),F(0)]
gs=[[Q('.001'),F(0),F(0)],[F(0),Q('.001'),F(0)],[F(0),F(0),F(0)]]
r=[F(0)]*3
box=hull(c,gs,r)
rows=[]
for n in range(41):
 if n in (0,5,10,20,30,40):
  A=hull(c,gs,r)
  rows.append({'time':float(n*h),'box_state_width':float(sum(x[1]-x[0] for x in box[:2])),'affine_state_width':float(sum(x[1]-x[0] for x in A[:2])),'box_J_upper':float(box[2][1]),'affine_J_upper':float(A[2][1]),'affine_remainder_p':float(r[0])})
 if n<40:
  box=step_box(box)
  c,gs,r=step_aff(c,gs,r)
with open('/mnt/data/lohner_generator_results.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
for row in rows:print(row)
