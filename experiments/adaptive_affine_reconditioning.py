from fractions import Fraction as Q
import csv
h=Q('0.025'); a0=Q('0.17'); d=Q('0.02'); k=Q('0.7'); ca=Q('0.1'); E0=Q('0.1'); floor=Q('0.08'); supply=Q('0.01'); basal=Q('0.005')
A=(a0-d,a0+d)
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[1],-x[0]
def sub(x,y):return add(x,neg(y))
def mul(x,y):
 p=[u*v for u in x for v in y];return min(p),max(p)
def sc(x,s):return mul(x,(s,s))
def pt(v):return v,v
def sym(r):return -r,r
def cube(x):return x[0]**3,x[1]**3
def sq(x):return (Q(0) if x[0]<=0<=x[1] else min(x[0]**2,x[1]**2),max(x[0]**2,x[1]**2))
def mag(x):return max(abs(x[0]),abs(x[1]))
def contains(x,y):return x[0]<=y[0] and y[1]<=x[1]
def dyad(x,b):
 D=1<<b;return Q((x[0]*D).__floor__(),D),Q((x[1]*D).__ceil__(),D)
def rhs(X):
 u,v,J=X; dif=sub(u,v)
 return [sub(sub(mul(A,u),cube(u)),sc(dif,k)),sub(sub(mul(A,v),cube(v)),sc(neg(dif),k)),sc(sq(dif),2*ca*k*k)]
def jac(X):
 u,v,_=X; dif=sub(u,v)
 return [[sub(sub(A,sc(sq(u),3)),pt(k)),pt(k),pt(Q(0))],[pt(k),sub(sub(A,sc(sq(v),3)),pt(k)),pt(Q(0))],[sc(dif,4*ca*k*k),sc(dif,-4*ca*k*k),pt(Q(0))]]
def picard_defect(X):
 F=rhs(X);rad=[h*(mag(v)+Q('0.2')) for v in F]
 for it in range(35):
  Y=[add(x,sym(r)) for x,r in zip(X,rad)]; FY=rhs(Y)
  if all(contains(y,add(x,sc(f,h))) for x,y,f in zip(X,Y,FY)):break
  rad=[2*r for r in rad]
 else:raise RuntimeError('Picard failure')
 D=jac(Y); accel=[]
 for i in range(3):
  s=pt(Q(0))
  for j in range(3):s=add(s,mul(D[i][j],FY[j]))
  accel.append(h*h*mag(s)/2)
 return accel,Y

def fpoint(c):
 u,v,_=c;return [a0*u-u**3-k*(u-v),a0*v-v**3-k*(v-u),2*ca*k*k*(u-v)**2]
def jpoint(c):
 u,v,_=c;return [[a0-3*u*u-k,k,Q(0)],[k,a0-3*v*v-k,Q(0)],[4*ca*k*k*(u-v),-4*ca*k*k*(u-v),Q(0)]]
def run(bits,steps,center=Q("0.17"),radius=Q("0.02")):
 global a0,d,A
 a0,d=center,radius;A=(a0-d,a0+d)
 c=[Q('0.45'),Q('0.21'),Q(0)]; s=[Q(0)]*3;r=[Q(0)]*3;box=[pt(v) for v in c]
 for n in range(steps):
  aff=[add(pt(c[i]),sym(abs(s[i])+r[i])) for i in range(3)]
  defect,Y=picard_defect(aff)
  f=fpoint(c);J=jpoint(c);U=[abs(s[i])+r[i] for i in range(3)]
  rem=[d*U[i]+3*abs(c[i])*U[i]**2+U[i]**3 for i in range(2)]
  rem.append(2*ca*k*k*(U[0]+U[1])**2)
  nr=[r[i]+h*(sum(abs(J[i][j])*r[j] for j in range(3))+rem[i])+defect[i] for i in range(3)]
  ns=[s[i]+h*(sum(J[i][j]*s[j] for j in range(3))+(d*c[i] if i<2 else Q(0))) for i in range(3)]
  nc=[c[i]+h*f[i] for i in range(3)]
  # Round affine coefficients outward via enclosing coefficient intervals.
  # Center and sensitivity rounding errors are absorbed in remainder.
  for i in range(3):
   cc=dyad(pt(nc[i]),bits);ss=dyad(pt(ns[i]),bits)
   nr[i]+=max(abs(cc[0]-nc[i]),abs(cc[1]-nc[i]))+max(abs(ss[0]-ns[i]),abs(ss[1]-ns[i]))
   nc[i]=(cc[0]+cc[1])/2;ns[i]=(ss[0]+ss[1])/2
   nr[i]=dyad(pt(nr[i]),bits)[1]
  bd,_=picard_defect(box);F=rhs(box)
  box=[dyad(add(add(box[i],sc(F[i],h)),sym(bd[i])),bits) for i in range(3)]
  c,s,r=nc,ns,nr
 aff=[add(pt(c[i]),sym(abs(s[i])+r[i])) for i in range(3)]
 H=steps*h;energy=E0+(supply-basal)*H-aff[2][1]
 # E(t) lower bounded by E0 - basal*t - J(t), as J is nondecreasing and supply>=0;
 # to certify every intermediate time use E0-basal*H - J_upper.
 prefix_lower=E0-basal*H-aff[2][1]
 return [bits,float(H),float(sum(x[1]-x[0] for x in box[:2])),float(sum(x[1]-x[0] for x in aff[:2])),float(box[2][1]),float(aff[2][1]),float(energy),float(prefix_lower),bool(prefix_lower>=floor)]

# Parameter-cell recentering: each subcell carries its own nominal trajectory and
# sensitivity, with its original parameter value held fixed for the full horizon.
# This is static cell-centered reconditioning, NOT a Lohner state-frame algorithm.
rows=[]
for steps in (4,8,12,16):
 for cells in (1,2,4):
  reports=[]
  for j in range(cells):
   center=Q('0.15')+Q('0.04')*Q(2*j+1,2*cells)
   radius=Q('0.02')/cells
   reports.append(run(40,steps,center,radius))
  row={'horizon':float(h*steps),'cells':cells,'max_cell_affine_width':max(r[3] for r in reports),'max_J_upper':max(r[5] for r in reports),'minimum_prefix_energy_lower':min(r[7] for r in reports),'all_cells_floor_certified':all(r[8] for r in reports),'max_direct_width':max(r[2] for r in reports)}
  rows.append(row);print(row,flush=True)
with open('/mnt/data/adaptive_affine_reconditioning_results.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
