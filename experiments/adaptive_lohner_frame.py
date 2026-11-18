exec(open('/mnt/data/long_horizon_affine_validation.py').read().split('c=[Q(')[0])
import csv

def matvec(M,v):return [sum(M[i][j]*v[j] for j in range(3)) for i in range(3)]
def mm(A,B):return [[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
def inv3(A):
 a=[list(A[i])+[F(i==j) for j in range(3)] for i in range(3)]
 for j in range(3):
  pivot=next(i for i in range(j,3) if a[i][j]);a[j],a[pivot]=a[pivot],a[j]
  v=a[j][j];a[j]=[z/v for z in a[j]]
  for i in range(3):
   if i!=j:
    v=a[i][j];a[i]=[a[i][k]-v*a[j][k] for k in range(6)]
 return [row[3:] for row in a]
def hull_frame(c,gs,B,r):
 return [I(c[i]-sum(abs(g[i]) for g in gs)-sum(abs(B[i][j])*r[j] for j in range(3)),c[i]+sum(abs(g[i]) for g in gs)+sum(abs(B[i][j])*r[j] for j in range(3))) for i in range(3)]
def step_frame(c,gs,B,r):
 aa=Q('.17');G=jac(c,aa);fc=[v[0] for v in rhs([I(v) for v in c],I(aa))]
 newc=[c[i]+h*fc[i] for i in range(3)]
 newgs=[]
 for idx,g in enumerate(gs):
  dg=matvec(G,g)
  if idx==2:dg=[dg[i]+(c[i]*Q('.02') if i<2 else 0) for i in range(3)]
  newgs.append([g[i]+h*dg[i] for i in range(3)])
 X=hull_frame(c,gs,B,r);T=tube(X);D=defect(T);JX=jac_interval(X,a)
 rad=[sum(abs(g[i]) for g in gs)+sum(abs(B[i][j])*r[j] for j in range(3)) for i in range(3)]
 forcing=[]
 for i in range(3):
  diff=sum(max(abs(JX[i][j][0]-G[i][j]),abs(JX[i][j][1]-G[i][j]))*rad[j] for j in range(3))
  if i<2:diff+=Q('.02')*rad[i]
  forcing.append(h*diff+D[i])
 M=[[F(i==j)+h*G[i][j] for j in range(3)] for i in range(3)]
 newB=mm(M,B)
 # Inverse frame transforms symmetric forcing into a new remainder box.
 inv=inv3(newB)
 newr=[r[i]+sum(abs(inv[i][j])*forcing[j] for j in range(3)) for i in range(3)]
 # Outward dyadic rounding: account for center and generator rounding as physical forcing.
 rounderr=[F(0)]*3
 for i in range(3):
  cc=dyad(I(newc[i]));rounderr[i]+=max(abs(cc[0]-newc[i]),abs(cc[1]-newc[i]));newc[i]=(cc[0]+cc[1])/2
  for g in newgs:
   gi=dyad(I(g[i]));rounderr[i]+=max(abs(gi[0]-g[i]),abs(gi[1]-g[i]));g[i]=(gi[0]+gi[1])/2
  for j in range(3):
   bj=dyad(I(newB[i][j]));rounderr[i]+=max(abs(bj[0]-newB[i][j]),abs(bj[1]-newB[i][j]))*newr[j];newB[i][j]=(bj[0]+bj[1])/2
 inv=inv3(newB)
 newr=[dyad(I(newr[i]+sum(abs(inv[i][j])*rounderr[j] for j in range(3))))[1] for i in range(3)]
 return newc,newgs,newB,newr

if __name__=='__main__':
 c=[Q('.33'),Q('.12'),F(0)];gs=[[Q('.001'),F(0),F(0)],[F(0),Q('.001'),F(0)],[F(0),F(0),F(0)]]
 r=[F(0)]*3;B=[[F(i==j) for j in range(3)] for i in range(3)]
 cb=c[:];gb=[g[:] for g in gs];rb=r[:];box=hull(c,gs,r);rows=[]
 for n in range(201):
  if n in (0,20,40,80,120,160,200):
   X=hull_frame(c,gs,B,r);Y=hull(cb,gb,rb);t=n*h
   rows.append(dict(time=float(t),direct_width=float(sum(x[1]-x[0] for x in box[:2])),fixed_affine_width=float(sum(x[1]-x[0] for x in Y[:2])),moving_frame_width=float(sum(x[1]-x[0] for x in X[:2])),fixed_J_upper=float(Y[2][1]),moving_J_upper=float(X[2][1]),moving_energy_lower=float(Q('.1')+Q('.005')*t-X[2][1]),floor_pass=Q('.1')+Q('.005')*t-X[2][1]>=Q('.08')))
  if n<200:
   box,_=step_box(box);cb,gb,rb,_=step_aff(cb,gb,rb);c,gs,B,r=step_frame(c,gs,B,r)
 with open('/mnt/data/adaptive_lohner_frame_results.csv','w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
 for row in rows:print(row)
