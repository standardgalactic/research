"""Dual continuous-time enclosures; rigorous intersection of expenditure projections."""
import runpy, csv
from fractions import Fraction as F
base=runpy.run_path('/mnt/data/long_horizon_affine_validation.py')
frame=runpy.run_path('/mnt/data/adaptive_lohner_frame.py')
Q=base['Q']; h=base['h']; I=base['I']; tube=base['tube']; rhs=base['rhs']; a=base['a']; sc=base['sc']; sq=base['sq']; ca=base['ca']; k=base['k']; dyad=base['dyad']
c=[Q('.33'),Q('.12'),F(0)];gs=[[Q('.001'),F(0),F(0)],[F(0),Q('.001'),F(0)],[F(0),F(0),F(0)]]
cf=c[:];gf=[g[:] for g in gs];rf=[F(0)]*3
B=[[F(i==j) for j in range(3)] for i in range(3)];r=[F(0)]*3
# A third, scalar-only bound integrates the nonnegative expenditure derivative using
# a validated trajectory tube from either state enclosure at each step.
Jquad=F(0); rows=[]
for n in range(201):
 X=frame['hull_frame'](c,gs,B,r);Y=base['hull'](cf,gf,rf)
 if n in (0,20,40,80,120,160,200):
  t=n*h; best=min(X[2][1],Y[2][1],Jquad)
  rows.append(dict(time=float(t),fixed_state_width=float(sum(z[1]-z[0] for z in Y[:2])),moving_state_width=float(sum(z[1]-z[0] for z in X[:2])),fixed_J_upper=float(Y[2][1]),moving_J_upper=float(X[2][1]),quadrature_J_upper=float(Jquad),dual_J_upper=float(best),dual_energy_lower=float(Q('.1')-best),resource_floor_pass=Q('.1')-best>=Q('.08')))
 if n==200:break
 # Both tubes enclose the same continuous-time state trajectories. Their
 # intersection is consequently a valid enclosure for the disagreement m.
 T1=tube(X);T2=tube(Y)
 mlo=max(T1[1][0],T2[1][0]);mhi=min(T1[1][1],T2[1][1]);assert mlo<=mhi
 dJ=sc(sq(I(mlo,mhi)),8*ca*k*k)
 Jquad=dyad(I(Jquad+h*dJ[1]))[1]
 cf,gf,rf,_=base['step_aff'](cf,gf,rf)
 c,gs,B,r=frame['step_frame'](c,gs,B,r)
with open('/mnt/data/observable_aware_dual_results.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
for row in rows:print(row)
