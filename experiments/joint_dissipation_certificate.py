#!/usr/bin/env python3
"""Exact-rational invariant-region and joint disagreement/expenditure certificate."""
from fractions import Fraction as F
import csv
from pathlib import Path
Q=lambda s:F(str(s))
k=Q('.7'); ca=Q('.1'); amin=Q('.15'); amax=Q('.19')
p0lo=Q('.329'); pstar=Q('.325'); mmax=Q('.121'); m0lo=Q('.119')
# At p=pstar the vector field points into p>=pstar for all a,m.
p_boundary_margin=amin-pstar**2-3*mmax**2
# At m=mmax the vector field points toward smaller m, for all p>=pstar.
Dmin=2*k-amax+3*pstar**2
assert p0lo>=pstar and p_boundary_margin>0 and Dmin>0
C=8*ca*k*k
Jbound=C*mmax*mmax/(2*Dmin)
old=Q('0.002081170556266443')
assert Jbound<old
# J(t)+(C/(2Dmin))m(t)^2 <= (C/(2Dmin))m(0)^2.
rows=[]
for horizon in [Q('.4'),Q('.8'),Q('1.2'),Q('1.6'),Q('2'),Q('10')]:
 rows.append({'horizon':str(horizon),'joint_J_upper':str(Jbound),'joint_J_upper_decimal':f'{float(Jbound):.12f}', 'uniform_energy_lower_decimal':f'{float(Q(".1")-Jbound):.12f}', 'floor_pass':str(Q('.1')-Jbound>Q('.08'))})
path=Path('/mnt/data/joint_dissipation_results.csv')
with path.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print('p boundary inward margin:',p_boundary_margin,'=',float(p_boundary_margin))
print('dissipation D lower:',Dmin,'=',float(Dmin))
print('J upper:',Jbound,'=',float(Jbound))
print('improvement over fixed affine at H=2:',float((old-Jbound)/old)*100,'%')
print('uniform energy lower:',float(Q('.1')-Jbound))
