"""Exact-rational moving affine frame benchmark for the tangent Stuart-Landau map.
This is a frozen-Jacobian discrete surrogate, NOT validated nonlinear flow.
"""
from fractions import Fraction as F
import csv
N=6
L=[[1,-1,0],[-1,2,-1],[0,-1,1]]
r=[F('0.45'),F('0.21'),F('0.33')]
im=[F('0.10'),F('-0.06'),F('0.04')]
a=F('0.17'); k=F('0.7'); s=F('0.16'); nu=[-1,0,1]
h=F('0.05')
# real state order r1,r2,r3,i1,i2,i3
J=[[F(0) for j in range(N)] for i in range(N)]
for i in range(3):
 q=r[i]**2+im[i]**2
 for j in range(3):
  d=int(i==j)
  J[i][j]=F(d)*(a-q-2*r[i]**2)-k*L[i][j]
  J[i+3][j+3]=F(d)*(a-q-2*im[i]**2)-k*L[i][j]
  J[i][j+3]=F(d)*(-2*r[i]*im[i]-s*nu[i])
  J[i+3][j]=F(d)*(-2*r[i]*im[i]+s*nu[i])
M=[[F(int(i==j))+h*J[i][j] for j in range(N)] for i in range(N)]
def mm(A,B):return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def mv(A,x):return [sum(A[i][j]*x[j] for j in range(len(x))) for i in range(len(A))]
# parameter uncertainty is a shared scalar; initial direction from growth sensitivity
v=[F('0.12'),F('-0.07'),F('0.04'),F('0.02'),F('0.06'),F('-0.03')]
d=F('0.01'); eps=F('0.00001')
# moving affine frame: exact M^n applied to initial independent remainder symbols and parameter symbol
A=[[F(int(i==j))*eps for j in range(N)]+[d*v[i]] for i in range(N)]
# axis-aligned interval propagation, initially hull(A)
b=[sum(abs(t) for t in row) for row in A]
rows=[]
for step in range(81):
 af=[sum(abs(x) for x in row) for row in A]
 box=sum(b); affine=sum(af)
 rows.append([step,float(h*step),float(box),float(affine),float(box/affine),str(box),str(affine)])
 if step<80:
  A=mm(M,A)
  b=mv([[abs(t) for t in row] for row in M],b)
with open('/mnt/data/moving_frame_results.csv','w',newline='') as f:
 w=csv.writer(f);w.writerow(['steps','time','box_total_halfwidth','affine_total_halfwidth','wrapping_factor','box_exact','affine_exact']);w.writerows(rows)
for i in [0,10,20,40,80]:print(rows[i][:5])
