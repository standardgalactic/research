"""Exact-rational enclosures for a NONLINEAR discrete Euler Stuart--Landau map.
No claim of validated continuous-time flow is made.
"""
from fractions import Fraction as Q
import csv
h=Q('0.025'); a0=Q('0.17'); d=Q('0.02'); k=Q('0.7'); cA=Q('0.1')
N=8

def add(I,J):return (I[0]+J[0],I[1]+J[1])
def mul(I,J):
 p=[x*y for x in I for y in J];return(min(p),max(p))
def scale(I,v):return (min(v*I[0],v*I[1]),max(v*I[0],v*I[1]))
def hull(c,r):return(c-r,c+r)
def cube(I):return(I[0]**3,I[1]**3)
def sq(I):
 lo=Q(0) if I[0]<=0<=I[1] else min(I[0]**2,I[1]**2)
 return (lo,max(I[0]**2,I[1]**2))
def f(x,a):
 x1,x2,J=x
 return [a*x1-x1**3-k*(x1-x2),a*x2-x2**3-k*(x2-x1),2*cA*k*k*(x1-x2)**2]
def jac(x):
 u,v,_=x
 return [[a0-3*u*u-k,k,Q(0)],[k,a0-3*v*v-k,Q(0)],[4*cA*k*k*(u-v),-4*cA*k*k*(u-v),Q(0)]]
def fb(X,A):
 u,v,J=X
 return [add(add(mul(A,u),scale(cube(u),-1)),scale(add(u,scale(v,-1)),-k)),
         add(add(mul(A,v),scale(cube(v),-1)),scale(add(v,scale(u,-1)),-k)),
         scale(sq(add(u,scale(v,-1))),2*cA*k*k)]
c=[Q('0.45'),Q('0.21'),Q(0)]; A=[Q(0),Q(0),Q(0)]; r=[Q(0)]*3
X=[hull(x,Q(0)) for x in c]
rows=[]
for step in range(N+1):
 widths=[2*(abs(A[i])+r[i]) for i in range(3)]
 boxwidth=[t[1]-t[0] for t in X]
 rows.append([step,float(step*h),*[float(v) for v in widths],*[float(v) for v in boxwidth],float(c[2]+abs(A[2])+r[2]),float(X[2][1]),float(sum(widths[:2])),float(sum(boxwidth[:2]))])
 if step==N:break
 fc=f(c,a0); M=jac(c)
 # shared parameter symbol xi in [-1,1], with a=a0+d xi
 newA=[A[i]+h*(sum(M[i][j]*A[j] for j in range(3))+(d*c[i] if i<2 else Q(0))) for i in range(3)]
 newc=[c[i]+h*fc[i] for i in range(3)]
 U=[abs(A[i])+r[i] for i in range(3)]
 # exact polynomial Taylor remainder in state/parameter deviations
 rem=[d*U[i]+3*abs(c[i])*U[i]**2+U[i]**3 for i in range(2)]
 rem.append(2*cA*k*k*(U[0]+U[1])**2)
 newr=[r[i]+h*(sum(abs(M[i][j])*r[j] for j in range(3))+rem[i]) for i in range(3)]
 # naive interval Euler of same nonlinear map
 FX=fb(X,hull(a0,d));X=[add(X[i],scale(FX[i],h)) for i in range(3)]
 c,A,r=newc,newA,newr
with open('/mnt/data/nonlinear_moving_affine_results.csv','w',newline='') as out:
 w=csv.writer(out);w.writerow(['step','time','affine_width_x1','affine_width_x2','affine_width_J','box_width_x1','box_width_x2','box_width_J','affine_J_upper','box_J_upper','affine_state_total_width','box_state_total_width']);w.writerows(rows)
for i in [0,2,4,6,8]:print(rows[i])
