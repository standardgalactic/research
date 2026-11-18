import csv
from mpmath import iv
iv.dps=35
N=3; k=iv.mpf('0.7'); sigma=iv.mpf('0.16'); ca=iv.mpf('0.1')
L=((1,-1,0),(-1,2,-1),(0,-1,1)); nu=(-1,0,1)
z0r=['0.45','0.21','0.33'];z0i=['0.10','-0.06','0.04']
B=iv.mpf('1.0') # ||z0||<1; a<=.22 => invariant ball B=1
# Norm bound: ||f|| <= (a_max+omega_max+B^2+k*lambda_max) B
# lambda_max(path3)=3
F=iv.mpf('0.22')+iv.mpf('0.16')+B*B+iv.mpf('2.1')
Lam=iv.mpf('0.22')+iv.mpf('0.16')+3*B*B+iv.mpf('2.1')
M=F*Lam
h=iv.mpf('0.002'); half=iv.mpf('0.5')
def lo(x):return float(x.a)
def hi(x):return float(x.b)
def width(x):return hi(x)-lo(x)
def sq(x):return x*x # interval product is an enclosure, including sign-crossing intervals
def f(x,a):
 r=x[:3];im=x[3:6];out=[]
 for i in range(3):
  rr=a*r[i]-sigma*nu[i]*im[i]-(sq(r[i])+sq(im[i]))*r[i]-k*sum(L[i][j]*r[j] for j in range(3))
  ii=a*im[i]+sigma*nu[i]*r[i]-(sq(r[i])+sq(im[i]))*im[i]-k*sum(L[i][j]*im[j] for j in range(3))
  out.append(rr)
  out.append(ii)
 return out[::2]+out[1::2]
def run(a0,a1,T):
 a=iv.mpf([a0,a1]);x=[iv.mpf(v) for v in z0r+z0i];J=iv.mpf('0');steps=round(T/.002)
 # |J''| <= 2 ca k^2 ||L||^2 B F, with lambda_max=3
 MJ=2*ca*k*k*iv.mpf('9')*B*F
 for _ in range(steps):
  fx=f(x,a)
  lr=[sum(L[i][j]*x[j] for j in range(3)) for i in range(3)]
  li=[sum(L[i][j]*x[j+3] for j in range(3)) for i in range(3)]
  power=ca*k*k*sum(sq(v) for v in lr+li)
  J=J+h*power+iv.mpf([-1,1])*half*h*h*MJ
  x=[x[i]+h*fx[i]+iv.mpf([-1,1])*half*h*h*M for i in range(6)]
 return x,J
rows=[]
for n in [8,16,32]:
 for T in [.02,.05,.1]:
  survives=0;maxw=0;true_kept=False;maxj=0
  for i in range(n):
   aa=.12+.1*i/n;bb=.12+.1*(i+1)/n
   x,J=run(str(aa),str(bb),T)
   # true parameter .173, observations from nominal high precision numerical trajectory unavailable;
   # report enclosures only, no unsupported exclusion claims
   maxw=max(maxw,max(width(v) for v in x));maxj=max(maxj,hi(J))
   if aa<=.173<=bb:true_kept=True
  rows.append([n,T,maxw,maxj,true_kept]);print(rows[-1],flush=True)
with open('/mnt/data/validated_joint_results.csv','w',newline='') as fp:
 w=csv.writer(fp);w.writerow(['parameter_cells','horizon','max_state_coordinate_interval_width','max_actuator_expenditure_upper','true_parameter_in_partition']);w.writerows(rows)
