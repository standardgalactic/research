from fractions import Fraction as F
import csv
A=(F('0.15'),F('0.19')); k=F('0.7'); ca=F('0.1'); h=F('0.01'); q=F('0.01'); basal=F('0.005'); E0=F('0.1'); floor=F('0.08')
def plus(x,y):return x[0]+y[0],x[1]+y[1]
def minus(x,y):return x[0]-y[1],x[1]-y[0]
def times(x,y):
 p=[a*b for a in x for b in y];return min(p),max(p)
def sc(x,a):return times(x,(a,a))
def sq(x):return (F(0) if x[0]<=0<=x[1] else min(x[0]**2,x[1]**2),max(x[0]**2,x[1]**2))
def cube(x):return x[0]**3,x[1]**3
def hull(x,y):return min(x[0],y[0]),max(x[1],y[1])
def incl(x,y):return x[0]<=y[0] and y[1]<=x[1]
def mag(x):return max(abs(x[0]),abs(x[1]))
def vector(X):
 u,v,J,E=X; dif=minus(u,v)
 return [minus(minus(times(A,u),cube(u)),sc(dif,k)),minus(minus(times(A,v),cube(v)),sc(minus(v,u),k)),sc(sq(dif),2*ca*k*k),minus((q,q),plus((basal,basal),sc(sq(dif),2*ca*k*k)))]
def dyad(x,bits):
 den=1<<bits
 return (F((x[0]*den).__floor__(),den),F((x[1]*den).__ceil__(),den))
def step(X,bits):
 # Picard tube Y contains X + [0,h] f(Y). Fixed point inclusion implies existence/uniqueness
 # on this step; smooth polynomial RHS ensures uniqueness.
 fx=vector(X); radius=[h*(mag(v)+F('0.01')) for v in fx]
 for it in range(50):
  Y=[(x[0]-r,x[1]+r) for x,r in zip(X,radius)]
  fy=vector(Y)
  if all(incl(Yi,hull(x,plus(x,sc(fi,h)))) for Yi,x,fi in zip(Y,X,fy)):break
  radius=[2*r for r in radius]
 else:raise ValueError('no Picard tube')
 # integrate derivative range over verified tube, round endpoints outward
 Z=[dyad(plus(x,sc(fi,h)),bits) for x,fi in zip(X,fy)]
 return Z,it+1,Y[3][0]

def run(bits,H):
 X=[(F('0.45'),F('0.45')),(F('0.21'),F('0.21')),(F(0),F(0)),(E0,E0)]
 maxiters=0; maxwidth=F(0); safe=True
 for i in range(int(H/h)):
  X,it,tube_energy_lower=step(X,bits);maxiters=max(it,maxiters)
  maxwidth=max(maxwidth,X[0][1]-X[0][0],X[1][1]-X[1][0]);safe=safe and X[3][0]>=floor and tube_energy_lower>=floor
 return [bits,float(H),float(X[0][1]-X[0][0]+X[1][1]-X[1][0]),float(X[2][1]),float(X[3][0]),float(maxwidth),maxiters,safe]
rows=[]
for bits in (24,40,64):
 for H in (F('0.1'),F('0.2'),F('0.4')):
  row=run(bits,H);rows.append(row);print(row,flush=True)
with open('/mnt/data/dyadic_viability_results.csv','w',newline='') as f:
 w=csv.writer(f);w.writerow(['dyadic_bits','horizon','state_width_sum','actuation_upper','energy_lower','max_state_width','max_picard_iterations','resource_floor_certified']);w.writerows(rows)
