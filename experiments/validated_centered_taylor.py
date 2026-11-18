"""Outward-rounded first-order parameter Taylor enclosures, short-horizon prototype.
All interval endpoints and parameter cells use decimal strings; mpmath.iv outward rounding.
"""
import csv
from mpmath import iv
from scipy.integrate import solve_ivp
import numpy as np
iv.dps=40
N=3; L=((1,-1,0),(-1,2,-1),(0,-1,1)); nu=(-1,0,1)
k=iv.mpf('0.7'); s=iv.mpf('0.16'); ca=iv.mpf('0.1')
z0=['0.45','0.21','0.33','0.10','-0.06','0.04']
B=iv.mpf('1'); F=iv.mpf('3.48'); Lam=iv.mpf('5.48'); h=iv.mpf('0.002'); HMAX=iv.mpf('0.2')
MJ=2*ca*k*k*9*B*F

def square(x):
    # intersection of generic interval product with nonnegative half-line
    return x*x

def f(z,a):
    out=[]
    for i in range(3):
        r=z[i]; y=z[i+3]; mod=square(r)+square(y)
        out.append(a*r-s*nu[i]*y-mod*r-k*sum(L[i][j]*z[j] for j in range(3)))
    for i in range(3):
        r=z[i]; y=z[i+3]; mod=square(r)+square(y)
        out.append(a*y+s*nu[i]*r-mod*y-k*sum(L[i][j]*z[j+3] for j in range(3)))
    return out

def sens(z,S,a):
    out=[]
    for i in range(3):
        r=z[i]; y=z[i+3]; u=S[i]; v=S[i+3]
        mod=square(r)+square(y); cross=2*(r*u+y*v)
        out.append(a*u-s*nu[i]*v-mod*u-cross*r-k*sum(L[i][j]*S[j] for j in range(3))+r)
    for i in range(3):
        r=z[i]; y=z[i+3]; u=S[i]; v=S[i+3]
        mod=square(r)+square(y); cross=2*(r*u+y*v)
        out.append(a*v+s*nu[i]*u-mod*v-cross*y-k*sum(L[i][j]*S[j+3] for j in range(3))+y)
    return out

def lap(z):return [sum(L[i][j]*z[j+o] for j in range(3)) for o in (0,3) for i in range(3)]
def power(z):return ca*k*k*sum(square(v) for v in lap(z))
def p_sens(z,S):return 2*ca*k*k*sum(x*y for x,y in zip(lap(z),lap(S)))
def pm(x):return iv.mpf([-1,1])*x

def boundU(t):return B*(iv.exp(Lam*t)-1)/Lam
def boundV(t):
    U=boundU(t)
    return (6*B*U*U+2*U)*(iv.exp(Lam*t)-1)/Lam
def boundW(t):
    U=boundU(t); V=boundV(t)
    return 2*ca*k*k*9*t*(U*U+B*V)

def center(a,T):
    z=[iv.mpf(v) for v in z0]; S=[iv.mpf('0') for _ in range(6)]; J=iv.mpf('0');K=iv.mpf('0')
    steps=round(T/.002)
    for j in range(steps):
        t=iv.mpf(j)*h
        U=boundU(t+h)
        MS=Lam*(Lam*U+B)+6*B*F*U+F
        MK=2*ca*k*k*9*(F*U+B*(Lam*U+B))
        dz=f(z,a);dS=sens(z,S,a);dJ=power(z);dK=p_sens(z,S)
        z=[v+h*dv+pm(h*h*Lam*F/2) for v,dv in zip(z,dz)]
        S=[v+h*dv+pm(h*h*MS/2) for v,dv in zip(S,dS)]
        J=J+h*dJ+pm(h*h*MJ/2)
        K=K+h*dK+pm(h*h*MK/2)
    return z,S,J,K

def endpoints(x):return float(x.a),float(x.b)
def width(x):
    l,u=endpoints(x);return u-l

def centered(a0,a1,T):
    mid=(iv.mpf(a0)+iv.mpf(a1))/2; d=(iv.mpf(a1)-iv.mpf(a0))/2
    z,S,J,K=center(mid,T)
    Rz=boundV(iv.mpf(str(T)))*d*d/2
    RJ=boundW(iv.mpf(str(T)))*d*d/2
    zz=[v+pm(d)*ss+pm(Rz) for v,ss in zip(z,S)]
    JJ=J+pm(d)*K+pm(RJ)
    return zz,JJ

def direct(a0,a1,T):
    z=[iv.mpf(v) for v in z0];J=iv.mpf('0');a=iv.mpf([a0,a1])
    for _ in range(round(T/.002)):
        dz=f(z,a);p=power(z)
        z=[v+h*dv+pm(h*h*Lam*F/2) for v,dv in zip(z,dz)]
        J=J+h*p+pm(h*h*MJ/2)
    return z,J

def true_state(a,T):
    Lnp=np.array(L,float); zinit=np.array([complex(float(z0[i]),float(z0[i+3])) for i in range(3)])
    def rhs(t,y):
        z=y[:3]+1j*y[3:6];dz=(a+1j*.16*np.array(nu))*z-np.abs(z)**2*z-.7*Lnp@z
        p=.1*.7**2*np.linalg.norm(Lnp@z)**2
        return np.r_[dz.real,dz.imag,p]
    y=solve_ivp(rhs,[0,T],np.r_[zinit.real,zinit.imag,0.],rtol=1e-12,atol=1e-14).y[:,-1]
    return y

rows=[]
for n in (8,32):
    for T in (.02,.05,.1,.2):
        maxdirect=maxcenter=0.; maxjdirect=maxjcenter=0.; inside=True
        for i in range(n):
            from decimal import Decimal
            a0=str(Decimal('0.12')+Decimal('0.10')*Decimal(i)/Decimal(n))
            a1=str(Decimal('0.12')+Decimal('0.10')*Decimal(i+1)/Decimal(n))
            zd,jd=direct(a0,a1,T);zc,jc=centered(a0,a1,T)
            maxdirect=max(maxdirect,max(map(width,zd)))
            maxcenter=max(maxcenter,max(map(width,zc)))
            maxjdirect=max(maxjdirect,endpoints(jd)[1]);maxjcenter=max(maxjcenter,endpoints(jc)[1])
            # Numerical sanity check only, not a proof of interval validity.
            for a in (float(a0), (float(a0)+float(a1))/2,float(a1)):
                yt=true_state(a,T)
                if not all(endpoints(v)[0]-1e-10<=val<=endpoints(v)[1]+1e-10 for v,val in zip(zc,yt[:6])) or not (endpoints(jc)[0]-1e-10<=yt[6]<=endpoints(jc)[1]+1e-10):inside=False
        row=[n,T,maxdirect,maxcenter,maxjdirect,maxjcenter,inside]
        rows.append(row);print(row,flush=True)
with open('/mnt/data/validated_centered_taylor_results.csv','w',newline='') as out:
    w=csv.writer(out);w.writerow(['cells','horizon','direct_max_coordinate_width','centered_max_coordinate_width','direct_max_J_upper','centered_max_J_upper','sampled_trajectories_contained']);w.writerows(rows)
