#!/usr/bin/env python3
"""Reproducible illustrative Stuart-Landau resource experiments; no biological data."""
import csv
import numpy as np
from scipy.integrate import solve_ivp

N=8
L=np.zeros((N,N))
for i in range(N):
    for j in ((i-1)%N,(i+1)%N):
        L[i,i]+=1; L[i,j]-=1
omega=1+0.06*(np.arange(N)-(N-1)/2)
e_star=r_star=0.1
z0=np.exp(1j*np.linspace(0,1.4,N))
y0=np.r_[z0.real,z0.imag,1.,1.]

def run(k0,q):
    def rhs(t,y):
        z=y[:N]+1j*y[N:2*N]; e,r=y[-2:]
        k=k0*np.clip(e/e_star,0,1)
        dz=(1+1j*omega-np.abs(z)**2)*z-k*(L@z)
        P=np.mean(np.abs(z)**2); D=np.real(np.vdot(z,L@z))
        s=.08*r
        de=s-.12-.02*P-.005*D-.005*e
        dr=q-s-.005*r
        return np.r_[dz.real,dz.imag,de,dr]
    def event(t,y):return min(y[-2]-e_star,y[-1]-r_star)
    event.terminal=True;event.direction=-1
    sol=solve_ivp(rhs,(0,80),y0,events=event,rtol=1e-8,atol=1e-10,max_step=.1,dense_output=True)
    t=sol.t[-1]; y=sol.y[:,-1]; z=y[:N]+1j*y[N:2*N]
    return dict(k0=k0,supply=q,exit_time=round(t,5),exited=bool(sol.t_events[0].size),transverse_norm=round(float(np.linalg.norm(z-z.mean())),6),min_amplitude=round(float(np.min(np.abs(z))),6),energy=round(float(y[-2]),6),reservoir=round(float(y[-1]),6))

rows=[run(k,q) for q in (.03,.25) for k in (.2,1.,4.)]
with open('/mnt/data/synchronization_results.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
for row in rows:print(row)
