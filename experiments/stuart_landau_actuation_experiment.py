#!/usr/bin/env python3
"""Full complex-amplitude Stuart--Landau network with a specified resistive actuation law.
All parameters are dimensionless illustrative values, not calibrated physiology.
"""
import csv
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent
N=8
L=np.zeros((N,N))
for i in range(N):
    j=(i+1)%N
    L[i,i]+=1;L[j,j]+=1;L[i,j]-=1;L[j,i]-=1
P=np.eye(N)-np.ones((N,N))/N
omega=np.array([-.16,-.10,-.05,-.015,.015,.05,.10,.16])
phase=np.linspace(0,2*np.pi,N,endpoint=False)
z0=(.90+.06*np.cos(phase))*np.exp(1j*phase)
E0=2.0; Emin=.35; H=30.0; tau=8.0
basal=.028; supply=.028; c_act=.045; c_osc=.012
# An explicit electrical analogue: actuator current I_i=u_i/chi,
# resistive power R sum |I_i|^2 = (R/chi^2) sum |u_i|^2.
# Here c_act represents R/chi^2, not an inferred biological coefficient.

def growth(E):
    return 0.8+0.12*np.tanh((E-1.0)/.8)

def rhs(k, hetero=True):
    w=omega if hetero else np.zeros(N)
    def f(t,y):
        z=y[:N]+1j*y[N:2*N]
        E=y[-1]
        u=-k*(L@z)
        dz=(growth(E)+1j*w-np.abs(z)**2)*z+u
        power=c_act*np.vdot(u,u).real+c_osc*np.vdot(z,z).real/N
        dE=supply-basal-power
        return np.r_[dz.real,dz.imag,dE]
    return f

def run(k,hetero=True):
    y0=np.r_[z0.real,z0.imag,E0]
    f=rhs(k,hetero)
    def exit_event(t,y):return y[-1]-Emin
    exit_event.terminal=True;exit_event.direction=-1
    sol=solve_ivp(f,(0,H),y0,rtol=1e-8,atol=1e-10,events=exit_event,dense_output=True,max_step=.08)
    t_end=float(sol.t[-1]); exit_time=float(sol.t_events[0][0]) if len(sol.t_events[0]) else None
    # Only report deadline coherence if the system remains operational to tau.
    if t_end>=tau-1e-8:
        y=sol.sol(tau); z=y[:N]+1j*y[N:2*N]
        disagreement=float(np.linalg.norm(P@z))
    else: disagreement=float('nan')
    z=sol.y[:N,-1]+1j*sol.y[N:2*N,-1]
    return dict(k=k,heterogeneous=hetero,exit_time=exit_time if exit_time is not None else '',
                survived_horizon=exit_time is None,energy_at_stop=float(sol.y[-1,-1]),
                disagreement_at_deadline=disagreement,
                disagreement_at_stop=float(np.linalg.norm(P@z)),
                stop_time=t_end)

if __name__=='__main__':
    rows=[]
    for scenario,e0,c in (('ample-reserve',2.0,.045),('scarce-reserve',.9,.12)):
        E0=e0; c_act=c
        for het in (False,True):
            for k in (.05,.15,.3,.6,1.,2.,3.,4.):
                rows.append(dict(scenario=scenario,**run(k,het)))
    with (HERE/'stuart_landau_actuation_results.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=rows[0].keys());wr.writeheader();wr.writerows(rows)
    for r in rows:print(r)
