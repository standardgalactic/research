#!/usr/bin/env python3
"""Exact linear-phase and energy trajectories; illustrative model only."""
import csv
import numpy as np
from scipy.linalg import expm
N=8
L=np.zeros((N,N))
for i in range(N):
    for j in ((i-1)%N,(i+1)%N):
        L[i,i]+=1; L[i,j]-=1
P=np.eye(N)-np.ones((N,N))/N
w0=P@np.linspace(0,1.4,N)
E0,Estar,q,d,c,H=1.,.1,.2,.1,.02,40.
rows=[]
for k in (.5,1.,2.,3.,4.):
    w=expm(-k*L*H)@w0
    slope=q-d-c*k*k
    exit_time=(E0-Estar)/(-slope) if slope<0 else float('inf')
    rows.append(dict(k=k,energy_at_H=round(E0+H*slope,6),viable_through_H=(E0+min(0,H*slope)>=Estar),exit_time=round(exit_time,6) if np.isfinite(exit_time) else 'none',transverse_norm_at_H=round(np.linalg.norm(w),9)))
with open('/mnt/data/coupling_viability_results.csv','w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
for row in rows:print(row)
