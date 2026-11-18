"""Young-parameter allocation for the direct spectral expenditure certificate.
Bound: J(T) <= cA k^2 W0/((1-th) beta) + cA k^2 F^2 T/(4 th (1-th) beta^2), 0<th<1.
Exact rational arithmetic except square roots, which are rounded upward at 2^-80."""
from fractions import Fraction as Fr
from math import isqrt
import csv, sys

def sqrt_up(x, bits=80):
    s = 1 << bits
    n = (x.numerator * s * s + x.denominator - 1) // x.denominator
    r = isqrt(n)
    if r * r < n: r += 1
    return Fr(r, s)

a, k, cA, R = Fr(19,100), Fr(7,10), Fr(1,10), Fr(45,100)
om = [Fr(-3,100), Fr(-1,100), Fr(1,100), Fr(3,100)]
z0 = [Fr(45,100), Fr(21,100), Fr(39,100), Fr(27,100)]
E0, s = Fr(1,10), Fr(5,1000)
graphs = {
 'path_4': ([(0,1),(1,2),(2,3)], Fr(117,200)),
 'cycle_4': ([(0,1),(1,2),(2,3),(3,0)], Fr(2)),
 'complete_4': ([(i,j) for i in range(4) for j in range(i+1,4)], Fr(4)),
}
mean = sum(om)/4
omt = [w-mean for w in om]
delta = max(abs(w) for w in omt)
F = R*sqrt_up(sum(w*w for w in omt))

def bound(W0, beta, T, th):
    return cA*k*k*W0/((1-th)*beta) + cA*k*k*F*F*T/(4*th*(1-th)*beta*beta)

rows = []
for name,(edges,l2) in graphs.items():
    W0 = sum((z0[i]-z0[j])**2 for i,j in edges)/2
    beta = k-(a+delta)/l2
    A = cA*k*k*F*F/(beta*beta)          # persistent rate at theta=1/2
    B = cA*k*k*W0/beta                  # transient coefficient
    # theta=1/2 baseline
    buf_half = 2*B
    # supply-matched allocation: smallest theta with rate(theta) <= s
    if s >= A:
        th_s = (1-sqrt_up(1-A/s))/2
        # round theta up slightly so the rate condition remains satisfied
        th_s = th_s + Fr(1,10**20)
    else:
        th_s = None
    buf_opt = B/(1-th_s) if th_s else None
    # finite-horizon optimum over theta
    fin = {}
    for T in (1,5,20,100):
        kap = cA*k*k*F*F*T/(4*beta*beta)
        al = B
        val = (sqrt_up(al+kap)+sqrt_up(kap))**2
        th = sqrt_up(kap)/(sqrt_up(kap)+sqrt_up(al+kap))
        fin[T] = (float(val), float(th), float(bound(W0,beta,T,Fr(1,2))), float(bound(W0,beta,T,th)))
    rows.append(dict(graph=name, beta=float(beta), W0=float(W0), rate_half=float(A), buffer_half=float(buf_half),
        theta_supply=float(th_s) if th_s else None, buffer_opt=float(buf_opt) if buf_opt else None,
        floor_half=float(E0-buf_half), floor_opt=float(E0-buf_opt) if buf_opt else None,
        buffer_ratio=float(buf_opt/buf_half) if buf_opt else None, fin=fin))
with open('young_parameter_results.csv','w',newline='') as f:
    w = csv.writer(f)
    w.writerow(['graph','beta','W0','rate_theta_half','buffer_theta_half','theta_supply_matched','buffer_optimized','floor_theta_half','floor_optimized','buffer_ratio'])
    for r in rows:
        w.writerow([r['graph'],r['beta'],r['W0'],r['rate_half'],r['buffer_half'],r['theta_supply'],r['buffer_opt'],r['floor_half'],r['floor_opt'],r['buffer_ratio']])
with open('young_finite_horizon_results.csv','w',newline='') as f:
    w = csv.writer(f)
    w.writerow(['graph','T','J_bound_optimal','theta_star','J_bound_theta_half','J_bound_check_at_theta_star'])
    for r in rows:
        for T,(v,th,h,c) in r['fin'].items():
            assert abs(v-c) < 1e-12*max(1,v), (v,c)
            w.writerow([r['graph'],T,v,th,h,c])
for r in rows: print({k_:v for k_,v in r.items() if k_!='fin'}); print(r['fin'])
