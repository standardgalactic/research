"""Single reproduction procedure for the Chapter 8 computations.
Usage: python3 reproduce.py          (run from the directory holding the five programs and reference CSVs)
Re-runs each program in a scratch directory, compares every regenerated table to the committed reference table
(tolerance 1e-9 relative, exact for integer columns), then checks the headline figures quoted in the monograph.
Dependencies are pinned in requirements.txt.  Random seeds are fixed inside the programs
(calibration 20261118, mode-resolved 7, amplitude envelope 11).  Exit status 0 means every check passed."""
import csv, os, shutil, subprocess, sys, tempfile, math
HERE=os.path.dirname(os.path.abspath(__file__))
PROGS=[('young_parameter_optimization.py',['young_parameter_results.csv','young_finite_horizon_results.csv']),
       ('certificate_calibration.py',['certificate_calibration_results.csv']),
       ('mode_resolved_certificate.py',['mode_resolved_results.csv']),
       ('amplitude_envelope_certificate.py',['amplitude_envelope_results.csv'])]
fails=[]
def check(ok,msg):
    print(('PASS ' if ok else 'FAIL ')+msg)
    if not ok: fails.append(msg)
def read(p): return list(csv.DictReader(open(p)))
def same(a,b):
    if a==b: return True
    try: x,y=float(a),float(b)
    except ValueError: return False
    return abs(x-y)<=1e-9*max(1.0,abs(x),abs(y))
tmp=tempfile.mkdtemp()
for f in os.listdir(HERE):
    if f.endswith('.py') and f!='reproduce.py': shutil.copy(os.path.join(HERE,f),tmp)
for prog,outs in PROGS:
    r=subprocess.run([sys.executable,prog],cwd=tmp,capture_output=True,text=True)
    check(r.returncode==0,f'{prog} runs')
    for o in outs:
        new=read(os.path.join(tmp,o)); ref=read(os.path.join(HERE,o))
        ok=len(new)==len(ref) and all(same(a[k],b[k]) for a,b in zip(new,ref) for k in b)
        check(ok,f'{o} matches reference ({len(new)} rows)')
def col(name,c): return [float(r[c]) for r in read(os.path.join(HERE,name))]
cal=read(os.path.join(HERE,'certificate_calibration_results.csv'))
check(len(cal)==216,'calibration has 216 rows')
check(all(float(r['bound_theta_half'])>=float(r['J_actual']) and float(r['bound_theta_optimal'])>=float(r['J_actual']) and float(r['bound_disagreement_route'])>=float(r['J_actual']) for r in cal),'every calibration bound >= integrated expenditure')
mr=read(os.path.join(HERE,'mode_resolved_results.csv'))
check(len(mr)==144 and all(r['valid']=='1' for r in mr),'144 modal comparisons, all valid')
check(abs(max(float(r['gain_pct']) for r in mr)-61.6)<0.05,'maximum modal gain 61.6 percent')
ae=read(os.path.join(HERE,'amplitude_envelope_results.csv'))
check(len(ae)==120 and all(r['valid']=='1' for r in ae),'120 envelope comparisons, all valid')
check(abs(max(float(r['gain_total_pct']) for r in ae)-99.4)<0.05,'maximum combined gain 99.4 percent')
import statistics
for state,T,exp in (('amp_0.45',1.0,0.12),('amp_0.05',5.0,61.24),('amp_0.01',5.0,96.19),('amp_0.20',100.0,1.99)):
    v=[float(r['gain_young_pct']) for r in ae if r['state']==state and float(r['T'])==T]
    check(abs(statistics.median(v)-exp)<0.006,f'envelope median gain {state}, T={T}: {exp}')
yp=read(os.path.join(HERE,'young_parameter_results.csv'))
p=[r for r in yp if r['graph']=='path_4'][0]
check(abs(float(p['buffer_theta_half'])-0.015792)<1e-6 and abs(float(p['buffer_optimized'])-0.007972)<1e-6,'path buffers 0.015792 and 0.007972')
shutil.rmtree(tmp)
print('ALL CHECKS PASSED' if not fails else f'{len(fails)} CHECK(S) FAILED'); sys.exit(1 if fails else 0)
