import csv
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors
p='/mnt/data/long_horizon_affine_report.pdf'
rows=list(csv.DictReader(open('/mnt/data/long_horizon_affine_results.csv')))
s=getSampleStyleSheet(); story=[]
def para(t,style='BodyText'):story.append(Paragraph(t,s[style]));story.append(Spacer(1,10))
para('Long-Horizon Continuous-Time Correlated-Generator Validation','Title')
para('Reduced real Stuart–Landau system: p\' = ap - p³ - 3pm², m\' = (a-2k)m - 3p²m - m³, J\' = 8 c_A k²m². Here a is persistent within [0.15,0.19], k = 0.7, c_A = 0.1, and initial intervals are p in [0.329,0.331], m in [0.119,0.121], J = 0.')
para('The experiment uses step h = 0.01, exact rational interval operations, outward 40-bit dyadic rounding, and a local Picard inclusion tube. A second-derivative interval bound encloses the continuous-time flow defect relative to a first-order step. Both direct intervals and a fixed correlated-generator representation are propagated from identical initial sets.')
data=[['H','Direct width','Affine width','Reduction','Affine J upper']]
for r in rows:
 t=float(r['time']); d=float(r['box_state_width']); a=float(r['affine_state_width']);j=float(r['affine_J_upper']);data.append([f'{t:.1f}',f'{d:.6f}',f'{a:.6f}',f'{100*(1-a/d):.1f}%',f'{j:.8f}'])
t=Table(data,repeatRows=1,colWidths=[43,93,93,75,104]);t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e6e9ed')),('GRID',(0,0),(-1,-1),.3,colors.grey),('FONTSIZE',(0,0),(-1,-1),9),('BOTTOMPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),7)]));story.append(t);story.append(Spacer(1,15))
para('Resource certificate. Let E(t) = 0.1 + (0.01 - 0.005)t - J(t), with floor E* = 0.08. Since J is nondecreasing and the net supply is nonnegative, for every t in [0,2], E(t) >= 0.1 - J(2) >= 0.0979188294 > 0.08. At the endpoint, E(2) >= 0.1079188294. These are distinct statements: the endpoint lower bound must not be presented as a uniform-in-time lower bound.')
para('Scope and interpretation. This is a longer-horizon stress test of the previously implemented fixed correlated-generator enclosure, not a matrix-valued moving-frame Lohner algorithm. The reported widths are sums of interval widths for p and m. Improved enclosure tightness does not establish optimality, and the numerical certificate applies only to this reduced mathematical model and its specified parameter range.')
para('Reproducibility. Run long_horizon_affine_validation.py using Python 3; it writes long_horizon_affine_results.csv. The experiment is deterministic and uses fractions.Fraction for interval endpoints.')
SimpleDocTemplate(p,pagesize=(595,842),leftMargin=48,rightMargin=48,topMargin=48,bottomMargin=48).build(story)
print(p)
