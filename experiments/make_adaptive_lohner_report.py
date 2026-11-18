from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
import csv
out='/mnt/data/adaptive_lohner_frame_report.pdf';styles=getSampleStyleSheet();styles.add(ParagraphStyle(name='SmallText2',parent=styles['BodyText'],fontSize=9.4,leading=14,spaceAfter=10))
D=SimpleDocTemplate(out,pagesize=(595,842),leftMargin=52,rightMargin=52,topMargin=50,bottomMargin=48)
story=[]
def p(s):story.append(Paragraph(s,styles['SmallText2']))
story.append(Paragraph('Matrix-Transported Remainders in a Validated Stuart–Landau Enclosure',styles['Title']))
story.append(Paragraph('Controlled comparison at identical initial uncertainty and arithmetic precision',styles['Italic']));story.append(Spacer(1,15))
p('<b>Model.</b> The real two-oscillator system is expressed in collective and disagreement coordinates: p\u2032 = ap - p³ - 3pm²; m\u2032 = (a - 2k)m - 3p²m - m³; J\u2032 = 8c<sub>A</sub>k²m². Here a belongs to [0.15, 0.19] and remains constant along each trajectory, k = 0.7 and c<sub>A</sub> = 0.1. Initial center is (p,m,J) = (0.33,0.12,0), with two independent initial state perturbations of radius 0.001. The parameter perturbation has radius 0.02.')
p('<b>Enclosure representation.</b> The fixed method uses x = c + Σ g<sub>i</sub> ξ<sub>i</sub> + [-r,r], with |ξ<sub>i</sub>| ≤ 1. The moving method uses x = c + Σ g<sub>i</sub> ξ<sub>i</sub> + B[-r,r], with an invertible rational 3×3 matrix B. Both use the same center and generator evolution. In the moving method, B<sub>n+1</sub> = (I + h Df(c<sub>n</sub>,a<sub>0</sub>))B<sub>n</sub>. The interval remainder is propagated in the new coordinates using the exact rational inverse of B<sub>n+1</sub>. No inverse is approximated numerically.')
p('<b>Continuous-time validation.</b> For each step, a strict Picard inclusion X + h f(T,A) contained in the interior of T provides a flow tube T. The exact-flow defect relative to Euler is bounded componentwise by (h²/2) sup<sub>T,A</sub>|D<sub>x</sub>f · f|. The uncertain parameter is constant in time, so no a\u2032 term appears. The Jacobian variation, parameter-affine residual, and defect are included as symmetric physical-coordinate forcing, which is mapped into the moving remainder frame using |B<sub>n+1</sub><sup>-1</sup>|. Dyadic rounding of centers, generators, frame entries and remainder radii is outward compensated. The construction is a conservative first-order validated moving-frame enclosure, not a high-order or optimized QR-based Lohner integrator.')
rows=list(csv.DictReader(open('/mnt/data/adaptive_lohner_frame_results.csv')))
data=[['H','Direct width','Fixed width','Moving width','Fixed J upper','Moving J upper']]
for r in rows:
 if float(r['time'])>0:
  data.append([r['time'],*[f'{float(r[k]):.6f}' for k in ['direct_width','fixed_affine_width','moving_frame_width','fixed_J_upper','moving_J_upper']]])
t=Table(data,colWidths=[40,86,82,88,92,94],repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e7edf2')),('GRID',(0,0),(-1,-1),.35,colors.grey),('FONTSIZE',(0,0),(-1,-1),7.5),('BOTTOMPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),7)]));story.append(t);story.append(Spacer(1,14))
p('<b>Interpretation.</b> At H = 2, matrix transport reduces combined state width from 0.0636043 to 0.0285487, a reduction of 55.1%. The expenditure upper bound increases from 0.00208117 to 0.00231344, however. This is a mixed outcome: a state-space wrapping improvement is not automatically an improvement in the accumulated nonlinear observable J. Direct interval propagation gives width 0.193859 at the same horizon.')
p('<b>Resource floor.</b> With E(0)=0.1, q=0.01, b=0.005 and E*=0.08, the nondecreasing expenditure J implies E(t) ≥ 0.1 - J(2) ≥ 0.09768656 for 0 ≤ t ≤ 2. At the endpoint, E(2) ≥ 0.10768656. These bounds omit any other energy demands beyond those explicitly represented in this reduced model.')
p('<b>Limitations.</b> The experiment uses a first-order center and generator transport with conservative residual forcing. It does not prove that matrix transport is optimal, nor does it establish a general improvement theorem. In particular, expenditure enclosures are sensitive to wrapping and accumulated defect in the J coordinate. Independent validation against higher-order Taylor-model flow integration would be a valuable next check. The script and CSV reproduce the reported calculations using Python Fraction arithmetic.')
D.build(story)
print(out)
