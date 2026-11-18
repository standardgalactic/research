import csv
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.units import inch
rows=list(csv.DictReader(open('/mnt/data/observable_aware_dual_results.csv')))
out='/mnt/data/observable_aware_dual_report.pdf'
doc=SimpleDocTemplate(out,pagesize=(612,792),rightMargin=54,leftMargin=54,topMargin=50,bottomMargin=50)
s=getSampleStyleSheet();s.add(ParagraphStyle(name='SmallX',parent=s['BodyText'],fontSize=9.5,leading=14,spaceAfter=10));s.add(ParagraphStyle(name='TitleX',parent=s['Title'],fontSize=16,leading=20,spaceAfter=15));story=[]
def p(t):story.append(Paragraph(t,s['SmallX']))
story.append(Paragraph('Observable-Aware Dual Enclosures for Nonlinear Oscillator Resource Bounds',s['TitleX']))
p('Technical supplement | 18 November 2026 | Reduced real Stuart–Landau model')
p('<b>Question.</b> Can combining a fixed correlated-generator enclosure with a matrix-transported remainder enclosure improve both physical-state widths and cumulative actuator-expenditure bounds? We compare the two already validated propagators under identical parameters, and add a third scalar expenditure certificate from their intersected continuous-time tubes.')
p('<b>Dynamics and uncertainty.</b> The collective coordinate p=(u+v)/2 and disagreement m=(u-v)/2 obey p\' = ap - p³ - 3pm², m\' = (a-2k)m - 3p²m - m³, and J\' = 8 c_A k² m². The persistent growth parameter is a in [0.15,0.19], k=0.7, c_A=0.1. Initial values are p=0.33±0.001, m=0.12±0.001, J=0. The time step is 0.01, and dyadic outward rounding uses 40 bits.')
p('<b>Third certificate.</b> For each step [t_n,t_n+h], let T_fixed and T_moving be Picard-validated tubes. Since both contain every admissible solution, the disagreement projection M_n = projection_m(T_fixed) intersect projection_m(T_moving) also contains every exact m(t) on the step. With C=8 c_A k², the scalar recurrence U_(n+1) = round_up(U_n + h C sup{m²:m in M_n}) gives J(t_n) ≤ U_n. Both other methods independently enclose J(t_n), so taking the minimum of their upper endpoints and U_n is valid. This construction uses the same fixed parameter in each trajectory, although interval evaluations overapproximate its influence.')
data=[['H','Fixed state width','Moving state width','Fixed J upper','Moving J upper','Quadrature J upper']]
for r in rows[1:]:data.append([f"{float(r['time']):.1f}",f"{float(r['fixed_state_width']):.6f}",f"{float(r['moving_state_width']):.6f}",f"{float(r['fixed_J_upper']):.7f}",f"{float(r['moving_J_upper']):.7f}",f"{float(r['quadrature_J_upper']):.7f}"])
t=Table(data,colWidths=[35,94,96,88,92,98],repeatRows=1,hAlign='CENTER');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e6e9ec')),('FONTSIZE',(0,0),(-1,-1),7.3),('BOTTOMPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),7),('LINEBELOW',(0,0),(-1,0),.6,colors.grey),('LINEBELOW',(0,1),(-1,-1),.2,colors.lightgrey)]));story.append(t);story.append(Spacer(1,12))
p('<b>Result.</b> At H=2, the moving frame decreases the summed physical-coordinate widths from 0.06360427 to 0.02854866 (55.1%). Yet its expenditure upper bound 0.00231344 exceeds the fixed affine bound 0.00208117. The independently accumulated tube bound, 0.00254582, is weaker than both. Taking the minimum recovers 0.00208117 without producing a new improvement. The conservative uniform resource bound E(t) ≥ 0.1 - 0.002081171 &gt; 0.08 holds throughout [0,2] when E\' = 0.01 - 0.005 - J\'.')
p('<b>Interpretation.</b> Narrower state intervals do not automatically yield a tighter bound on a path-dependent quadratic observable. Tube-based interval quadrature discards temporal dependence, and intersecting projected tubes does not recover correlations between m(t) and accumulated J(t). A stronger observable-aware reconditioning would need to transport a joint correlated representation of (p,m,J), potentially with additional generators representing integrated quadratic terms. Merely choosing the narrowest state hull at each step is insufficient.')
p('<b>Scope and auditability.</b> The Python script computes the comparison using exact fractions, Picard tube inclusions, second-derivative flow-defect bounds and outward dyadic rounding. This is a finite-horizon certificate for the reduced mathematical model, not a biological viability claim. The underlying propagators are reused from the preceding fixed-affine and adaptive-frame implementations. The scalar quadrature certificate is independently derived from their validated tubes.')
doc.build(story)
print(out)
