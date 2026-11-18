from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib import colors
import csv
rows=list(csv.DictReader(open('/mnt/data/moving_frame_reconditioning_results.csv')))
p='/mnt/data/moving_frame_reconditioning_report.pdf'
doc=SimpleDocTemplate(p,pagesize=(612,792),leftMargin=62,rightMargin=62,topMargin=62,bottomMargin=58)
s=getSampleStyleSheet();s.add(ParagraphStyle(name='Body2',parent=s['BodyText'],leading=15,spaceAfter=11));s.add(ParagraphStyle(name='Title2',parent=s['Title'],fontSize=16,leading=20,spaceAfter=16))
story=[Paragraph('Fixed and Time-Dependent Frames for Validated Stuart–Landau Enclosures',s['Title2']),Paragraph('A controlled nonlinear comparison | November 2026',s['Italic'])]
paras=[
'<b>System and common uncertainty.</b> Consider the reduced real two-oscillator Stuart–Landau system in collective and disagreement coordinates p=(u+v)/2 and m=(u-v)/2. Its dynamics are p\u2032=ap-p³-3pm², m\u2032=(a-2k)m-3p²m-m³, and J\u2032=8 c_A k² m². The uncertain growth parameter is fixed in a in [0.15,0.19], with k=0.7, c_A=0.1, p(0) in [0.329,0.331], m(0) in [0.119,0.121], and J(0)=0. Both calculations use these identical correlated initial conditions.',
'<b>Moving coordinate construction.</b> The fixed method propagates (p,m,J). The moving method writes m=s_n v on each step, with a positive rational scale s_n. During a step, s_n is constant and v\u2032=(a-2k)v-3p²v-s_n²v³. At the next checkpoint, s_(n+1)=(1-2kh)s_n, and the interval for v is multiplied by s_n/s_(n+1). This is an invertible change of coordinates; it does not discard or add admissible physical states except through outward rounding.',
'<b>Continuous-time inclusion.</b> For each interval X and parameter box A, an inflated interval Y is constructed until X+[0,h]f(Y,A) is contained in Y. The next enclosure is outward-rounded from X+h f(Y,A). This encloses the continuous flow for each fixed parameter in A, conditional on the standard local Picard inclusion theorem. All interval endpoints, containment comparisons, and rescalings are computed with exact rational arithmetic; endpoints are rounded outward to 40 dyadic bits.',
'<b>Interpretation.</b> The moving coordinate is aligned with linear disagreement decay, but does not change the fundamental coordinatewise enclosure of nonlinear remainder terms. The observed improvement is very small, so the experiment does not establish a practically significant advantage for this elementary reconditioning scheme. More general matrix-valued frames and explicit remainder correlations remain untested.',
'<b>Resource-floor statement.</b> With E(t)=0.1-0.005t-J(t) and E*=0.08, monotonicity of J and the interval expenditure upper bound at the horizon imply E(t) is at least the reported conservative lower bound for all earlier times. The model is synthetic, with no calibrated biological expenditure interpretation.'
]
for para in paras[:3]:story.append(Paragraph(para,s['Body2']))
data=[['H','Fixed width','Moving width','Fixed J upper','Moving J upper']]
for r in rows[1:]:data.append([f"{float(r['horizon']):.1f}",f"{float(r['fixed_physical_width']):.6f}",f"{float(r['moving_physical_width']):.6f}",f"{float(r['fixed_J_upper']):.7f}",f"{float(r['moving_J_upper']):.7f}"])
t=Table(data,colWidths=[38,105,105,112,112],repeatRows=1);t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#eeeeee')),('GRID',(0,0),(-1,-1),.3,colors.grey),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),('FONTSIZE',(0,0),(-1,-1),8),('BOTTOMPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),7)]));story.append(t);story.append(Spacer(1,14))
for para in paras[3:]:story.append(Paragraph(para,s['Body2']))
story.append(Paragraph('<b>Reproducibility.</b> The companion Python script generates the CSV using Fraction-based exact rational intervals. The reported results cover 80 steps of length h=0.01.',s['Body2']))
doc.build(story)
print(p)
