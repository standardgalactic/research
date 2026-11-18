from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import inch
from pathlib import Path
out='/mnt/data/piecewise_riccati_report.pdf'
doc=SimpleDocTemplate(out,pagesize=(612,792),leftMargin=64,rightMargin=64,topMargin=60,bottomMargin=55)
s=getSampleStyleSheet();s.add(ParagraphStyle(name='TitleX',parent=s['Title'],fontSize=16,leading=20,spaceAfter=18));s.add(ParagraphStyle(name='BodyX',parent=s['BodyText'],fontSize=10.5,leading=16,spaceAfter=12));s.add(ParagraphStyle(name='HeadX',parent=s['Heading2'],spaceBefore=10,spaceAfter=9))
story=[]
def p(t):story.append(Paragraph(t,s['BodyX']))
story.append(Paragraph('Piecewise Riccati Dissipation for Oscillator Expenditure',s['TitleX']))
p('Exact-rational comparison certificate for the reduced real Stuart–Landau oscillator pair. This note extends the prior constant-Riccati estimate by combining it with a previously established stronger lower bound on collective-state dissipation after time t = 1.')
story.append(Paragraph('Model and assumptions',s['HeadX']))
p('With p = (u+v)/2, m = (u-v)/2, a in [0.15, 0.19] fixed on each trajectory, k = 0.7 and c_A = 0.1, the equations are p\' = ap - p³ - 3pm², m\' = -(2k-a+3p²+m²)m, and J\' = 8c_Ak²m². The initial squared disagreement satisfies V(0) = m(0)² &lt;= (0.121)². The previously established invariant-region and comparison arguments give p &gt;= 0.325 for all t &gt;= 0 and p &gt;= 0.329 for t &gt;= 1. These two premises are inherited, not reproved here.')
story.append(Paragraph('Comparison and integration',s['HeadX']))
p('Set V = m² and C = 8c_Ak² = 0.392. Then V\' &lt;= -2(D(t)+V)V, where D(t) = D0 = 1.526875 for 0 &lt;= t &lt; 1, and D(t) = D1 = 1.534723 for t &gt;= 1. For constant D and starting value v, the maximal comparison solution is V_D(t;v) = D v exp(-2Dt)/(D+v[1-exp(-2Dt)]). Its time integral from zero to T equals one half times log((D+v)/(D+V_D(T;v))).')
p('Let v1 = V_D0(1;V0), where V0 = (0.121)². The total all-horizon expenditure is bounded by J_infty &lt;= (C/2) [log(1+V0/D0) - log(1+v1/D0) + log(1+v1/D1)]. The last two terms have a negative sum because D1 &gt; D0. This is strictly sharper than the constant-D0 Riccati bound. The comparison is valid on the entire trajectory, not only at the grid points.')
story.append(Paragraph('Numerical certification',s['HeadX']))
data=[['Method','J upper bound','Uniform E lower bound'],['Constant linear','0.001879417765043','0.098120582234957'],['Constant Riccati','0.001870464213121','0.098129535786879'],['Piecewise Riccati','0.001870015073372','0.098129984926628']]
t=Table(data,colWidths=[140,145,170],repeatRows=1);t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#eaeef2')),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),('FONTSIZE',(0,0),(-1,-1),9),('BOTTOMPADDING',(0,0),(-1,-1),8),('TOPPADDING',(0,0),(-1,-1),8),('LINEBELOW',(0,0),(-1,0),.6,colors.grey),('LINEBELOW',(0,-1),(-1,-1),.4,colors.grey)]));story.append(t);story.append(Spacer(1,13))
p('The exact-rational verification script bounds exp(-2D0) by alternating Taylor sums and bounds log(1+x) with alternating odd/even truncations. It uses the lower exponential enclosure for v1, which maximizes the upper-bound expression because its correction decreases with v1. The resulting bound improves the constant-Riccati certificate by approximately 0.024012%.')
story.append(Paragraph('Resource-floor implication and scope',s['HeadX']))
p('If E(t) = 0.1 + (q-b)t - J(t) with q-b = 0.005 &gt;= 0 and J nondecreasing, then E(t) &gt;= 0.1-J_infty &gt;= 0.098129984926628 &gt; 0.08 for all t &gt;= 0. This is an analytic comparison certificate for the reduced real-state model. It is not a calibrated biological prediction, and it relies on the invariant-region premises above. The displayed decimal bounds are rounded for presentation; all comparisons in the accompanying code use exact fractions.')
doc.build(story)
print(out,Path(out).stat().st_size)
