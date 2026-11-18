from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from fractions import Fraction as F
from pathlib import Path
out='/mnt/data/joint_dissipation_report.pdf'
styles=getSampleStyleSheet(); styles.add(ParagraphStyle(name='TitleCustom',parent=styles['Title'],fontSize=15,leading=19,spaceAfter=15));styles.add(ParagraphStyle(name='BodyCustom',parent=styles['BodyText'],fontSize=10.4,leading=15,spaceAfter=9));styles.add(ParagraphStyle(name='EquationCustom',parent=styles['BodyText'],fontName='Courier',fontSize=9,leading=14,leftIndent=17,spaceAfter=10))
doc=SimpleDocTemplate(out,pagesize=(595,842),leftMargin=58,rightMargin=58,topMargin=55,bottomMargin=55)
story=[]
def p(t,style='BodyCustom'):story.append(Paragraph(t,styles[style]))
p('A Joint Dissipation Certificate for Oscillator Disagreement and Expenditure','TitleCustom')
p('Research supplement | 18 November 2026')
p('<b>Model and objective.</b> For the real two-oscillator Stuart–Landau subsystem in collective/disagreement coordinates, with persistent a in [0.15, 0.19], k=0.7 and c_A=0.1, consider')
p('p\' = ap - p^3 - 3pm^2','EquationCustom')
p('m\' = -(2k-a+3p^2+m^2)m,   J\' = 8c_A k^2 m^2.','EquationCustom')
p('Initial uncertainty is p(0) in [0.329, 0.331], m(0) in [0.119, 0.121], and J(0)=0. The goal is an upper bound on cumulative expenditure J(t), uniformly in time and across the persistent parameter uncertainty.')
p('<b>Invariant region.</b> The half-strip p >= 0.325 and 0 <= m <= 0.121 is positively invariant for the relevant initial states. At p=0.325, using a >= 0.15 and m <= 0.121,')
p('p\'/p >= 0.15 - 0.325^2 - 3(0.121)^2 = 0.000452 > 0.','EquationCustom')
p('At m=0.121, the disagreement derivative is negative because 2k-a+3p^2+m^2 is positive. At m=0 it vanishes. Thus the region cannot be exited through either boundary. Polynomial vector fields are locally Lipschitz; the dissipative cubic terms and boundedness of m prevent finite-time escape of p. The same argument applies to every fixed a in the prescribed interval.')
p('<b>Joint differential inequality.</b> Put V=m^2 and D=2k-a+3p^2+m^2. The invariant region implies')
p('D >= 1.4 - 0.19 + 3(0.325)^2 = 1.526875 = 2443/1600.','EquationCustom')
p('Since V\'=-2DV and J\'=C V, where C=8c_A k^2=0.392, the combined observable W=J+(C/(2D_min))V satisfies')
p('W\' = C V(1 - D/D_min) <= 0.','EquationCustom')
p('Hence, for every t >= 0, J(t)+(C/(2D_min))m(t)^2 <= (C/(2D_min))m(0)^2. Dropping the nonnegative terminal term and maximizing over the initial interval gives the exact rational bound')
p('J(t) <= 0.392(0.121)^2 / (2 x 1.526875) = 0.001879417765043.','EquationCustom')
p('<b>Comparison and resource floor.</b> The previous validated fixed-affine enclosure at H=2 gave J(2) <= 0.002081170556266443. The joint dissipation certificate improves that upper bound by about 9.69%. If E\'=q-b-J\' with E(0)=0.1, q=0.01 and b=0.005, then E(t)=0.1+0.005t-J(t) >= 0.098120582234957 for all t >= 0, exceeding the floor E*=0.08.')
data=[['Method','Expenditure upper bound'],['Fixed affine, H=2','0.002081170556266'],['Joint dissipation, all t >= 0','0.001879417765043']]
t=Table(data,colWidths=[255,190]);t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e8edf1')),('GRID',(0,0),(-1,-1),0.4,colors.grey),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),('FONTSIZE',(0,0),(-1,-1),9),('BOTTOMPADDING',(0,0),(-1,-1),8),('TOPPADDING',(0,0),(-1,-1),8)]));story.append(t)
story.append(Spacer(1,14))
p('<b>Scope.</b> This is an analytic continuous-time certificate, not a numerical Lohner integration and not a claim about a calibrated biological system. It exploits the exact quadratic coupling between disagreement decay and actuator expenditure. The uniform bound is conservative because it replaces the state-dependent dissipation D(t) by a constant lower bound. It can be sharpened with certified lower bounds on p(t) or piecewise estimates of D(t).')
doc.build(story)
print(out)
