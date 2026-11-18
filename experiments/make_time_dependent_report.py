from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors
p='/mnt/data/time_dependent_dissipation_report.pdf';doc=SimpleDocTemplate(p,pagesize=(595,842),leftMargin=57,rightMargin=57,topMargin=55,bottomMargin=55)
s=getSampleStyleSheet(); story=[]
def para(x,style='BodyText'):story.append(Paragraph(x,s[style]));story.append(Spacer(1,10))
para('Time-Dependent Dissipation and a Quadratic Comparison Certificate','Title')
para('Reduced real Stuart–Landau dynamics; exact-rational inequality verification','Italic')
para('For p=(u+v)/2 and m=(u-v)/2, the equations are p\u2032=p(a-p²-3m²), m\u2032=-(2k-a+3p²+m²)m, J\u2032=Cm², with k=0.7, C=8c_Ak²=0.392, c_A=0.1, and a in [0.15,0.19] fixed on each trajectory. Initial intervals are p in [0.329,0.331] and m in [0.119,0.121].')
para('Invariant region and first dissipation bound','Heading2')
para('The region p >= 0.325, 0 <= m <= 0.121 is forward invariant. On p=0.325, a-p²-3m² >= 0.000452 > 0; m decreases for positive m. Thus D=2k-a+3p²+m² >= D0=1.526875. Writing V=m², V\u2032=-2(D0+V+R)V with R>=0, hence V(t) <= M exp(-2D0 t), M=0.121².')
para('A certified later lower bound on p','Heading2')
para('For t >= 0.5 and p in [0.325,0.329], p\u2032 >= 0.325[0.15-0.329²-3M exp(-D0)] > 0.0104709591. The minimum increase during the following half-unit of time is greater than 0.00523547, exceeding the strip width 0.004. A first-crossing argument therefore proves p(t) >= 0.329 for t >= 1. At its lower boundary the vector field points inward. Accordingly D >= D1=1.534723 for t >= 1.')
para('Two-stage expenditure certificate','Heading2')
para('For t >= 0, integrate V(t) <= M exp(-2D0 t) on [0,1], and V(t) <= M exp(-2D0) exp(-2D1(t-1)) afterward. This yields J(infinity) <= (CM/2)[(1-exp(-2D0))/D0 + exp(-2D0)/D1]. The implementation uses alternating exact-rational Taylor bounds for the exponential, rather than floating-point comparisons.')
para('Riccati sharpening','Heading2')
para('Because V\u2032 <= -2(D0+V)V, comparison with the scalar Riccati equation yields J(infinity) <= (C/2) log(1+M/D0). The script bounds the logarithm above by the 15-term alternating rational series; M/D0 < 1. This bound is independent of a finite integration horizon.')
rows=[['Certificate','J upper','Energy lower'],['Constant D0','0.001879417765043','0.098120582234957'],['Two-stage D(t)','0.001878964319073','0.098121035680927'],['Riccati V','0.001870464213121','0.098129535786879']]
t=Table(rows,colWidths=[155,150,150]);t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e6ebf0')),('GRID',(0,0),(-1,-1),.35,colors.grey),('FONTSIZE',(0,0),(-1,-1),9),('BOTTOMPADDING',(0,0),(-1,-1),8),('TOPPADDING',(0,0),(-1,-1),8)]));story.append(t);story.append(Spacer(1,14))
para('Resource interpretation and scope','Heading2')
para('For E\u2032=0.01-0.005-J\u2032 and E(0)=0.1, positivity of the net supply implies E(t) >= 0.1-J(infinity). The strongest bound here gives E(t) >= 0.098129535786879 > 0.08 for all t >= 0. This is a certificate for the reduced real oscillator system, not a claim about biological calibration. The two-stage and Riccati certificates are separately valid; their numerical bounds should not be added together.')
para('Reproduction: python time_dependent_dissipation.py. All mathematical comparisons and transcendental upper bounds are evaluated using fractions.Fraction. Decimal values in the table are display approximations.','Italic')
doc.build(story)
print(p)
