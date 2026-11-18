from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,KeepTogether
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
import csv
p='/mnt/data/continuous_lohner_report.pdf';doc=SimpleDocTemplate(p,pagesize=(595,842),rightMargin=58,leftMargin=58,topMargin=55,bottomMargin=50)
s=getSampleStyleSheet();s.add(ParagraphStyle(name='BodyX',parent=s['BodyText'],fontSize=10.2,leading=15,spaceAfter=11));s.add(ParagraphStyle(name='TitleX',parent=s['Title'],fontSize=16,leading=20,spaceAfter=18));s.add(ParagraphStyle(name='SubX',parent=s['Heading2'],fontSize=12,spaceBefore=12,spaceAfter=7));story=[]
def P(t,sty='BodyX'):story.append(Paragraph(t, s[sty]))
P('Continuous-Time Validation of Correlated-Generator Stuart–Landau Enclosures','TitleX')
P('A reproducible comparison of direct interval and affine-generator propagation for a nonlinear oscillator flow. All numerical interval endpoints are represented by exact rational fractions and rounded outward to a dyadic grid of 40 fractional bits.')
P('Model and uncertainty','SubX')
P('The real two-oscillator system is expressed in collective and disagreement coordinates p=(u+v)/2 and m=(u-v)/2. Its equations are p\' = ap - p³ - 3pm²; m\' = (a-2k)m - 3p²m - m³; J\' = 8c<sub>A</sub>k²m². The parameter a is constant along each trajectory and ranges over [0.15, 0.19]. We use k=0.7, c<sub>A</sub>=0.1, p(0) in [0.329, 0.331], m(0) in [0.119, 0.121], J(0)=0, and h=0.01.')
P('Validated flow step','SubX')
P('For an enclosure X of the current state, a tube T is found with the strict inclusion X + h f(T,[a]) inside the interior of T. The interval Picard existence argument then encloses every solution over the complete step. Since a is constant in time, x\'\' = D<sub>x</sub>f(x,a) f(x,a). Thus the exact flow endpoint belongs to X + h f(X,[a]) + [-D,D], with D<sub>i</sub> = (h²/2) sum<sub>j</sub> sup |(D<sub>x</sub>f)<sub>ij</sub> f<sub>j</sub>| evaluated over T. The same validated defect is added to the correlated affine propagator.')
P('Correlated affine enclosure','SubX')
P('Represent x=c+sum g<sub>l</sub> ξ<sub>l</sub>+[-r,r], where independent ξ<sub>l</sub> range over [-1,1]. Two generators encode initial-state uncertainty and a third encodes the persistent parameter a=0.17+0.02 ξ<sub>3</sub>. The center follows a first-order flow step, while generators follow the Jacobian linearization with the parameter-forcing term. A componentwise remainder includes the old-remainder transport, an interval bound on Jacobian variation over X, the parameter-state cross term, the validated flow defect, and outward rounding. Consequently the resulting enclosure contains the exact flow endpoint, although it is not necessarily a minimal affine enclosure.')
P('Numerical results','SubX')
rows=list(csv.DictReader(open('/mnt/data/continuous_lohner_results.csv')))
data=[['Horizon','Interval width','Affine width','J upper (affine)']]
for r in rows[1:]:data.append([f"{float(r['time']):.1f}",f"{float(r['box_state_width']):.6f}",f"{float(r['affine_state_width']):.6f}",f"{float(r['affine_J_upper']):.8f}"])
t=Table(data,colWidths=[75,118,118,130],hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#eaeef2')),('GRID',(0,0),(-1,-1),.3,colors.grey),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),('FONTSIZE',(0,0),(-1,-1),9),('BOTTOMPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),7)]));story.append(t)
P('Resource-floor certificate','SubX')
P('For E(t)=0.1+(0.01-0.005)t-J(t), the computed affine expenditure upper bound at t=0.8 is 0.0017792040. Because J\' = 8c<sub>A</sub>k²m² is nonnegative, the bound J(t) ≤ J(0.8) holds throughout [0,0.8]. Since the net supply rate is positive, E(t) ≥ 0.1-J(0.8) &gt; 0.09822079 &gt; 0.08 on the entire interval. This conservative bound applies between grid points as well as at the sampled horizons.')
P('Scope and limitations','SubX')
P('This is a continuous-time validated affine-generator benchmark for the reduced real-state system. It does not certify the full complex Stuart–Landau network, uncertain time-dependent parameters, or a biological model. The frame itself is not adaptively reconditioned; the result tests correlation retention plus validated nonlinear remainder propagation. A full Lohner implementation would additionally evolve and recondition an invertible state-frame matrix and validate its transformed remainder.')
doc.build(story)
print(p)
