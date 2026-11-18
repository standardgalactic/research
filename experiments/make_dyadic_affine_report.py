from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,KeepTogether
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.units import inch
import csv
p='/mnt/data/dyadic_affine_combined.pdf';s=getSampleStyleSheet();s.add(ParagraphStyle(name='TitleCustom',parent=s['Title'],fontSize=16,leading=20,spaceAfter=12));s.add(ParagraphStyle(name='BodyCustom',parent=s['BodyText'],fontSize=10,leading=15,spaceAfter=9));s.add(ParagraphStyle(name='HeadingCustom',parent=s['Heading2'],fontSize=12,spaceBefore=12,spaceAfter=6));doc=SimpleDocTemplate(p,pagesize=(612,792),leftMargin=60,rightMargin=60,topMargin=54,bottomMargin=50)
story=[]
def P(t):story.append(Paragraph(t,s['BodyCustom']))
def H(t):story.append(Paragraph(t,s['HeadingCustom']))
story.append(Paragraph('Dyadic Compression with Parameter-Affine Continuous-Time Enclosures',s['TitleCustom']))
P('Flyxion | Mathematical and computational research supplement | November 2026')
P('<b>Abstract.</b> We combine rational interval integration, outward dyadic rounding, and a parameter-affine decomposition for a real invariant two-oscillator Stuart–Landau subsystem. A common growth parameter is fixed across time but belongs to a continuous prior interval. The comparison measures enclosure widths, cumulative actuator expenditure, and a conservative all-time resource-floor condition. The experiment is restricted to a reduced real subsystem; it is not a certificate for the full complex network.')
H('Model and enclosure')
P('The reduced dynamics are u\u0307 = au - u³ - k(u-v), v\u0307 = av - v³ - k(v-u), J\u0307 = 2 c<sub>A</sub> k² (u-v)². We use a in [0.15,0.19], k=0.7, c<sub>A</sub>=0.1, (u(0),v(0),J(0))=(0.45,0.21,0), and step size h=0.025. The shared uncertain parameter is represented by a=a<sub>0</sub>+d xi with a<sub>0</sub>=0.17, d=0.02 and |xi|≤1.')
P('At each step the parameter-affine representation is y=c+s xi+e, with coordinatewise |e|≤r. The center and first-order coefficient are transported by the nominal vector field and Jacobian. A remainder majorant accounts for nonlinear terms, parameter–error products, and the continuous-time Taylor defect. All inclusion comparisons are performed using rational intervals. Center and sensitivity coefficients are dyadically rounded, and their rounding errors are added to the remainder before its upper endpoint is rounded outward.')
P('For the continuous-time defect, a local tube Y is accepted only after checking X+[0,h]f(Y,A) is contained in Y. On that tube the second derivative is enclosed by Df(Y,A) f(Y,A), yielding the step remainder h² sup|Df f|/2. This is a sufficient local inclusion construction, not an assertion of minimality.')
H('Numerical comparison')
rows=list(csv.DictReader(open('/mnt/data/dyadic_affine_combined_results.csv')))
data=[['H','Direct width','Affine width','Reduction','Affine J upper']]
for r in rows:
 if r['bits']!='40':continue
 a=float(r['affine_state_width']);b=float(r['box_state_width']);data.append([f"{float(r['horizon']):.1f}",f'{b:.6f}',f'{a:.6f}',f'{100*(1-a/b):.1f}%',f"{float(r['affine_J_upper']):.8f}"])
t=Table(data,colWidths=[45,105,105,85,115]);t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e8edf2')),('LINEBELOW',(0,0),(-1,0),0.6,colors.grey),('BOTTOMPADDING',(0,0),(-1,-1),8),('TOPPADDING',(0,0),(-1,-1),8),('FONTSIZE',(0,0),(-1,-1),9),('ALIGN',(1,1),(-1,-1),'RIGHT')]));story.append(t)
H('Resource-floor condition')
P('With E(0)=0.1, supply q=0.01, basal demand b=0.005, and E\u0307=q-b-J\u0307, a sufficient certificate for every t in [0,H] is E(0)-bH-J<sub>upper</sub>(H) ≥ E<sub>*</sub>, since supply is nonnegative and J is nondecreasing. At H=0.4, the 40-bit affine computation gives this conservative lower bound as 0.09657608, exceeding the stipulated floor E<sub>*</sub>=0.08. The stronger endpoint-only lower bound, which credits supply, is 0.10057608.')
H('Scope and limitations')
P('The reported guarantee is conditional on the algebraic remainder enclosure and local Picard hypotheses used by the implementation. The comparison is performed with exact rational arithmetic for interval endpoints and checks; its output CSV is decimal formatting of those rational results. This implementation uses a simple axis-aligned remainder and does not yet transport a fully general Lohner frame. It also does not establish parameter identification, full-network viability, or calibrated physiological energy expenditure. The natural next test is a longer horizon with adaptive reconditioning of the affine frame and explicit checks of the remainder growth.')
doc.build(story)
print(p)
