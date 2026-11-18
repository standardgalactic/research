from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
import csv
p='/mnt/data/lohner_generator_report.pdf'; doc=SimpleDocTemplate(p,pagesize=(595,842),leftMargin=62,rightMargin=62,topMargin=55,bottomMargin=55)
s=getSampleStyleSheet();s.add(ParagraphStyle(name='T',parent=s['Title'],fontSize=15,leading=19,spaceAfter=16));s.add(ParagraphStyle(name='B',parent=s['BodyText'],fontSize=10,leading=15,spaceAfter=10));story=[]
def para(t):story.append(Paragraph(t,s['B']))
story.append(Paragraph('Correlated-Generator Propagation for a Nonlinear Oscillator Euler Map',s['T']))
para('Flyxion | Research supplement | November 2026')
para('<b>Scope.</b> We compare a conventional interval extension of a nonlinear Euler map against a correlated-generator affine enclosure, with a bounded nonlinear residual. The comparison is for the discrete Euler map only. It does not certify the continuous-time differential equation.')
para('<b>Model.</b> In collective/disagreement coordinates, p=(u+v)/2 and m=(u-v)/2, the real two-oscillator equations are p\' = ap - p³ - 3pm², m\' = (a-2k)m - 3p²m - m³, J\' = 8 cA k²m². We take k=0.7, cA=0.1, a in [0.15,0.19], initial p in [0.329,0.331], initial m in [0.119,0.121], and J(0)=0.')
para('<b>Affine construction.</b> Each discrete state is represented by y=c+G ξ+e, where the same three independent generators encode initial p, initial m, and the fixed uncertain a. The remainder e belongs to a coordinate interval. The Jacobian at the center transports G, while a mean-value bound on Jacobian variation bounds nonlinear remainder accumulation. Exact rational coefficients and outward dyadic rounding at 40 bits enclose all updates of this discrete map.')
rows=list(csv.DictReader(open('/mnt/data/lohner_generator_results.csv')))
data=[['Euler time','Interval width','Affine width','Interval J upper','Affine J upper']]
for r in rows:data.append([f"{float(r['time']):.2f}",f"{float(r['box_state_width']):.6f}",f"{float(r['affine_state_width']):.6f}",f"{float(r['box_J_upper']):.7f}",f"{float(r['affine_J_upper']):.7f}"])
t=Table(data,colWidths=[72,92,92,108,106]);t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e9edf2')),('GRID',(0,0),(-1,-1),.3,colors.grey),('FONTSIZE',(0,0),(-1,-1),8),('ALIGN',(0,0),(-1,-1),'CENTER'),('BOTTOMPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),7)]));story.append(t);story.append(Spacer(1,15))
para('<b>Result.</b> At t=0.4, the affine enclosure is 37.7% narrower in combined state-coordinate width than the independent interval enclosure. Its upper expenditure bound is also smaller. Both methods begin with the same initial set; the improvement is therefore not due to a different initial uncertainty hull.')
para('<b>Proof principle.</b> For a differentiable map T and a centered uncertain set, T(c+v)=T(c)+DT(c)v+R, where |R| is bounded componentwise by a verified Jacobian-variation bound times |v|. The affine update retains the linear generator contribution, and adds the bounded remainder to e. Every coefficient-rounding error is included in e. Induction gives inclusion for all iterates of the discrete map.')
para('<b>Limitations.</b> This is not a fully implemented moving-frame Lohner algorithm: no adaptive matrix reconditioning is used, and the enclosure is not yet a validated continuous-time flow. The reported J values bound the discrete accumulator, not continuous integral expenditure. A Picard-validated step with a time Taylor remainder is required before inferring resource-floor viability for the differential equations.')
doc.build(story)
print(p)
