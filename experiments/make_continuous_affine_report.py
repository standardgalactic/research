from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
P='/mnt/data/continuous_affine_extension.pdf'
doc=SimpleDocTemplate(P,pagesize=(595,842),leftMargin=62,rightMargin=62,topMargin=60,bottomMargin=55)
s=getSampleStyleSheet();s.add(ParagraphStyle(name='MyTitle',parent=s['Title'],fontSize=16,leading=21,spaceAfter=17));s.add(ParagraphStyle(name='BodyX',parent=s['BodyText'],fontSize=10.2,leading=15,spaceAfter=11));s.add(ParagraphStyle(name='SectionX',parent=s['Heading2'],fontSize=12,spaceBefore=11,spaceAfter=7))
story=[]
def p(x,style='BodyX'):story.append(Paragraph(x,s[style]))
p('Picard-Validated Continuous-Time Moving-Affine Enclosures for Nonlinear Oscillator Actuation','MyTitle')
p('Flyxion | Research supplement | November 2026')
p('Abstract','SectionX')
p('We extend a dependency-preserving exact-rational enclosure of a nonlinear discrete oscillator map to the continuous-time differential equations by adding a locally verified Picard tube and a second-derivative truncation enclosure. The calculation includes accumulated actuator expenditure. The resulting affine enclosure is narrower than direct interval propagation on a short benchmark, but does not yet establish a long-horizon or full complex-network certificate.')
p('Model and parameter dependence','SectionX')
p('The real two-oscillator subsystem is x1\' = a x1 - x1^3 - k(x1-x2), x2\' = a x2 - x2^3 - k(x2-x1), J\' = 2 cA k^2 (x1-x2)^2. The initial state is (0.45, 0.21, 0); k = 0.7, cA = 0.1, and a = 0.17 + 0.02 xi for a shared xi in [-1,1]. The interval of possible a values is [0.15,0.19]. The same xi is retained at every step.')
p('Certified local time remainder','SectionX')
p('Let X enclose the initial state at one step and Y be an interval tube satisfying X + h f(Y) contained in Y, with strict inflation available in the implementation. A Picard argument encloses every solution over that step inside Y. Since the system is autonomous, y\'\' = Df(y) f(y). Therefore the exact flow from any initial y in X lies in y + h f(y) + (h^2/2) Df(Y) f(Y), where all interval operations enclose their corresponding pointwise operations. The implementation uses exact rational endpoints for these operations and verifies the tube inclusion before advancing.')
p('Moving-affine propagation','SectionX')
p('Write y = c + A xi + e, with coordinatewise |e| bounded by r. At each step, the center follows a nominal Euler update, the shared parameter coefficient follows the first-order sensitivity update, and the remainder encloses nonlinear cross terms, existing remainder transport, and the certified continuous-time truncation error. For the state coordinates, the nonlinear residual includes d U + 3|c| U^2 + U^3, where U = |A| + r and d = 0.02. For expenditure, the quadratic residual includes 2 cA k^2 (U1+U2)^2. These are inclusion bounds rather than fitted error estimates.')
p('Numerical result','SectionX')
data=[['Horizon','Direct state width','Affine state width','Direct J upper','Affine J upper'],['0.025','0.00109885','0.00109885','0.00014794','0.00014794'],['0.050','0.00223689','0.00219577','0.00028517','0.00028493'],['0.075','0.00341667','0.00329203','0.00041253','0.00041180'],['0.100','0.00464085','0.00438884','0.00053082','0.00052934']]
t=Table(data,colWidths=[62,100,100,103,103],repeatRows=1);t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e9edf2')),('GRID',(0,0),(-1,-1),0.3,colors.grey),('FONTSIZE',(0,0),(-1,-1),8),('BOTTOMPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),7)]));story.append(t);story.append(Spacer(1,12))
p('Resource implication and limitations','SectionX')
p('If resource is E(t) = E0 - J(t), the certified upper bound on J gives the sufficient condition E0 - E_floor greater than or equal to J_upper(H) for this reduced model, since J is nondecreasing. The bound is not a biological energy calibration. The test is limited to four steps of length 0.025, a real invariant subsystem, and a single uncertain growth parameter. It does not cover complex phase dynamics, frequency heterogeneity, or long-horizon viability. Exact-rational arithmetic is costly because polynomial denominators grow rapidly; practical extension requires rigorously rounded finite-precision arithmetic and more sophisticated wrapping control.')
p('Reproducibility','SectionX')
p('The accompanying continuous_affine_validated.py implements rational interval operations, Picard-tube verification, second-derivative enclosures, and a direct-box comparator. The CSV records all five time checkpoints. This supplement claims validated mathematical enclosures under the implemented model and arithmetic, not a general-purpose formally verified ODE solver.')
doc.build(story);print(P)
