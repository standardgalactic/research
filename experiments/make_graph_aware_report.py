from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib import colors
import csv
p='/mnt/data/graph_aware_frequency_report.pdf'; rows=list(csv.DictReader(open('/mnt/data/graph_aware_frequency_results.csv')))
s=getSampleStyleSheet();s.add(ParagraphStyle(name='Body2',parent=s['BodyText'],fontSize=10.2,leading=15,spaceAfter=10));s.add(ParagraphStyle(name='Head2',parent=s['Heading2'],spaceBefore=12,spaceAfter=7));doc=SimpleDocTemplate(p,pagesize=(595,842),rightMargin=54,leftMargin=54,topMargin=54,bottomMargin=54)
story=[]
def para(x):story.append(Paragraph(x,s['Body2']))
def head(x):story.append(Paragraph(x,s['Head2']))
story.append(Paragraph('Graph-Aware Frequency Forcing in Stuart–Landau Networks',s['Title']));para('Flyxion | Mathematical research supplement | November 2026')
para('This supplement refines the forced-disagreement bound for complex Stuart–Landau oscillators on a connected undirected graph. It uses the spatial pattern of frequency mismatch rather than only its largest component. The result is a sufficient bound, not an estimate of observed biological expenditure.')
head('Model and assumptions')
para('Let z_i be complex oscillator states satisfying z-dot = (a+i Omega)z - g(z) - kLz, where g(z)_i = |z_i|² z_i. Assume a is constant with a ≤ a_+, k > 0, L is the symmetric graph Laplacian, and |z_i(t)| ≤ R for every node and time. For a_+ ≥ 0, the last assumption follows if R ≥ max(max_i |z_i(0)|, sqrt(a_+)), by the scalar maximum principle. Set P=I-11ᵀ/N, y=Pz, omega-tilde = P omega, and delta = max_i |omega_i - mean(omega)|.')
head('Proposition: graph-resolved disagreement certificate')
para('Suppose gamma = k lambda_2 - a_+ - delta > 0. On the mean-zero subspace, define A = kL - (a_+ + delta)P, which is positive definite there. Let h = omega-tildeᵀ A† omega-tilde, with A† the inverse on that subspace. Then for all t ≥ 0,')
para('||Pz(t)||² ≤ exp(-gamma t)||Pz(0)||² + (R² h/gamma)(1-exp(-gamma t)).')
para('In particular, the asymptotic disagreement ceiling is R sqrt(h/gamma).')
head('Proof')
para('Write z = mean(z)1+y and remove the common rotating frequency. The cubic map g is monotone on the real two-dimensional representation of each complex coordinate, hence Re〈y,Pg(z)〉 ≥ 0. The heterogeneous-frequency term decomposes as iP Omega y + i mean(z) omega-tilde. Its first contribution to the disagreement-energy derivative is at most delta ||y||², while |mean(z)| ≤ R. Therefore, for V=||y||²/2, V-dot ≤ -〈y,Ay〉 + R|〈y,i omega-tilde〉|. Weighted Young inequality gives V-dot ≤ -〈y,Ay〉/2 + R² h/2 ≤ -gamma V + R² h/2. Scalar comparison proves the assertion. All inner products here use their real parts.')
head('Numerical comparison')
para('The four-node examples use a_+=0.19, k=0.7, R=0.45, frequencies (-0.03,-0.01,0.01,0.03), and conservative rational lower bounds for lambda_2. The quadratic form h is evaluated exactly by rational matrix inversion on the mean-zero subspace. Reported ceilings are rounded upward to 40-bit dyadic rationals before decimal display.')
data=[['Graph','Original ceiling','Graph-aware ceiling','Reduction']]
for r in rows:data.append([r['graph'].replace('_4',''),f"{float(r['original_disagreement_ceiling']):.8f}",f"{float(r['graph_aware_disagreement_ceiling']):.8f}",f"{float(r['graph_vs_original_pct']):.2f}%"])
t=Table(data,colWidths=[95,115,140,90],repeatRows=1);t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e9edf1')),('LINEBELOW',(0,0),(-1,0),.5,colors.grey),('BOTTOMPADDING',(0,0),(-1,-1),8),('TOPPADDING',(0,0),(-1,-1),8),('FONTSIZE',(0,0),(-1,-1),9)]));story.append(t)
head('Resource implication and limitations')
para('Since ||Lz||² ≤ lambda_max² ||Pz||², actuator power obeys P_A ≤ c_A k² lambda_max² ||Pz||². Put C=c_A k² lambda_max², U=R²h/gamma, and Y0=||Pz(0)||². Then J(H) ≤ C[U H+(Y0-U)(1-exp(-gamma H))/gamma]. This is a finite-horizon upper bound; its persistent-rate ceiling is C U. For a resource equation E-dot = s-P_A with constant net supply s, the sufficient all-horizon condition s ≥ C U yields E(t) ≥ E(0)-C max(Y0-U,0)/gamma. Other consumption terms must be included separately.')
para('The graph-aware ceiling is not guaranteed to improve on every other correlation-preserving certificate, and its advantage depends on frequency-mode alignment. The comparison uses conservative spectral bounds, not verified nonlinear trajectories. It assumes a common real growth rate, undirected coupling, bounded amplitudes, and a fixed frequency vector. The pseudoinverse is a graph-aware forcing metric; it is not an exact steady-state solution.')
doc.build(story);print(p)
