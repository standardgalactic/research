import csv
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,KeepTogether
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.units import inch
from pathlib import Path
P=Path('/mnt/data');rows=list(csv.DictReader(open(P/'direct_spectral_expenditure_results.csv')))
styles=getSampleStyleSheet();styles.add(ParagraphStyle(name='Title2',parent=styles['Title'],fontSize=15,leading=19,spaceAfter=14));styles.add(ParagraphStyle(name='Body2',parent=styles['BodyText'],fontSize=9.5,leading=14,spaceAfter=10));styles.add(ParagraphStyle(name='Eq2',parent=styles['BodyText'],fontName='Courier',fontSize=9,leading=14,leftIndent=15,spaceAfter=12));styles.add(ParagraphStyle(name='Small2',parent=styles['BodyText'],fontSize=8,leading=11))
doc=SimpleDocTemplate(str(P/'direct_spectral_expenditure_report.pdf'),pagesize=(612,792),leftMargin=56,rightMargin=56,topMargin=50,bottomMargin=50)
story=[]
def p(t):story.append(Paragraph(t,styles['Body2']))
def eq(t):story.append(Paragraph(t,styles['Eq2']))
story.append(Paragraph('Direct Spectral Certification of Heterogeneous-Network Actuator Expenditure',styles['Title2']))
p('A network-level energy estimate avoids the conversion of a disagreement-norm ceiling into an actuator-power ceiling using the square of the largest Laplacian eigenvalue. The result applies to the full complex Stuart–Landau equations, rather than to a two-oscillator reduction.')
p('<b>Model and assumptions.</b> Consider z-dot = (a + i Omega)z - g(z) - kLz, with g(z)_i = |z_i|^2 z_i, connected undirected graph Laplacian L, uniform real growth rate a, k &gt; 0, and fixed intrinsic frequencies omega_i. Let omega-tilde = omega - mean(omega)1, delta = max |omega-tilde_i|. Assume |z_i(0)| &lt;= R and R^2 &gt;= max(a,0). The maximum-amplitude comparison principle then ensures |z_i(t)| &lt;= R. Let lambda2_lower be a certified lower bound on the positive Laplacian spectral gap.')
p('<b>Proposition (direct expenditure certificate).</b> Set W = z*Lz/2, Q = ||Lz||^2, beta = k - (a+delta)/lambda2_lower, and F = R||omega-tilde||. If beta &gt; 0, then for every T &gt;= 0,')
eq('W-dot &lt;= -beta Q + F sqrt(Q),')
eq('J(T) = c_A k^2 integral_0^T Q(t)dt')
eq('     &lt;= 2 c_A k^2 W(0)/beta + c_A k^2 F^2 T/beta^2.')
p('<b>Proof.</b> Differentiating W gives Re[(Lz)* z-dot]. The real growth term is a z*Lz, bounded by (a/lambda2_lower)Q. The coupling contributes -kQ. The cubic term is nonpositive because, for each edge, Re[(z_i-z_j)* (g(z_i)-g(z_j))] &gt;= 0 by monotonicity of the gradient of |z|^4/4 on the real plane. Subtract the common mean frequency, whose contribution to W-dot vanishes. Writing z = Pz + mean(z)1, the remaining frequency term is at most delta||Lz||||Pz|| + R||omega-tilde||||Lz||, bounded by (delta/lambda2_lower)Q + F sqrt(Q). Young’s inequality F sqrt(Q) &lt;= beta Q/2 + F^2/(2beta) and W(T) &gt;= 0 yield the claimed integral inequality.')
p('<b>Resource consequence.</b> Suppose E-dot = s - c_A k^2 Q with constant net supply s = q-b. If s &gt;= c_A k^2 F^2/beta^2, then E(t) &gt;= E(0) - 2 c_A k^2 W(0)/beta for every t &gt;= 0. This is a sufficient all-time floor certificate; it does not assert convergence of actuator power to the displayed rate.')
data=[['Topology','Transient buffer','New rate','Prior rate','Energy floor']]
for r in rows:data.append([r['graph'].replace('_4',''),f"{float(r['direct_transient_expenditure_buffer']):.6f}",f"{float(r['direct_persistent_power_ceiling']):.7f}",f"{float(r['previous_graph_aware_power_ceiling']):.7f}",f"{float(r['resource_floor_certified']):.6f}"])
t=Table(data,colWidths=[90,100,95,95,110],repeatRows=1);t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e5eaf0')),('GRID',(0,0),(-1,-1),.35,colors.grey),('VALIGN',(0,0),(-1,-1),'MIDDLE'),('FONTSIZE',(0,0),(-1,-1),8),('BOTTOMPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),7)]));story.append(t)
story.append(Spacer(1,12))
p('<b>Reproducibility and scope.</b> Calculations use Python fractions for all displayed rational bounds and SymPy rational Laplacians. Inputs are n=4, a=0.19, k=0.7, c_A=0.1, R=0.45, omega=(-0.03,-0.01,0.01,0.03), initial real z=(0.45,0.21,0.39,0.27), s=0.005 and E(0)=0.1. The path uses conservative lambda2_lower=0.585, cycle=2, complete=4. The table compares certified upper bounds, not measured actuator expenditure. The transient bound is deliberately conservative and can be optimized by adjusting the Young parameter. Frequency mismatch remains an external work term and no exact phase-locked state is assumed.')
doc.build(story)
print(P/'direct_spectral_expenditure_report.pdf')
