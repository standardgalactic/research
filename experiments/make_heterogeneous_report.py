from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,KeepTogether
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
import csv
P='/mnt/data/'
rows=list(csv.DictReader(open(P+'heterogeneous_network_results.csv')))
doc=SimpleDocTemplate(P+'heterogeneous_network_report.pdf',pagesize=(612,792),rightMargin=54,leftMargin=54,topMargin=48,bottomMargin=45)
s=getSampleStyleSheet();s.add(ParagraphStyle(name='Title2',parent=s['Title'],fontSize=16,leading=19,spaceAfter=15));s.add(ParagraphStyle(name='Body2',parent=s['BodyText'],fontSize=9.6,leading=14,spaceAfter=10));s.add(ParagraphStyle(name='Small2',parent=s['BodyText'],fontSize=8.6,leading=12,spaceAfter=7));s.add(ParagraphStyle(name='Head2',parent=s['Heading2'],fontSize=11,spaceBefore=11,spaceAfter=6))
story=[]
def p(t,sty='Body2'):story.append(Paragraph(t,s[sty]))
p('Persistent Frequency Heterogeneity and Network Resource Viability','Title2')
p('A spectral dissipation and resource-supply certificate for complex Stuart–Landau oscillators. November 2026.','Small2')
p('<b>Model and assumptions.</b> Let L be the symmetric graph Laplacian of a connected undirected network, with eigenvalues 0 = lambda_1 &lt; lambda_2 and lambda_max. Consider dz_i/dt = (a + i omega_i)z_i - |z_i|²z_i - k(Lz)_i, with common a, k &gt; 0 and complex z_i. Define P = I - 11ᵀ/N, y = ||Pz||_2, and J(t) = c_A k² integral_0^t ||Lz(s)||_2² ds. Resource evolves as E(t) = E_0 + (q-b)t - J(t).')
p('<b>Proposition 1 (bounded amplitude).</b> Set R = max(max_i |z_i(0)|, sqrt(max(a,0))). Then |z_i(t)| ≤ R for all i and t ≥ 0. At a site attaining maximal modulus, the real part of the diffusive coupling is nonpositive in the radial derivative. Hence the maximum modulus cannot cross R upward.')
p('<b>Proposition 2 (forced disagreement).</b> Choose omega_ref and let delta = max_i |omega_i - omega_ref|. Put gamma = k lambda_2 - a and B = delta sqrt(N)R. The cubic map z to |z|²z is monotone over real two-dimensional coordinates. Since Pz is orthogonal to constants, the common-frequency rotation makes zero contribution to d||Pz||²/dt, and the heterogeneous-frequency term has magnitude at most B y. Therefore D⁺y ≤ -gamma y + B. If gamma &gt; 0, y(t) ≤ r + (y_0-r) exp(-gamma t), r = B/gamma. This bound applies to the full complex network, without assuming phase locking.')
p('<b>Proposition 3 (persistent expenditure and viability).</b> The spectral inequality ||Lz|| ≤ lambda_max y gives J’(t) ≤ C y(t)², C = c_A k² lambda_max². For y_0 ≥ r and s = y_0-r ≥ 0, integration yields J(t) ≤ C[r²t + (2rs/gamma)(1-exp(-gamma t)) + (s²/(2gamma))(1-exp(-2gamma t))]. If q-b ≥ C r², then E(t) ≥ E_0 - C[2rs/gamma + s²/(2gamma)] for all t ≥ 0. If y_0 &lt; r, the simpler envelope y(t) ≤ r yields E(t) ≥ E_0 whenever q-b ≥ C r². These are sufficient, not necessary, conditions.')
p('Numerical certificate','Head2')
p('The four initial real states are (0.45, 0.21, 0.39, 0.27), with zero initial imaginary components. Frequencies are (-0.03, -0.01, 0.01, 0.03), a = 0.19, k = 0.7, c_A = 0.1, E_0 = 0.1, and q-b = 0.005. We use R = 0.45, delta = 0.03, B = 0.027, and exact rational spectral bounds: path lambda_2 ≥ 0.585, lambda_max ≤ 3.415; cycle (2,4); complete (4,4). The rational script rounds sqrt(y_0²) upward at 40 binary places.')
data=[['Graph','Rate upper','Transient buffer','All-time E floor']]
for row in rows:data.append([row['graph'],f"{float(row['steady_expenditure_rate_upper']):.8f}",f"{float(row['transient_expenditure_buffer']):.8f}",f"{float(row['all_time_energy_floor']):.8f}" if row['all_time_energy_floor'] else 'not certified'])
t=Table(data,colWidths=[112,105,119,160],repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e8edf2')),('GRID',(0,0),(-1,-1),.35,colors.grey),('FONTSIZE',(0,0),(-1,-1),8),('BOTTOMPADDING',(0,0),(-1,-1),8),('TOPPADDING',(0,0),(-1,-1),8)]));story.append(t)
p('<b>Interpretation.</b> The path certificate fails the supply-rate condition because its conservative persistent expenditure ceiling is 0.00864641 &gt; 0.005. This does not establish actual resource exhaustion. The cycle and complete-graph sufficient conditions certify E(t) &gt; 0.08 for every t ≥ 0. More connectivity improves the disagreement decay margin in these examples, but the lambda_max² actuator-cost factor can counteract that improvement in other graph families.')
p('<b>Scope and limitations.</b> The result assumes undirected diffusive coupling, identical real growth parameter, a global amplitude bound, no additional oscillator costs, and constant net resource supply. It is a comparison theorem, not a numerical trajectory simulation or a claim of biological calibration. The persistent expenditure bound is an upper bound, not proof of nonzero actual expenditure. The finite-horizon formula can still be used if net supply is smaller than the asymptotic expenditure ceiling.')
doc.build(story)
print(P+'heterogeneous_network_report.pdf')
