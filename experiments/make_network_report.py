from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,KeepTogether
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.units import inch
from pathlib import Path
P=Path('/mnt/data'); s=getSampleStyleSheet();s.add(ParagraphStyle(name='TitleCustom',parent=s['Title'],fontSize=16,leading=20,spaceAfter=15));s.add(ParagraphStyle(name='BodyCustom',parent=s['BodyText'],fontSize=10,leading=15,spaceAfter=10));s.add(ParagraphStyle(name='MathCustom',parent=s['BodyText'],fontName='Courier',fontSize=9,leading=14,leftIndent=15,spaceAfter=10))
doc=SimpleDocTemplate(str(P/'network_dissipation_report.pdf'),pagesize=(612,792),rightMargin=55,leftMargin=55,topMargin=54,bottomMargin=52)
story=[]
def p(t,style='BodyCustom'):story.append(Paragraph(t.replace('<','&lt;').replace('>','&gt;'),s[style]))
p('A Spectral Expenditure Certificate for Stuart–Landau Networks','TitleCustom')
p('Scope. Consider N real oscillators on an undirected, connected, weighted graph with symmetric positive-semidefinite Laplacian L, satisfying x_i\' = a x_i - x_i^3 - k(Lx)_i. Here a is a common fixed real growth parameter with a ≤ a_max, k > 0, and c_A > 0. Define J\' = c_A k² ||Lx||² and J(0)=0. The result is valid for arbitrary real initial states and all future times; global existence follows from the dissipative cubic term.')
p('Let P=I-11ᵀ/N, y=Px, and V=||y||²/2. Since L1=0 and the scalar cubic is monotone, the identity Σ_i (x_i−x̄)x_i³ = (1/N)Σ_{i<j}(x_i−x_j)(x_i³−x_j³) ≥ 0 gives')
p('V\' = a||y||² − Σ_i y_i x_i³ − k yᵀLy ≤ 2a_max V − k yᵀLy.','MathCustom')
p('The spectral Poincaré inequality yᵀLy ≥ 2λ₂V implies V(t) ≤ V(0) exp[−2(kλ₂−a_max)t] whenever kλ₂>a_max. Because L² ≤ λ_max L as quadratic forms, ||Lx||² ≤ λ_max yᵀLy. Integrating the preceding differential inequality yields')
p('k ∫₀ᵀ yᵀLy dt ≤ V(0) + 2a_max ∫₀ᵀ V dt ≤ V(0) kλ₂/(kλ₂−a_max),','MathCustom')
p('provided a_max ≥ 0. Therefore, for every T ≥ 0,')
p('J(T) ≤ c_A k² λ_max λ₂ V(0)/(kλ₂−a_max).','MathCustom')
p('The alternative direct exponential bound is J(T) ≤ c_A k² λ_max² V(0)/(kλ₂−a_max). The edge-energy estimate is no larger when λ₂ ≤ λ_max. The proof also works with upper and lower rational eigenvalue enclosures, replacing λ₂ by a certified lower bound and λ_max by a certified upper bound in the monotone expression c_A k² λ_max V(0)·λ₂/(kλ₂−a_max). Since this expression decreases with λ₂ for positive a_max, the replacement is conservative.')
p('Reproducible four-node example')
p('Use x(0)=(0.45, 0.21, 0.39, 0.27), V(0)=9/500=0.018, k=0.7, a_max=0.19, and c_A=0.1. The path eigenvalues are 0, 2−√2, 2, 2+√2; the script uses rational bounds λ₂≥0.585 and λ_max≤3.415. The cycle has eigenvalues 0,2,2,4, and the complete graph has 0,4,4,4.')
t=Table([['Graph','λ₂ lower','λmax upper','J upper','E floor lower'],['Path P4','0.585','3.415','0.00802751','0.09197249'],['Cycle C4','2','4','0.00583140','0.09416860'],['Complete K4','4','4','0.00540690','0.09459310']],colWidths=[92,82,88,108,110]);t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e7edf3')),('GRID',(0,0),(-1,-1),.35,colors.grey),('FONTSIZE',(0,0),(-1,-1),8.7),('BOTTOMPADDING',(0,0),(-1,-1),8),('TOPPADDING',(0,0),(-1,-1),8)]));story.append(t);story.append(Spacer(1,12))
p('Resource guarantee. For E\' = q−b−J\', with E(0)=0.1 and q−b≥0, E(t)≥0.1−J(∞). Thus all three four-node networks certify E(t)>0.08 for every t≥0. If the full model also contains additional nonnegative oscillator operating costs, those costs must be separately bounded; the displayed E guarantee covers only basal demand and the actuator expenditure.')
p('Limitations. The theorem assumes real states, a common growth parameter, undirected symmetric coupling, and no additive forcing or heterogeneous intrinsic frequencies. For complex oscillators with common frequency, the cubic monotonicity argument extends in real inner-product form, and the common rotational term contributes no disagreement work. Heterogeneous frequencies introduce an additional signed work term that must be bounded. No claim about biological calibration or uncontrolled network dynamics follows from this synthetic certificate.')
doc.build(story)
print(P/'network_dissipation_report.pdf')
