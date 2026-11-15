"""Render the accepted native release's actual results as a two-page PDF summary."""
import csv
from html import escape
import json
from pathlib import Path
import sys
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate,Paragraph,Table,TableStyle,Spacer,PageBreak
from runtime import file_hash


def table(rows, widths):
    """Create a compact scientific results table with repeated, high-contrast headers."""
    result=Table(rows,colWidths=widths,repeatRows=1)
    result.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#173D5E')),
                               ('TEXTCOLOR',(0,0),(-1,0),colors.white),('FONTSIZE',(0,0),(-1,-1),8),
                               ('TOPPADDING',(0,0),(-1,-1),8),('BOTTOMPADDING',(0,0),(-1,-1),8),
                               ('LINEBELOW',(0,1),(-1,-1),.3,colors.lightgrey)]))
    return result


def build(bundle):
    """Report accepted baseline/local-start checks, uncertainty and exact input/source identity."""
    acceptance=json.loads((bundle/'acceptance.json').read_text());fit=acceptance['baseline']
    if not acceptance['passed']:raise ValueError('Cannot label a failed run as accepted')
    styles=getSampleStyleSheet();story=[Paragraph('Quadra: native swordfish release',styles['Title']),
        Paragraph('Reproducibility and local-start acceptance: PASS',styles['Heading2'])]
    for text in [
        'A clean C++ build ran from a frozen source snapshot using native normalized JSON. '
        'The baseline initial values match all 74 native parameters and 57 historical '
        'deviations preserved from Jon\'s SS3 .par file. No raw SS3 files or generated '
        'observation factories were included in the execution snapshot.',
        'The baseline fit calculated fresh joint-likelihood curvature and reference-point '
        'uncertainty. Three deterministic perturbations use 0.10, 0.25 and 0.50 local '
        'standard-deviation scales along the full covariance directions. They start '
        'from the original .par-derived vector plus perturbations. The baseline covariance '
        'also scales optimizer steps; independent convergence checks use the original coordinates.',
    ]:story.append(Paragraph(text,styles['BodyText']))
    rows=[['Run','Scale','Initial objective','Final objective','Gradient norm','Max SE distance'],
          ['Baseline','--',f"{fit['initial_objective']:.9f}",f"{fit['objective']:.9f}",f"{fit['gradient_norm']:.3g}",'--']]
    for r in acceptance['perturbed_starts']:
        rows.append([r['name'],str(r['scale']),f"{r['initial_objective']:.9f}",f"{r['objective']:.9f}",f"{r['gradient_norm']:.3g}",f"{r['maximum_coordinate_difference_in_SE']:.3g}"])
    story += [Spacer(1,12),table(rows,[65,40,108,108,85,90]),Spacer(1,12)]
    points=json.loads((bundle/'reference_point_comparison.json').read_text());ses=json.loads((bundle/'reference_SE_comparison.json').read_text())
    statements=[
        f"SS3 objective recorded with .par: {acceptance['ss3_objective']:.12f}. Quadra difference: {acceptance['ss3_objective_difference']:.3g} (absolute tolerance 1e-4).",
        f"All ten SS3 reference-point checks pass (relative tolerance 1e-5); maximum error {100*max(r['relative_error'] for r in points['checks']):.5g}%.",
        f"All ten SS3 reference-point SE checks pass (relative tolerance 1%); maximum difference {100*max(r['relative_difference'] for r in ses['checks']):.5g}%.",
        'Every perturbed run satisfies gradient norm at most 1e-4, objective recovery within 1e-5 and maximum coordinate difference below 0.01 local SE. These are local robustness checks, not proof of global optimality.',
    ]
    for text in statements:story.append(Paragraph(text,styles['BodyText']))
    story += [Paragraph('Reproduction and provenance',styles['Heading2']),Paragraph('Input SHA-256:',styles['BodyText']),Paragraph(file_hash(bundle/'inputs/model.json'),styles['Code'])]
    snapshot=json.loads((bundle/'source_snapshot.json').read_text())
    story.append(Paragraph(f"The bundle preserves {len(snapshot['files'])} source/target hashes, exact input JSON, perturbation seeds, compiler/build logs, executable hashes, fitted outputs and the two-step curvature/Jacobian checks. See release_manifest.json and acceptance.json. The bundle contains observations and is intended for deliberate sharing.",styles['BodyText']))
    story += [PageBreak(),Paragraph('Reference-point uncertainty',styles['Title'])]
    with (bundle/'results/baseline/management_reference_point_uncertainty.csv').open() as stream:uncertainty=list(csv.DictReader(stream))
    rows=[['Benchmark / quantity','Estimate','SE','CV','Lower 95%','Upper 95%']]
    names={'SSB_tonnes':'SSB (t)','dead_catch_tonnes':'Catch (t/year)','annual_F_SS3_units4_ages1_10':'Annual F'}
    for row in uncertainty:
        if row['quantity'] not in names or (row['benchmark']=='unfished' and row['quantity']!='SSB_tonnes'):continue
        fmt='.6f' if row['quantity'].startswith('annual_F') else ',.2f'
        rows.append([row['benchmark']+' / '+names[row['quantity']],format(float(row['estimate']),fmt),format(float(row['standard_error']),fmt),f"{100*float(row['CV']):.2f}%",format(float(row['lower95']),fmt),format(float(row['upper95']),fmt)])
    story += [table(rows,[150,75,65,50,78,78]),Spacer(1,14)]
    for text in [
        'SEs use the full correlated local delta method, J V J-transpose. Curvature is '
        'the inverse observed Hessian of the joint penalized likelihood in the 84 '
        'independent bounded-logit coordinates, including historical recruitment deviations.',
        'Both Hessian and reference-point Jacobians are checked at steps 1e-4 and 5e-5. '
        'Curvature is positive definite without regularization. Positive quantities use '
        'pointwise log-scale Wald intervals; nonconstant SPR proportions use logit-scale '
        'intervals. Fixed structural quantities have zero SE. All 24 reference quantities '
        'and their covariance are included in the CSV outputs.',
        'Annual F uses SS3 reporting units 4, pooled ages 1-10. Reference points use '
        'terminal-year biology and fleet/season fishing mixture. Fixed biology and '
        'hyperparameters remain fixed. These are local parameter-confidence intervals; '
        'future process variation, bootstrap/model uncertainty and projection validation '
        'are outside this release acceptance scope.',
    ]:story.append(Paragraph(text,styles['BodyText']))
    diagnostic=acceptance['uncertainty'];story.append(Paragraph(f"Maximum relative reference SE change across numerical steps: {diagnostic['maximum_reference_SE_step_change']:.3g}. Minimum Hessian eigenvalue: {diagnostic['minimum_eigenvalue']:.6g}.",styles['BodyText']))
    destination=bundle/'output/pdf/native_release_summary.pdf';destination.parent.mkdir(parents=True,exist_ok=True)
    def footer(canvas,document):
        """Identify the case and number pages without overlapping report content."""
        canvas.setFont('Helvetica',8);canvas.setFillColor(colors.HexColor('#56636B'))
        canvas.drawString(44,24,'QUADRA / NORTH PACIFIC SWORDFISH / NATIVE RELEASE ACCEPTANCE')
        canvas.drawRightString(letter[0]-44,24,str(document.page))
    SimpleDocTemplate(str(destination),pagesize=letter,leftMargin=38,rightMargin=38,topMargin=34,bottomMargin=42,title='Quadra native swordfish release acceptance').build(story,onFirstPage=footer,onLaterPages=footer)
    print('Release PDF: '+str(destination))


def main():
    """Build the PDF for one already numerically accepted staging bundle."""
    if len(sys.argv)!=2:raise SystemExit('Usage: build_native_release_report.py BUNDLE')
    build(Path(sys.argv[1]).resolve())


if __name__=='__main__':main()
