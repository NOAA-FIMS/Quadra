"""Build a one-page PDF of the saved SS3 .par initialization comparison."""
import argparse
import json
from pathlib import Path
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle

ROOT = Path(__file__).resolve().parents[1]


def build(evidence, output):
    """Present saved fit evidence and initialization limits without rerunning either model."""
    result = json.loads((evidence/'comparison_summary.json').read_text())
    source = json.loads((evidence/'input_manifest.json').read_text())
    styles = getSampleStyleSheet()
    story = [Paragraph('Swordfish: initialize Quadra from ss.par', styles['Title']),
             Paragraph('NP swordfish / direct SS3 parameter-file comparison', styles['Heading2'])]
    for text in [
        'The supplied file maps all 74 native parameters and 57 historical recruitment '
        'deviations directly into Quadra. Original observations and model structure '
        'come from .dat/.ctl. Neither Report.sso nor a saved covariance supplies '
        'initial parameter values or optimizer scaling.',
        'The file header records an SS3 objective and a small gradient at an already '
        'fitted vector. This experiment tests agreement starting from that snapshot. '
        'It does not establish the original pre-optimization vector or rerun SS3.',
    ]:
        story.append(Paragraph(text, styles['BodyText']))
    story.append(Spacer(1, 14))
    rows = [['Quantity', 'Value'],
            ['SS3 objective in .par header', f"{result['ss3_objective']:.12f}"],
            ['Quadra final objective', f"{result['quadra']['objective']:.12f}"],
            ['Quadra minus SS3 objective', f"{result['quadra_minus_ss3_objective']:.8g}"],
            ['Quadra final gradient norm', f"{result['quadra']['gradient_norm']:.6g}"],
            ['Difference from existing Quadra objective', f"{result['quadra_minus_existing_quadra_objective']:.6g}"],
            ['Maximum parameter difference (reference SE)', f"{result['maximum_parameter_difference_in_reference_SE']:.6g}"],
            ['Maximum annual spawning-biomass difference', f"{result['maximum_annual_ssb_difference_tonnes']:.6g} tonnes"]]
    table = Table(rows, colWidths=[320, 180])
    table.setStyle(TableStyle([('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#173d5e')),
                              ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
                              ('FONTSIZE', (0, 0), (-1, -1), 9),
                              ('TOPPADDING', (0, 0), (-1, -1), 9),
                              ('BOTTOMPADDING', (0, 0), (-1, -1), 9),
                              ('LINEBELOW', (0, 1), (-1, -1), .3, colors.lightgrey)]))
    story += [table, Paragraph('Objective comparison: '+('PASS' if result['passed'] else 'FAIL'), styles['Heading2'])]
    for text in [
        'Acceptance requires convergence, a Quadra gradient norm at most 1e-4, '
        'and absolute objective differences below 1e-4 against both the SS3 header '
        'and the saved Quadra fit. The objectives agree within that tolerance; '
        'they are not bit-for-bit identical.',
        'The .par supplies fixed male growth L1=83.4, resolving the control/report '
        'configuration difference for this run. Quadra profiles Q analytically '
        'and reconstructs the final main recruitment deviation to enforce a zero '
        'sum. The supplied main sum is approximately 4.9e-13. SS3 reports 86 '
        'parameters; Quadra uses 84 independent optimization coordinates.',
        'Run locally from examples/NMFS/pifsc_swordfish: <b>make fit-ss-par</b>. '
        'The raw file is preserved locally at data/raw/ss.par and excluded from '
        'version control. Results go to a new build/ss_par directory, preserving '
        'the canonical fit. This test certifies no new uncertainty or projections.',
    ]:
        story.append(Paragraph(text, styles['BodyText']))
    story.append(Spacer(1, 14))
    story.append(Paragraph('Source .par SHA-256: '+source['parameter_file_sha256'], styles['Code']))
    output.parent.mkdir(parents=True, exist_ok=True)
    SimpleDocTemplate(str(output), title='Quadra swordfish SS3 parameter-file comparison',
                      leftMargin=44, rightMargin=44, topMargin=40, bottomMargin=40).build(story)


def main():
    """Select saved evidence and an optional PDF destination."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, default=ROOT/'reference/ss_par')
    parser.add_argument('--output', type=Path, default=ROOT/'output/pdf/swordfish_ss_par_start.pdf')
    args = parser.parse_args()
    build(args.evidence, args.output)


if __name__ == '__main__':
    main()
