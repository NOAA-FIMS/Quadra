"""Create a concise PDF describing the saved control-input initialization experiment."""
import argparse
import csv
import json
from pathlib import Path
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.graphics.shapes import Drawing, Line, String, PolyLine

ROOT = Path(__file__).resolve().parents[1]


def build(evidence, output):
    """Render saved numerical evidence without running or changing either fit.

    Parameters are paths to a comparison directory and a destination PDF. The
    comparison JSON records acceptance checks; the PDF preserves failures as
    failures. Biomass lines use the independently saved annual comparison CSV.
    """
    result = json.loads((evidence / 'comparison_summary.json').read_text())
    styles = getSampleStyleSheet()
    story = [Paragraph('Swordfish: starting from SS3 inputs', styles['Title'])]
    story.append(Paragraph('Independent Quadra joint-likelihood initialization experiment', styles['Heading2']))
    text = [
        'The optimizer reads the original observation (.dat) and control (.ctl) inputs. '
        'All estimated native parameters start at their control initial values; all 57 '
        'historical recruitment deviations start at zero. Neither Report.sso estimates '
        'nor a saved covariance matrix supplies optimization values.',
        'The literal control point cannot evaluate the full likelihood: catch depletion '
        'produces an empty predicted length composition. A declared catch-only stage '
        'first adjusts recruitment scale to make the population feasible. Control '
        'activation phases then introduce 1, 8, 64 and 84 independent coordinates, '
        'using the final weighted joint objective. This does not replay SS3 phase-specific weights.',
        'Refinement learns a linear optimization metric from differences of Quadra AD '
        'gradients along this cold trajectory. It is not an inferential covariance. '
        'Numerically saturated bounded logits are moved inward before refinement. '
        'The recorded experiment resumes its own saved cold checkpoints; it never '
        'resumes the SS3-initialized reference fit.',
    ]
    for item in text:
        story.append(Paragraph(item, styles['BodyText']))
    story.append(Spacer(1, 12))
    rows = [['Configuration', 'Male L1', 'Objective', 'Full gradient']]
    for label, key in [('Strict controls', 'strict'), ('Matched fixed growth', 'matched_fixed')]:
        row = result[key]
        rows.append([label, f"{row['fixed_male_L1']:.1f}", f"{row['objective']:.9f}", f"{row['gradient_norm']:.3g}"])
    rows.append(['SS3-initialized Quadra', '83.4', f"{result['ss3_initialized']['objective']:.9f}", 'Reference'])
    table = Table(rows, colWidths=[195, 62, 112, 100])
    table.setStyle(TableStyle([('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#173d5e')),
                              ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
                              ('FONTSIZE', (0, 0), (-1, -1), 9),
                              ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
                              ('TOPPADDING', (0, 0), (-1, -1), 8),
                              ('LINEBELOW', (0, 1), (-1, -1), .3, colors.lightgrey)]))
    story.append(table)
    story.append(Paragraph('Comparison checks: ' + ('PASS' if result['passed'] else 'FAIL'), styles['Heading2']))
    story.append(Paragraph(
        f"Matched objective difference: {result['matched_objective_difference']:.3g}; "
        f"largest parameter difference: {result['maximum_parameter_difference_in_reference_SE']:.3g} reference SE; "
        f"largest annual spawning-biomass difference: {result['maximum_matched_annual_SSB_difference_tonnes']:.3g} tonnes. "
        'Strict controls use male L1=83.2, while the supplied fitted report uses 83.4; '
        'their objective values describe different fixed configurations.', styles['BodyText']))
    with (evidence / 'annual_ssb_comparison.csv').open() as stream:
        annual = [r for r in csv.DictReader(stream) if r['configuration'] == 'matched_fixed']
    drawing = Drawing(470, 180)
    years = [int(r['year']) for r in annual]
    values = [[float(r[k]) for r in annual] for k in ['ss3_initialized_ssb_tonnes', 'control_start_ssb_tonnes']]
    ymax = max(max(v) for v in values) * 1.05
    drawing.add(String(48, 165, 'Female spawning biomass (tonnes)', fontSize=10))
    drawing.add(Line(48, 30, 460, 30)); drawing.add(Line(48, 30, 48, 150))
    for fraction in [0, .5, 1]:
        y = 30 + 120 * fraction
        drawing.add(String(43, y - 3, f'{ymax*fraction:,.0f}', fontSize=8, textAnchor='end'))
    for year in [years[0], 2000, years[-1]]:
        x = 48 + 412 * (year - years[0]) / (years[-1] - years[0])
        drawing.add(String(x, 15, str(year), fontSize=8, textAnchor='middle'))
    for series, color, dash in zip(values, [colors.HexColor('#173d5e'), colors.HexColor('#d8662b')], [None, [4, 3]]):
        points = []
        for year, value in zip(years, series):
            points.extend([48 + 412 * (year-years[0])/(years[-1]-years[0]), 30 + 120*value/ymax])
        drawing.add(PolyLine(points, strokeColor=color, strokeWidth=1.5, strokeDashArray=dash))
    story.append(drawing)
    story.append(Paragraph('Blue: SS3-initialized Quadra. Orange dashed: matched-fixed control start. '
                           'Overlapping lines indicate agreement. This experiment tests initialization; '
                           'it does not validate missing SS3 projection/bootstrap files or the separate striped marlin model.', styles['BodyText']))
    output.parent.mkdir(parents=True, exist_ok=True)
    SimpleDocTemplate(str(output), title='Quadra swordfish control-input starts',
                      leftMargin=44, rightMargin=44, topMargin=36, bottomMargin=36).build(story)


def main():
    """Accept explicit saved comparison evidence and an optional PDF destination."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=ROOT/'output/pdf/swordfish_control_start.pdf')
    args = parser.parse_args()
    build(args.evidence, args.output)


if __name__ == '__main__':
    main()
