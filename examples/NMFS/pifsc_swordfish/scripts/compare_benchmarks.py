"""Compare independent equilibrium calculations with rounded SS3 report targets."""
import csv
import json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

def compare(points_path, report):
    """Compare a model-generated reference-point CSV with ten printed SS3 targets.

    Expected values are supplied only after model calculation; no SS3 values
    enter equilibrium calculations. The 1e-5 relative tolerance accommodates
    rounded report targets. Missing reference rows or fields fail explicitly.
    """
    rows = {r['benchmark']: r for r in csv.DictReader(points_path.open())}
    checks = []
    for benchmark, suffix in [('unfished', 'unfished'), ('Btgt20', 'Btgt'), ('SPR20', 'SPR'), ('MSY', 'MSY')]:
        fields = [('SSB_tonnes', 'SSB_')]
        if benchmark != 'unfished':
            fields += [('dead_catch_tonnes', 'Dead_Catch_'), ('annual_F_SS3_units4_ages1_10', 'annF_')]
        for field, prefix in fields:
            actual, expected = float(rows[benchmark][field]), report[prefix + suffix]
            error = abs(actual - expected) / max(abs(expected), 1e-12)
            checks.append(dict(quantity=prefix + suffix, quadra=actual, ss3=expected, relative_error=error, passed=error < 1e-5))
    return dict(passed=all(r['passed'] for r in checks), tolerance_relative=1e-5, checks=checks)

def main():
    """Execute the compare benchmarks command, write its documented outputs and exit nonzero if validation fails.

    Notes
    -----
    Reference points are calculated independently before comparison. Terminal-2021 and 2019-2021 mean fishing patterns define different benchmark calculations and must not be silently interchanged.
    """
    report = {}
    for line in (ROOT / 'data/raw/Report.sso').read_text().splitlines():
        cells = line.split()
        if len(cells) == 3:
            try:
                report[cells[0]] = float(cells[1])
            except ValueError:
                pass
    result = compare(ROOT/'build/management_reference_points.csv', report)
    (ROOT / 'build/benchmark_comparison.json').write_text(json.dumps(result, indent=2) + '\n')
    if not result['passed']:
        raise SystemExit('Management benchmark comparison failed')
    print(f"All {len(result['checks'])} management benchmark comparisons pass.")
if __name__ == '__main__':
    main()
