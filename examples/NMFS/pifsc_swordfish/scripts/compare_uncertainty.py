"""Audit local joint-likelihood standard errors against SS3 reported values."""
import csv
import json
from pathlib import Path
from import_ss3 import read_lines, extract
ROOT = Path(__file__).resolve().parents[1]

def main():
    """Execute the compare uncertainty command, write its documented outputs and exit nonzero if validation fails.

    Notes
    -----
    The covariance is a local joint-likelihood approximation conditional on fixed biology/hyperparameters. Two finite-difference steps test numerical stability; comparison to rounded report SEs does not establish repeated-sampling interval coverage.
    """
    _, rows = extract(read_lines(ROOT/'data/raw/Report.sso'), 'PARAMETERS', 'Num Label Value', r'^\d+\s')
    coarse = list(csv.DictReader((ROOT/'build/derived_uncertainty_step1e-4.csv').open()))
    fine = list(csv.DictReader((ROOT/'build/derived_uncertainty.csv').open()))
    if [r['quantity'] for r in coarse] != [r['quantity'] for r in fine]:
        raise ValueError('Derived uncertainty labels differ across steps')
    change = max(abs(float(a['standard_error'])/float(b['standard_error'])-1) for a, b in zip(coarse, fine))
    (ROOT/'build/derived_step_check.json').write_text(json.dumps(dict(steps=[1e-4, 5e-5], maximum_relative_SE_change=change, passed=change < 1e-4), indent=2)+'\n')
    if change >= 1e-4:
        raise ValueError('Derived standard errors are sensitive to differentiation step')
    report = {r[1]: float(r[10]) for r in rows if r[10] != '_'}
    mapping = {r['canonical_name']: r['report_label'] for r in csv.DictReader((ROOT/'build/biology_parameter_mapping.csv').open())}
    comparisons = []
    for row in csv.DictReader((ROOT/'build/parameter_uncertainty.csv').open()):
        label = mapping.get(row['parameter'])
        if label in report:
            actual, expected = float(row['standard_error']), report[label]
            comparisons.append(dict(quantity=row['parameter'], quadra_se=actual, ss3_se=expected, relative_difference=abs(actual/expected-1)))
    derived_report = {}
    for line in read_lines(ROOT/'data/raw/Report.sso'):
        fields = line.split()
        if len(fields) == 3:
            try:
                derived_report[fields[0]] = float(fields[2])
            except ValueError:
                pass
    for row in csv.DictReader((ROOT/'build/derived_uncertainty.csv').open()):
        expected = derived_report[row['quantity']]
        actual = float(row['standard_error'])
        comparisons.append(dict(quantity=row['quantity'], quadra_se=actual, ss3_se=expected, relative_difference=abs(actual/expected-1)))
    result = dict(method='local joint likelihood, fixed biology and hyperparameters; delta method', tolerance_relative=.01, comparisons=comparisons)
    result['passed'] = len(comparisons) == 39 and all(r['relative_difference'] < .01 for r in comparisons)
    (ROOT/'build/uncertainty_comparison.json').write_text(json.dumps(result, indent=2)+'\n')
    print(f"Compared {len(comparisons)} standard errors; maximum relative difference {max(r['relative_difference'] for r in comparisons):.6g}")
    if not result['passed']:
        raise SystemExit('Uncertainty comparison failed')
if __name__ == '__main__':
    main()
