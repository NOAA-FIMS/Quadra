"""Check point estimates, intervals and full covariance propagation; save input provenance."""
import csv
import hashlib
import json
import math
from pathlib import Path
from runtime import verify_covariance
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'build/projections'

def read(path):
    """Read a CSV with column labels into dictionaries for independent uncertainty artifact validation.

    Parameters
    ----------
    path
        Filesystem input/output path in the schema described by the function summary.

    Notes
    -----
    Full covariance propagation is checked independently for selected cross-output entries. Bounds are pointwise local log-scale Wald intervals. They exclude future recruitment process variation and cannot be interpreted as simultaneous trajectory coverage or bootstrap quantiles.
    """
    with path.open() as f:
        return list(csv.DictReader(f))

def main():
    """Execute the check projection uncertainty command, write its documented outputs and exit nonzero if validation fails.

    Notes
    -----
    Full covariance propagation is checked independently for selected cross-output entries. Bounds are pointwise local log-scale Wald intervals. They exclude future recruitment process variation and cannot be interpreted as simultaneous trajectory coverage or bootstrap quantiles.
    """
    verify_covariance(ROOT)
    for source in json.loads((ROOT/'reference/source_manifest.json').read_text()):
        raw = ROOT/'data/raw'/source['file']
        if hashlib.sha256(raw.read_bytes()).hexdigest() != source['sha256']:
            raise ValueError('Source input changed since the reference covariance was prepared')
    estimates = {(int(r['scenario']), int(r['year'])): r for r in read(OUT/'deterministic_projections.csv')}
    rows = read(OUT/'projection_uncertainty.csv')
    quantities = {'spawning_biomass_tonnes', 'recruitment_thousands', 'biomass_january_tonnes', 'dead_catch_tonnes', 'annual_F_SS3_units4_ages1_10'}
    keys = {(int(r['scenario']), int(r['year']), r['quantity']) for r in rows}
    expected_keys = {(s, y, q) for s in range(1, 6) for y in range(2022, 2032) for q in quantities}
    if keys != expected_keys or len(rows) != 250:
        raise ValueError('Expected five quantities for every scenario and year')
    maximum_point_difference = 0.
    for r in rows:
        estimate, se, low, high, cv = [float(r[k]) for k in ['estimate', 'standard_error', 'lower95', 'upper95', 'coefficient_of_variation']]
        reference = float(estimates[int(r['scenario']), int(r['year'])][r['quantity']])
        maximum_point_difference = max(maximum_point_difference, abs(estimate/reference-1))
        if not all(map(math.isfinite, [estimate, se, low, high, cv])) or not 0 < low < estimate < high or se <= 0:
            raise ValueError('Invalid uncertainty interval')
        if abs(math.sqrt(low*high)/estimate-1) > 1e-12 or abs(se/estimate-cv) > 1e-12:
            raise ValueError('Positive interval or coefficient of variation inconsistent')
    if maximum_point_difference > 1e-10:
        raise ValueError('Uncertainty point estimates differ from deterministic forecasts')
    with (OUT/'projection_uncertainty_2031.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0]);writer.writeheader()
        writer.writerows(r for r in rows if r['year'] == '2031')
    labels = [f"S{r['scenario']}_Y{r['year']}_{r['quantity']}" for r in rows]
    covariance_rows = read(OUT/'projection_output_covariance.csv')
    jacobian_rows = read(OUT/'projection_jacobian.csv')
    parameter_covariance = read(ROOT/'build/coordinate_covariance.csv')
    if [r['output'] for r in covariance_rows] != labels or [r['output'] for r in jacobian_rows] != labels:
        raise ValueError('Output covariance or Jacobian order differs')
    parameters = [r['parameter'] for r in parameter_covariance]
    if len(parameters) != 84 or list(jacobian_rows[0])[1:] != parameters or list(covariance_rows[0])[1:] != labels:
        raise ValueError('Covariance column order differs')
    for i, row in enumerate(covariance_rows):
        if abs(float(row[labels[i]])/float(rows[i]['standard_error'])**2-1) > 1e-12:
            raise ValueError('Variance does not equal squared reported SE')
        for j in range(i):
            a, b = float(row[labels[j]]), float(covariance_rows[j][labels[i]])
            if abs(a-b) > 1e-10*max(1., abs(a), abs(b)):
                raise ValueError('Asymmetric output covariance')
    # Independent scalar recomputation checks selected cross-year and cross-quantity covariances.
    for i, j in [(0, 1), (0, 45), (49, 249), (102, 187), (0, 0), (249, 249)]:
        g, h = jacobian_rows[i], jacobian_rows[j]
        expected = sum(float(g[p])*sum(float(parameter_covariance[k][q])*float(h[q]) for q in parameters) for k, p in enumerate(parameters))
        actual = float(covariance_rows[i][labels[j]])
        if abs(actual-expected) > 1e-9*max(abs(expected), 1.):
            raise ValueError('Full covariance propagation inconsistent')
    numerical = json.loads((OUT/'projection_uncertainty_numerics.json').read_text())
    if not numerical['passed']:
        raise ValueError('Numerical validation failed')
    paths = ['build/quadra_fitted_parameters.csv', 'build/coordinate_covariance.csv', 'build/uncertainty_coordinates.csv', 'data/projection_scenarios.csv', 'reference/source_manifest.json']
    result = dict(passed=True, intervals=250, point_estimate_maximum_relative_difference=maximum_point_difference,
                  interval='pointwise 95% log-scale Wald; not simultaneous bands or future outcome prediction intervals',
                  parameter_uncertainty='full 84-coordinate local joint-likelihood covariance',
                  target_rule='recalibrate each scenario and 2019-2021 fishing mixture under each parameter perturbation',
                  omitted=['future recruitment process variation', 'model structural uncertainty', 'original assessment bootstrap resampling'],
                  input_sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}, numerical=numerical)
    (OUT/'projection_uncertainty_validation.json').write_text(json.dumps(result, indent=2)+'\n')
    print(f"All 250 pointwise projection intervals pass; point estimates agree within {maximum_point_difference:.3g} relative.")
if __name__ == '__main__':
    main()
