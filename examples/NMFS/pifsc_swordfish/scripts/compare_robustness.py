"""Require convergence to the same local optimum from every configured perturbation."""
import csv
import json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

def main():
    """Execute the compare robustness command, write its documented outputs and exit nonzero if validation fails.

    Notes
    -----
    The three seeded perturbations test a local neighborhood of the fitted solution. Convergence to the same objective is not a proof of global optimality or recovery from unfitted control initial values.
    """
    rows = list(csv.DictReader((ROOT/'build/robustness.csv').open()))
    baseline = json.loads((ROOT/'build/quadra_fit_summary.json').read_text())['objective']
    passed = len(rows) == 3 and all(r['converged'] == '1' and float(r['gradient_norm']) < 1e-4 and abs(float(r['final_objective'])-baseline) < 1e-6 for r in rows)
    result = dict(passed=passed, starts=len(rows), baseline_objective=baseline, maximum_objective_difference=max(abs(float(r['final_objective'])-baseline) for r in rows), coordinate_rescaling='constant Cholesky factor of local covariance; original-coordinate objective and gradient chain rule verified', scope='local perturbed starts; no global optimum claim')
    (ROOT/'build/robustness_summary.json').write_text(json.dumps(result, indent=2)+'\n')
    if not passed:
        raise SystemExit('Perturbed-start check failed')
    print('All three perturbed starts converged to the same objective.')
if __name__ == '__main__':
    main()
