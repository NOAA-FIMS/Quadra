#!/usr/bin/env bash
set -euo pipefail

F="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/quadra/bigeye_level21_age_based_natural_mortality_diagnostic.cpp"
cp "$F" "$F.before_descent_probe_compile_fix.$(date +%Y%m%d_%H%M%S)"

python3 - <<'PY'
from pathlib import Path

p = Path("examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/quadra/bigeye_level21_age_based_natural_mortality_diagnostic.cpp")
s = p.read_text()

start = s.find('    {\n      std::cout << "\\nBigeye Level 21 Optimizer Descent Probe\\n";')
end = s.find('pifsc_bigeye_tuna::write_bigeye_report_suite(report_paths, observations,', start)

if start < 0 or end < 0:
    raise SystemExit("Could not find probe block bounds")

replacement = r'''    {
      std::cout << "\nBigeye Level 21 Optimizer Descent Probe\n";
      std::cout << "======================================\n";
      std::cout << "alpha,fx_minus_alpha_g,delta_minus,fx_plus_alpha_g,delta_plus,max_abs_step\n";

      Eigen::VectorXd x(static_cast<Eigen::Index>(fit.par.size()));
      Eigen::VectorXd g(static_cast<Eigen::Index>(fit.fixed_gradient.size()));

      for (std::size_t i = 0; i < fit.par.size(); ++i) {
        x[static_cast<Eigen::Index>(i)] = fit.par[i];
      }
      for (std::size_t i = 0; i < fit.fixed_gradient.size(); ++i) {
        g[static_cast<Eigen::Index>(i)] = fit.fixed_gradient[i];
      }

      const double fx0 = fit.value;

      std::vector<int> fixed_idx;
      std::vector<int> random_idx;
      fixed_idx.reserve(fit.par.size());
      random_idx.reserve(objective.n_years());

      for (std::size_t i = 0; i < fit.par.size(); ++i) {
        fixed_idx.push_back(static_cast<int>(i));
      }
      for (std::size_t t = 0; t < objective.n_years(); ++t) {
        random_idx.push_back(static_cast<int>(fit.par.size() + t));
      }

      quadra::LBFGSObjective<pifsc_bigeye_tuna::BigeyeQuadraObjective>
          eval_obj(objective, params, fixed_idx, random_idx, opts);

      const std::array<double, 7> alphas = {
          1.0e-2, 3.0e-3, 1.0e-3, 3.0e-4, 1.0e-4, 3.0e-5, 1.0e-5};

      for (const double alpha : alphas) {
        Eigen::VectorXd x_minus = x - alpha * g;
        Eigen::VectorXd x_plus = x + alpha * g;

        Eigen::VectorXd g_tmp_minus;
        Eigen::VectorXd g_tmp_plus;

        const double f_minus = eval_obj(x_minus, g_tmp_minus);
        const double f_plus = eval_obj(x_plus, g_tmp_plus);
        const double max_abs_step = (alpha * g).cwiseAbs().maxCoeff();

        std::cout << std::setprecision(12)
                  << alpha << ","
                  << f_minus << ","
                  << (f_minus - fx0) << ","
                  << f_plus << ","
                  << (f_plus - fx0) << ","
                  << max_abs_step << "\n";
      }
    }

'''
p.write_text(s[:start] + replacement + s[end:])
PY

echo "patched: $F"
echo "run: ./run_bigeye_level21_age_based_m_check.sh"
