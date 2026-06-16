#!/usr/bin/env bash
set -euo pipefail

python3 - <<'PY'
from pathlib import Path

p = Path("examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit.cpp")
s = p.read_text()

if "ResidualDiagnostics" not in s:
    helper = r'''
struct ResidualDiagnostics {
  int n = 0;
  double catch_rmse_log = 0.0;
  double index_rmse_log = 0.0;
  double catch_mean_log_residual = 0.0;
  double index_mean_log_residual = 0.0;
  double max_abs_catch_log_residual = 0.0;
  double max_abs_index_log_residual = 0.0;
};

void write_residual_diagnostics(
    const std::string& path,
    const std::vector<sefsc_red_snapper::Observation>& observations,
    const quadra::OptResult& fit) {
  sefsc_red_snapper::AgeStructuredParams params;
  params.log_r0 = fit.par[0];
  params.log_fbar = fit.par[1];
  params.log_q = fit.par[2];

  const auto rows =
      sefsc_red_snapper::run_deterministic_age_structured_model(observations,
                                                                params);

  ResidualDiagnostics d;
  d.n = static_cast<int>(rows.size());

  double catch_sum = 0.0, catch_ss = 0.0;
  double index_sum = 0.0, index_ss = 0.0;

  for (const auto& row : rows) {
    const double cr = std::log(std::max(row.catch_obs, 1.0e-12)) -
                      std::log(std::max(row.catch_hat, 1.0e-12));
    const double ir = std::log(std::max(row.index_obs, 1.0e-12)) -
                      std::log(std::max(row.index_hat, 1.0e-12));

    catch_sum += cr;
    catch_ss += cr * cr;
    index_sum += ir;
    index_ss += ir * ir;

    d.max_abs_catch_log_residual =
        std::max(d.max_abs_catch_log_residual, std::abs(cr));
    d.max_abs_index_log_residual =
        std::max(d.max_abs_index_log_residual, std::abs(ir));
  }

  if (d.n > 0) {
    d.catch_mean_log_residual = catch_sum / d.n;
    d.index_mean_log_residual = index_sum / d.n;
    d.catch_rmse_log = std::sqrt(catch_ss / d.n);
    d.index_rmse_log = std::sqrt(index_ss / d.n);
  }

  std::ofstream out(path);
  out << "metric,value,note\n";
  out << std::setprecision(12);
  out << "n," << d.n << ",number of fitted years\n";
  out << "catch_rmse_log," << d.catch_rmse_log << ",root mean squared log catch residual\n";
  out << "index_rmse_log," << d.index_rmse_log << ",root mean squared log index residual\n";
  out << "catch_mean_log_residual," << d.catch_mean_log_residual << ",mean log observed minus predicted catch\n";
  out << "index_mean_log_residual," << d.index_mean_log_residual << ",mean log observed minus predicted index\n";
  out << "max_abs_catch_log_residual," << d.max_abs_catch_log_residual << ",maximum absolute log catch residual\n";
  out << "max_abs_index_log_residual," << d.max_abs_index_log_residual << ",maximum absolute log index residual\n";
}
'''
    s = s.replace("\nint main()", "\n" + helper + "\nint main()", 1)

s = s.replace(
    '  const std::string trajectory_path =\n'
    '      "examples/NMFS/sefsc_red_snapper/outputs/quadra_fitted_trajectory.csv";',
    '  const std::string trajectory_path =\n'
    '      "examples/NMFS/sefsc_red_snapper/outputs/quadra_fitted_trajectory.csv";\n'
    '  const std::string residual_diagnostics_path =\n'
    '      "examples/NMFS/sefsc_red_snapper/outputs/quadra_fit_residual_diagnostics.csv";',
    1,
)

s = s.replace(
    "  write_fitted_trajectory(trajectory_path, observations, fit);\n",
    "  write_fitted_trajectory(trajectory_path, observations, fit);\n"
    "  write_residual_diagnostics(residual_diagnostics_path, observations, fit);\n",
    1,
)

s = s.replace(
    '  std::cout << "wrote:      " << trajectory_path << "\\n";\n',
    '  std::cout << "wrote:      " << trajectory_path << "\\n";\n'
    '  std::cout << "wrote:      " << residual_diagnostics_path << "\\n";\n',
    1,
)

p.write_text(s)
PY
