#!/usr/bin/env bash
set -euo pipefail

echo "== Add fitted trajectory output to SEFSC red snapper Quadra fit =="

python3 - <<'PY'
from pathlib import Path

p = Path("examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit.cpp")
s = p.read_text()

marker = "}  // namespace\n\nint main()"
if "write_fitted_trajectory" not in s:
    helper = r'''

void write_fitted_trajectory(
    const std::string& path,
    const std::vector<sefsc_red_snapper::Observation>& observations,
    const quadra::OptResult& fit) {
  if (fit.par.size() < 3) {
    throw std::runtime_error("Cannot write fitted trajectory: expected at least 3 fixed parameters");
  }

  sefsc_red_snapper::AgeStructuredParams params;
  params.log_r0 = fit.par[0];
  params.log_fbar = fit.par[1];
  params.log_q = fit.par[2];

  const auto rows =
      sefsc_red_snapper::run_deterministic_age_structured_model(observations,
                                                                params);

  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Could not open fitted trajectory CSV: " + path);
  }

  out << "year,recruitment,total_biomass,ssb_proxy,depletion,Fbar,"
      << "catch_obs,catch_hat,catch_log_residual,index_obs,index_hat,"
      << "index_log_residual\n";

  out << std::fixed << std::setprecision(6);

  for (const auto& row : rows) {
    const double catch_log_residual =
        std::log(std::max(row.catch_obs, 1.0e-12)) -
        std::log(std::max(row.catch_hat, 1.0e-12));
    const double index_log_residual =
        std::log(std::max(row.index_obs, 1.0e-12)) -
        std::log(std::max(row.index_hat, 1.0e-12));

    out << row.year << "," << row.recruitment << "," << row.total_biomass
        << "," << row.ssb_proxy << "," << row.depletion << ","
        << row.fbar << "," << row.catch_obs << "," << row.catch_hat
        << "," << catch_log_residual << "," << row.index_obs << ","
        << row.index_hat << "," << index_log_residual << "\n";
  }
}
'''
    if marker not in s:
        raise SystemExit("Could not find helper insertion marker")
    s = s.replace(marker, helper + "\n\n" + marker)

old = '''  const std::string summary_path =
      "examples/NMFS/sefsc_red_snapper/outputs/quadra_fit_summary.csv";
'''
new = '''  const std::string summary_path =
      "examples/NMFS/sefsc_red_snapper/outputs/quadra_fit_summary.csv";
  const std::string trajectory_path =
      "examples/NMFS/sefsc_red_snapper/outputs/quadra_fitted_trajectory.csv";
'''
if new not in s:
    if old not in s:
        raise SystemExit("Could not find summary_path block")
    s = s.replace(old, new)

old = '''  sefsc_red_snapper::write_fit_summary(summary_path, fit);

  std::cout << "SEFSC red-snapper-style Quadra fixed-effect fit\\n";
'''
new = '''  sefsc_red_snapper::write_fit_summary(summary_path, fit);
  sefsc_red_snapper::write_fitted_trajectory(trajectory_path, observations, fit);

  std::cout << "SEFSC red-snapper-style Quadra fixed-effect fit\\n";
'''
if new not in s:
    if old not in s:
        raise SystemExit("Could not find write_fit_summary call")
    s = s.replace(old, new)

old = '''  std::cout << "wrote:      " << summary_path << "\\n";
'''
new = '''  std::cout << "wrote:      " << summary_path << "\\n";
  std::cout << "wrote:      " << trajectory_path << "\\n";
'''
if new not in s:
    if old not in s:
        raise SystemExit("Could not find summary print")
    s = s.replace(old, new)

p.write_text(s)
PY

cat > examples/NMFS/sefsc_red_snapper/validation/fitted_trajectory_checklist.md <<'MD'
# Fitted Trajectory Checklist

- [x] fixed-effect fit summary written
- [x] fitted deterministic trajectory written
- [x] observed catch and predicted catch included
- [x] observed index and predicted index included
- [x] log residuals included
- [ ] residual diagnostics summary
- [ ] fitted trajectory plotted
- [ ] age-composition likelihood added
MD

echo
echo "Patched fitted trajectory output."
echo
echo "Run:"
echo "  ./examples/NMFS/sefsc_red_snapper/run_red_snapper_quadra_fit.sh"
echo "  head examples/NMFS/sefsc_red_snapper/outputs/quadra_fitted_trajectory.csv"
