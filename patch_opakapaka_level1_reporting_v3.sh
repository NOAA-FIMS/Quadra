#!/usr/bin/env bash
set -euo pipefail

echo "== Anchor-flexible patch: Opakapaka Level-1 reporting outputs =="

stamp="$(date +%Y%m%d_%H%M%S)"
mkdir -p .quadra_patch_backups
cpp="examples/opakapaka_projection/opakapaka_projection.cpp"
[[ -f "$cpp" ]] || { echo "ERROR: missing $cpp" >&2; exit 1; }
cp "$cpp" ".quadra_patch_backups/opakapaka_projection.cpp.level1_reporting_v3_${stamp}.bak"

python3 - <<'PY'
from pathlib import Path
import re

p = Path("examples/opakapaka_projection/opakapaka_projection.cpp")
s = p.read_text()

if "QUADRA_LEVEL1_UNCERTAINTY_REPORTING_V3" in s:
    print("Level-1 reporting already installed.")
    raise SystemExit(0)

for inc in ["#include <algorithm>\n", "#include <limits>\n"]:
    if inc not in s and "#include <fstream>\n" in s:
        s = s.replace("#include <fstream>\n", "#include <fstream>\n" + inc, 1)

helper = r'''
// QUADRA_LEVEL1_UNCERTAINTY_REPORTING_V3
struct LogQUncertaintyReport {
  double objective = std::numeric_limits<double>::quiet_NaN();
  double fd_step = std::numeric_limits<double>::quiet_NaN();
  double fd_gradient = std::numeric_limits<double>::quiet_NaN();
  double fd_hessian = std::numeric_limits<double>::quiet_NaN();
  double covariance_log_q = std::numeric_limits<double>::quiet_NaN();
  double se_log_q = std::numeric_limits<double>::quiet_NaN();
  double log_q = std::numeric_limits<double>::quiet_NaN();
  double q = std::numeric_limits<double>::quiet_NaN();
  double se_q = std::numeric_limits<double>::quiet_NaN();
  double log_q_lwr_95 = std::numeric_limits<double>::quiet_NaN();
  double log_q_upr_95 = std::numeric_limits<double>::quiet_NaN();
  double q_lwr_95 = std::numeric_limits<double>::quiet_NaN();
  double q_upr_95 = std::numeric_limits<double>::quiet_NaN();
};

template <class Model>
LogQUncertaintyReport compute_log_q_uncertainty_report(
    Model& model,
    quadra::ParameterVector& params,
    quadra::LaplaceOptions& opts,
    const quadra::OptResult& fit) {
  LogQUncertaintyReport out;
  if (fit.par.size() != 1) return out;

  const std::vector<int> fixed_idx = {0};
  std::vector<int> random_idx;
  for (std::size_t i = 1; i < params.size(); ++i) {
    random_idx.push_back(static_cast<int>(i));
  }

  auto eval_at = [&](double theta) {
    auto tmp = params;
    tmp.params.at(0).value = theta;
    Eigen::VectorXd x(1);
    x[0] = theta;
    had::ADGraph graph;
    auto u_hat = quadra::solve_random_effects_laplace(
        model, tmp, x, fixed_idx, random_idx, graph);
    auto res = quadra::laplace_eval_at_u_star(
        model, tmp, fixed_idx, random_idx, x, u_hat, graph, opts);
    return res.value;
  };

  out.objective = fit.value;
  out.log_q = fit.par.at(0);
  out.q = std::exp(out.log_q);
  out.fd_step = std::max(1.0e-5, 1.0e-4 * (1.0 + std::abs(out.log_q)));

  const double fm = eval_at(out.log_q - out.fd_step);
  const double fp = eval_at(out.log_q + out.fd_step);
  if (!std::isfinite(fm) || !std::isfinite(fp) || !std::isfinite(out.objective)) return out;

  out.fd_gradient = (fp - fm) / (2.0 * out.fd_step);
  out.fd_hessian = (fp - 2.0 * out.objective + fm) / (out.fd_step * out.fd_step);

  if (std::isfinite(out.fd_hessian) && out.fd_hessian > 0.0) {
    out.covariance_log_q = 1.0 / out.fd_hessian;
    out.se_log_q = std::sqrt(out.covariance_log_q);
    out.se_q = out.q * out.se_log_q;
    out.log_q_lwr_95 = out.log_q - 1.96 * out.se_log_q;
    out.log_q_upr_95 = out.log_q + 1.96 * out.se_log_q;
    out.q_lwr_95 = std::exp(out.log_q_lwr_95);
    out.q_upr_95 = std::exp(out.log_q_upr_95);
  }
  return out;
}

inline void write_uncertainty_summary_csv(const std::string& path, const LogQUncertaintyReport& u) {
  std::ofstream out(path);
  out << "field,value\n";
  out << "objective," << u.objective << "\n";
  out << "fd_step," << u.fd_step << "\n";
  out << "fd_gradient_log_q," << u.fd_gradient << "\n";
  out << "fd_hessian_log_q," << u.fd_hessian << "\n";
  out << "covariance_log_q," << u.covariance_log_q << "\n";
  out << "se_log_q," << u.se_log_q << "\n";
  out << "se_q," << u.se_q << "\n";
  out << "hessian_positive," << (u.fd_hessian > 0.0 ? "yes" : "no") << "\n";
}

inline void write_covariance_matrix_csv(const std::string& path, const LogQUncertaintyReport& u) {
  std::ofstream out(path);
  out << "row,col,value\n";
  out << "log_q,log_q," << u.covariance_log_q << "\n";
}

inline void write_correlation_matrix_csv(const std::string& path) {
  std::ofstream out(path);
  out << "row,col,value\n";
  out << "log_q,log_q,1\n";
}

inline void write_standard_errors_csv(const std::string& path, const LogQUncertaintyReport& u) {
  std::ofstream out(path);
  out << "parameter,scale,estimate,se\n";
  out << "log_q,log," << u.log_q << "," << u.se_log_q << "\n";
  out << "q,natural," << u.q << "," << u.se_q << "\n";
}

inline void write_confidence_intervals_csv(const std::string& path, const LogQUncertaintyReport& u) {
  std::ofstream out(path);
  out << "parameter,scale,estimate,se,lwr_95,upr_95\n";
  out << "log_q,log," << u.log_q << "," << u.se_log_q << "," << u.log_q_lwr_95 << "," << u.log_q_upr_95 << "\n";
  out << "q,natural," << u.q << "," << u.se_q << "," << u.q_lwr_95 << "," << u.q_upr_95 << "\n";
}

inline void write_random_effect_uncertainty_csv(const std::string& path, const std::vector<double>& u_hat) {
  std::ofstream out(path);
  out << "effect,mode,conditional_se,conditional_variance,note\n";
  for (std::size_t i = 0; i < u_hat.size(); ++i) {
    out << "log_B[" << i << "]," << u_hat[i] << ",,,pending selected-inverse/random-effect covariance extraction\n";
  }
}

inline void write_derived_quantities_csv(const std::string& path, const std::vector<Observation>& data, const std::vector<double>& u_hat, double q_hat) {
  std::ofstream out(path);
  out << "year,biomass,index_hat,depletion,F_proxy\n";
  const double b0 = u_hat.empty() ? std::numeric_limits<double>::quiet_NaN() : std::exp(u_hat.front());
  for (std::size_t i = 0; i < data.size() && i < u_hat.size(); ++i) {
    const double biomass = std::exp(u_hat[i]);
    const double depletion = b0 > 0.0 ? biomass / b0 : std::numeric_limits<double>::quiet_NaN();
    const double f_proxy = biomass > 0.0 ? data[i].catch_mt / biomass : std::numeric_limits<double>::quiet_NaN();
    out << data[i].year << "," << biomass << "," << q_hat * biomass << "," << depletion << "," << f_proxy << "\n";
  }
}

inline void write_pending_quantity_uncertainty_csv(const std::string& path, const std::vector<Observation>& data) {
  std::ofstream out(path);
  out << "year,quantity,estimate,se,lwr_95,upr_95,note\n";
  for (const auto& obs : data) {
    out << obs.year << ",biomass,,,,,pending delta-method propagation\n";
    out << obs.year << ",depletion,,,,,pending delta-method propagation\n";
    out << obs.year << ",F_proxy,,,,,pending delta-method propagation\n";
  }
}

inline void write_projection_uncertainty_csv(const std::string& path, const std::vector<ProjectionRow>& rows) {
  std::ofstream out(path);
  out << "scenario,year,quantity,estimate,se,lwr_95,upr_95,note\n";
  for (const auto& row : rows) {
    out << row.scenario << "," << row.year << ",biomass," << row.biomass << ",,,,pending projection covariance/simulation envelope\n";
    out << row.scenario << "," << row.year << ",index," << row.index << ",,,,pending projection covariance/simulation envelope\n";
  }
}

inline void write_runtime_memory_summary_csv(const std::string& path, double runtime_ms, std::size_t random_effects, std::size_t hessian_nonzeros) {
  std::ofstream out(path);
  out << "field,value\n";
  out << "fit_runtime_ms," << runtime_ms << "\n";
  out << "random_effects," << random_effects << "\n";
  out << "hessian_nonzeros," << hessian_nonzeros << "\n";
  out << "peak_rss_mb,\n";
  out << "note,peak RSS is captured by benchmark runner rather than model executable\n";
}

'''

m = re.search(r'(?m)^int\s+main\s*\(', s)
if not m:
    raise SystemExit("ERROR: could not find int main(")
s = s[:m.start()] + helper + "\n" + s[m.start():]

# Use the fit summary write as the anchor because projection writers have drifted.
anchor_patterns = [
    r'(write_fit_summary_csv\s*\([^;]*synthetic_fit_summary\.csv[^;]*\)\s*;\s*)',
    r'(write_summary_csv\s*\([^;]*synthetic_fit_summary\.csv[^;]*\)\s*;\s*)',
    r'(std::ofstream\s+summary[^;]*synthetic_fit_summary\.csv[\s\S]*?summary\.close\s*\(\s*\)\s*;\s*)',
]
m = None
for pat in anchor_patterns:
    m = re.search(pat, s, re.S)
    if m:
        break

if not m:
    # Last resort: put calls immediately before the "Wrote outputs" print.
    m = re.search(r'(\s*std::cout\s*<<\s*"\\nWrote outputs:)', s)
    if not m:
        raise SystemExit("ERROR: could not find fit-summary or Wrote outputs anchor")

calls = r'''
  const auto logq_uncertainty =
      compute_log_q_uncertainty_report(model, params, opts, fit);

  write_uncertainty_summary_csv("examples/opakapaka_projection/outputs/uncertainty_summary.csv", logq_uncertainty);
  write_covariance_matrix_csv("examples/opakapaka_projection/outputs/covariance_matrix.csv", logq_uncertainty);
  write_correlation_matrix_csv("examples/opakapaka_projection/outputs/correlation_matrix.csv");
  write_standard_errors_csv("examples/opakapaka_projection/outputs/standard_errors.csv", logq_uncertainty);
  write_confidence_intervals_csv("examples/opakapaka_projection/outputs/confidence_intervals.csv", logq_uncertainty);
  write_random_effect_uncertainty_csv("examples/opakapaka_projection/outputs/random_effect_uncertainty.csv", fit.u_hat);
  write_derived_quantities_csv("examples/opakapaka_projection/outputs/derived_quantities.csv", data, fit.u_hat, std::exp(fit.par.at(0)));
  write_pending_quantity_uncertainty_csv("examples/opakapaka_projection/outputs/derived_quantity_uncertainty.csv", data);
  write_projection_uncertainty_csv("examples/opakapaka_projection/outputs/projection_uncertainty.csv", projection_rows);
  write_runtime_memory_summary_csv("examples/opakapaka_projection/outputs/runtime_memory_summary.csv", fit.runtime_ms, fit.u_hat.size(), 58);

'''

s = s[:m.end()] + calls + s[m.end():]
p.write_text(s)
print("Patched opakapaka_projection.cpp")
PY

cat > inspect_opakapaka_level1_reporting_v3.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail
cpp="examples/opakapaka_projection/opakapaka_projection.cpp"

echo "== Reporting markers =="
grep -n "QUADRA_LEVEL1_UNCERTAINTY_REPORTING_V3\\|write_uncertainty_summary_csv\\|write_projection_uncertainty_csv" "$cpp"

echo
echo "== Run existing Opakapaka runner =="
./run_opakapaka_projection.sh

echo
echo "== New outputs =="
ls -1 examples/opakapaka_projection/outputs | grep -E 'uncertainty|covariance|correlation|standard_errors|confidence|derived|runtime' || true

echo
echo "== uncertainty_summary.csv =="
cat examples/opakapaka_projection/outputs/uncertainty_summary.csv

echo
echo "== standard_errors.csv =="
cat examples/opakapaka_projection/outputs/standard_errors.csv
SH
chmod +x inspect_opakapaka_level1_reporting_v3.sh

echo
echo "Backups saved with suffix: level1_reporting_v3_${stamp}.bak"
echo "Run:"
echo "  ./inspect_opakapaka_level1_reporting_v3.sh"
