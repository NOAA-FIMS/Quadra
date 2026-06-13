#!/usr/bin/env bash
set -euo pipefail

echo "== Add SEFSC red snapper objective function scaffold =="

BASE="examples/NMFS/sefsc_red_snapper"
mkdir -p "$BASE"/{quadra,outputs,validation}

cat > "$BASE/quadra/red_snapper_objective.hpp" <<'CPP'
#pragma once

#include "red_snapper_age_structured.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace sefsc_red_snapper {

struct ObjectiveOptions {
  double sigma_log_index = 0.20;
  double sigma_log_catch = 0.15;
  double min_positive = 1.0e-12;
};

struct ObjectiveBreakdown {
  double total = 0.0;
  double index_nll = 0.0;
  double catch_nll = 0.0;
  int n_index = 0;
  int n_catch = 0;
};

inline double square(double x) { return x * x; }

inline double lognormal_nll_no_constant(double observed, double predicted,
                                        double sigma, double min_positive) {
  const double obs = std::max(observed, min_positive);
  const double pred = std::max(predicted, min_positive);
  const double z = (std::log(obs) - std::log(pred)) / sigma;
  return 0.5 * square(z);
}

inline ObjectiveBreakdown evaluate_objective_breakdown(
    const std::vector<Observation>& observations,
    const AgeStructuredParams& params,
    const ObjectiveOptions& options = ObjectiveOptions{}) {
  ObjectiveBreakdown out;

  const auto rows = run_deterministic_age_structured_model(observations, params);
  if (rows.size() != observations.size()) {
    throw std::runtime_error("Objective trajectory/observation size mismatch");
  }

  for (std::size_t i = 0; i < observations.size(); ++i) {
    const auto& obs = observations[i];
    const auto& pred = rows[i];

    if (std::isfinite(obs.index) && obs.index > 0.0) {
      const double nll = lognormal_nll_no_constant(
          obs.index, pred.index_hat, options.sigma_log_index,
          options.min_positive);
      out.index_nll += nll;
      ++out.n_index;
    }

    if (std::isfinite(obs.catch_mt) && obs.catch_mt > 0.0) {
      const double nll = lognormal_nll_no_constant(
          obs.catch_mt, pred.catch_hat, options.sigma_log_catch,
          options.min_positive);
      out.catch_nll += nll;
      ++out.n_catch;
    }
  }

  out.total = out.index_nll + out.catch_nll;
  return out;
}

inline double evaluate_objective(
    const std::vector<Observation>& observations,
    const AgeStructuredParams& params,
    const ObjectiveOptions& options = ObjectiveOptions{}) {
  return evaluate_objective_breakdown(observations, params, options).total;
}

}  // namespace sefsc_red_snapper
CPP

# Split reusable age-structured model logic into a header if it does not exist yet.
if [[ ! -f "$BASE/quadra/red_snapper_age_structured.hpp" ]]; then
  python3 - <<'PY'
from pathlib import Path

cpp = Path("examples/NMFS/sefsc_red_snapper/quadra/red_snapper_age_structured.cpp")
hpp = Path("examples/NMFS/sefsc_red_snapper/quadra/red_snapper_age_structured.hpp")

s = cpp.read_text()
marker = "}  // namespace sefsc_red_snapper\n\nint main()"
idx = s.find(marker)
if idx < 0:
    raise SystemExit("Could not split red_snapper_age_structured.cpp into header")

header_body = s[: idx + len("}  // namespace sefsc_red_snapper\n")]
header_body = header_body.replace('#include <iostream>\n', '')
header = "#pragma once\n\n" + header_body
hpp.write_text(header)

main_part = s[idx + len("}  // namespace sefsc_red_snapper\n\n"):]
new_cpp = '#include "red_snapper_age_structured.hpp"\n\n#include <iostream>\n\n' + main_part
cpp.write_text(new_cpp)
print("created", hpp)
print("rewrote", cpp)
PY
fi

cat > "$BASE/quadra/evaluate_red_snapper_objective.cpp" <<'CPP'
#include "red_snapper_objective.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void write_objective_summary(
    const std::string& path,
    const sefsc_red_snapper::ObjectiveBreakdown& obj,
    const sefsc_red_snapper::AgeStructuredParams& params) {
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Could not open objective summary CSV: " + path);
  }

  out << "field,value\n";
  out << std::setprecision(12);
  out << "objective_total," << obj.total << "\n";
  out << "index_nll," << obj.index_nll << "\n";
  out << "catch_nll," << obj.catch_nll << "\n";
  out << "n_index," << obj.n_index << "\n";
  out << "n_catch," << obj.n_catch << "\n";
  out << "log_r0," << params.log_r0 << "\n";
  out << "r0," << std::exp(params.log_r0) << "\n";
  out << "log_m," << params.log_m << "\n";
  out << "m," << std::exp(params.log_m) << "\n";
  out << "log_fbar," << params.log_fbar << "\n";
  out << "fbar," << std::exp(params.log_fbar) << "\n";
  out << "log_q," << params.log_q << "\n";
  out << "q," << std::exp(params.log_q) << "\n";
  out << "sel_a50," << params.sel_a50 << "\n";
  out << "sel_slope," << params.sel_slope << "\n";
}

}  // namespace

int main() {
  const std::string input_path =
      "examples/NMFS/sefsc_red_snapper/data/synthetic_red_snapper_observations.csv";
  const std::string summary_path =
      "examples/NMFS/sefsc_red_snapper/outputs/objective_summary.csv";

  const auto observations = sefsc_red_snapper::read_observations(input_path);

  sefsc_red_snapper::AgeStructuredParams params;
  sefsc_red_snapper::ObjectiveOptions options;

  const auto breakdown =
      sefsc_red_snapper::evaluate_objective_breakdown(observations, params,
                                                      options);

  write_objective_summary(summary_path, breakdown, params);

  std::cout << "SEFSC red-snapper-style objective scaffold\n";
  std::cout << "objective_total: " << breakdown.total << "\n";
  std::cout << "index_nll:       " << breakdown.index_nll << "\n";
  std::cout << "catch_nll:       " << breakdown.catch_nll << "\n";
  std::cout << "wrote:           " << summary_path << "\n";

  return 0;
}
CPP

cat > "$BASE/run_red_snapper_objective.sh" <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p examples/NMFS/sefsc_red_snapper/outputs

c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  -Icore \
  -o examples/NMFS/sefsc_red_snapper/quadra/evaluate_red_snapper_objective \
  examples/NMFS/sefsc_red_snapper/quadra/evaluate_red_snapper_objective.cpp

./examples/NMFS/sefsc_red_snapper/quadra/evaluate_red_snapper_objective
SH
chmod +x "$BASE/run_red_snapper_objective.sh"

cat > "$BASE/validation/objective_checklist.md" <<'MD'
# Objective Function Checklist

- [x] lognormal index likelihood
- [x] lognormal catch likelihood
- [x] objective breakdown output
- [x] reusable objective header
- [ ] parameter optimization
- [ ] age-composition likelihood
- [ ] recruitment-deviation prior
- [ ] TMB objective parity
MD

echo
echo "Added objective scaffold."
echo
echo "Run:"
echo "  ./examples/NMFS/sefsc_red_snapper/run_red_snapper_objective.sh"
echo "  cat examples/NMFS/sefsc_red_snapper/outputs/objective_summary.csv"
