#!/usr/bin/env bash
set -euo pipefail

mkdir -p examples/surplus_production .quadra_patch_backups

bench="benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp"
if [[ -f "$bench" ]]; then
  cp "$bench" ".quadra_patch_backups/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.cleanup_debug_prints.$(date +%Y%m%d_%H%M%S).bak"

  python3 - <<'PYEOF'
from pathlib import Path

p = Path("benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp")
s = p.read_text()

for fn in [
    "PrintFlatIntermediateDirectionalCounters",
    "PrintIntermediateEdgeSlotRegistryDiagnostic",
    "PrintBatchEdgeSlotCoverageDiagnostic",
]:
    while fn in s:
        idx = s.find(fn)
        start = s.rfind("        if (m == 500", 0, idx)
        if start < 0:
            break
        end = s.find("\n\n", idx)
        if end < 0:
            break
        s = s[:start] + s[end + 2:]

p.write_text(s)
PYEOF
fi

cat > examples/surplus_production/surplus_production.hpp <<'EOF'
#ifndef QUADRA_EXAMPLES_SURPLUS_PRODUCTION_HPP
#define QUADRA_EXAMPLES_SURPLUS_PRODUCTION_HPP

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace quadra_examples {
namespace surplus_production {

struct Data {
  std::vector<double> catch_observed;
  std::vector<double> index_observed;
};

struct Parameters {
  double log_r = std::log(0.35);
  double log_K = std::log(1200.0);
  double log_q = std::log(0.0015);
  double log_sigma_index = std::log(0.18);
  double logit_B0_frac = std::log(0.85 / 0.15);
};

struct Derived {
  double r = 0.0;
  double K = 0.0;
  double q = 0.0;
  double sigma_index = 0.0;
  double B0_frac = 0.0;
  double MSY = 0.0;
  double B_MSY = 0.0;
  double F_MSY = 0.0;
  double depletion_terminal = 0.0;
  std::vector<double> biomass;
  std::vector<double> index_predicted;
  std::vector<double> log_index_residuals;
};

inline double inv_logit(const double x) {
  if (x >= 0.0) {
    const double z = std::exp(-x);
    return 1.0 / (1.0 + z);
  }
  const double z = std::exp(x);
  return z / (1.0 + z);
}

inline void validate_data(const Data& data) {
  if (data.catch_observed.empty()) {
    throw std::runtime_error("catch_observed must not be empty");
  }
  if (data.index_observed.size() != data.catch_observed.size()) {
    throw std::runtime_error("index_observed and catch_observed must have same length");
  }
  for (const double c : data.catch_observed) {
    if (!(c >= 0.0)) throw std::runtime_error("catch_observed must be nonnegative");
  }
  for (const double i : data.index_observed) {
    if (!(i > 0.0)) throw std::runtime_error("index_observed must be positive");
  }
}

inline Derived evaluate_derived(const Data& data, const Parameters& par) {
  validate_data(data);

  Derived out;
  out.r = std::exp(par.log_r);
  out.K = std::exp(par.log_K);
  out.q = std::exp(par.log_q);
  out.sigma_index = std::exp(par.log_sigma_index);
  out.B0_frac = inv_logit(par.logit_B0_frac);

  const int n = static_cast<int>(data.catch_observed.size());
  out.biomass.assign(static_cast<std::size_t>(n + 1), 0.0);
  out.index_predicted.assign(static_cast<std::size_t>(n), 0.0);
  out.log_index_residuals.assign(static_cast<std::size_t>(n), 0.0);

  out.biomass[0] = out.B0_frac * out.K;

  for (int t = 0; t < n; ++t) {
    const double B_t = std::max(out.biomass[static_cast<std::size_t>(t)], 1e-12);
    out.index_predicted[static_cast<std::size_t>(t)] = out.q * B_t;
    out.log_index_residuals[static_cast<std::size_t>(t)] =
        std::log(data.index_observed[static_cast<std::size_t>(t)]) -
        std::log(out.index_predicted[static_cast<std::size_t>(t)]);

    const double production = out.r * B_t * (1.0 - B_t / out.K);
    const double next_B = B_t + production - data.catch_observed[static_cast<std::size_t>(t)];
    out.biomass[static_cast<std::size_t>(t + 1)] = std::max(next_B, 1e-9);
  }

  out.MSY = out.r * out.K / 4.0;
  out.B_MSY = out.K / 2.0;
  out.F_MSY = out.r / 2.0;
  out.depletion_terminal = out.biomass.back() / out.K;

  return out;
}

inline double negative_log_likelihood(const Data& data, const Parameters& par) {
  const Derived d = evaluate_derived(data, par);
  double nll = 0.0;
  const double sigma = d.sigma_index;
  const double log_norm = std::log(sigma) + 0.5 * std::log(2.0 * M_PI);

  for (const double residual : d.log_index_residuals) {
    const double z = residual / sigma;
    nll += log_norm + 0.5 * z * z;
  }

  return nll;
}

inline Data make_demo_data() {
  Data data;
  data.catch_observed = {
      80.0, 88.0, 95.0, 105.0, 115.0, 125.0, 130.0, 128.0,
      120.0, 110.0, 100.0, 90.0, 82.0, 78.0, 75.0};
  data.index_observed = {
      1.55, 1.50, 1.43, 1.34, 1.22, 1.10, 0.98, 0.88,
      0.81, 0.78, 0.77, 0.80, 0.84, 0.88, 0.92};
  return data;
}

inline Parameters make_demo_parameters() {
  Parameters par;
  par.log_r = std::log(0.38);
  par.log_K = std::log(1250.0);
  par.log_q = std::log(0.00145);
  par.log_sigma_index = std::log(0.16);
  par.logit_B0_frac = std::log(0.90 / 0.10);
  return par;
}

inline void print_report(const Data& data, const Parameters& par) {
  const Derived d = evaluate_derived(data, par);
  const double nll = negative_log_likelihood(data, par);

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "Surplus production example\n";
  std::cout << "--------------------------\n";
  std::cout << "negative log likelihood = " << nll << "\n\n";

  std::cout << "Parameters\n";
  std::cout << "  r            = " << d.r << "\n";
  std::cout << "  K            = " << d.K << "\n";
  std::cout << "  q            = " << d.q << "\n";
  std::cout << "  sigma_index  = " << d.sigma_index << "\n";
  std::cout << "  B0 / K       = " << d.B0_frac << "\n\n";

  std::cout << "Reference points\n";
  std::cout << "  MSY          = " << d.MSY << "\n";
  std::cout << "  B_MSY        = " << d.B_MSY << "\n";
  std::cout << "  F_MSY        = " << d.F_MSY << "\n";
  std::cout << "  B_terminal/K = " << d.depletion_terminal << "\n\n";

  std::cout << std::setw(6) << "year"
            << std::setw(14) << "catch"
            << std::setw(14) << "biomass"
            << std::setw(14) << "index obs"
            << std::setw(14) << "index pred"
            << std::setw(14) << "log resid"
            << "\n";

  for (std::size_t t = 0; t < data.catch_observed.size(); ++t) {
    std::cout << std::setw(6) << t
              << std::setw(14) << data.catch_observed[t]
              << std::setw(14) << d.biomass[t]
              << std::setw(14) << data.index_observed[t]
              << std::setw(14) << d.index_predicted[t]
              << std::setw(14) << d.log_index_residuals[t]
              << "\n";
  }

  std::cout << std::setw(6) << data.catch_observed.size()
            << std::setw(14) << "-"
            << std::setw(14) << d.biomass.back()
            << std::setw(14) << "-"
            << std::setw(14) << "-"
            << std::setw(14) << "-"
            << "\n";
}

}  // namespace surplus_production
}  // namespace quadra_examples

#endif
EOF

cat > examples/surplus_production/run_surplus_production.cpp <<'EOF'
#include "surplus_production.hpp"

int main() {
  const auto data = quadra_examples::surplus_production::make_demo_data();
  const auto par = quadra_examples::surplus_production::make_demo_parameters();
  quadra_examples::surplus_production::print_report(data, par);
  return 0;
}
EOF

cat > examples/surplus_production/README.md <<'EOF'
# Surplus production example

This is a deterministic Schaefer surplus production model:

```text
B[t+1] = B[t] + r B[t] (1 - B[t] / K) - C[t]
I[t]   = q B[t] exp(epsilon[t])
```

The example reports the objective, biomass trajectory, fitted index, residuals, and reference points:

```text
MSY = rK / 4
B_MSY = K / 2
F_MSY = r / 2
terminal depletion = B_terminal / K
```

Run from the repository root:

```bash
./run_surplus_production_example.sh
```

Suggested next steps:

1. Add finite-difference gradient checks.
2. Add a simple optimizer wrapper.
3. Add process error.
4. Add random effects and Laplace evaluation.
EOF

cat > run_surplus_production_example.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} -Iexamples/surplus_production   examples/surplus_production/run_surplus_production.cpp   -o build/examples/run_surplus_production

./build/examples/run_surplus_production
EOF

chmod +x run_surplus_production_example.sh

cat <<'EOF'

Installed surplus production example and cleaned temporary benchmark prints.

Added:
  examples/surplus_production/surplus_production.hpp
  examples/surplus_production/run_surplus_production.cpp
  examples/surplus_production/README.md
  run_surplus_production_example.sh

Run:
  ./run_surplus_production_example.sh

Optional sanity check:
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

EOF
