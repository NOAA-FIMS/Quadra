#!/usr/bin/env bash
set -euo pipefail

echo "== Add SEFSC red snapper Quadra optimizer adapter =="

BASE="examples/NMFS/sefsc_red_snapper"
mkdir -p "$BASE"/{quadra,outputs,validation}

cat > "$BASE/quadra/red_snapper_quadra_fit.cpp" <<'CPP'
#include "red_snapper_age_structured.hpp"

#include "../../../core/optimizer.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace sefsc_red_snapper {

template <class T>
T exp_t(const T& x) {
  using std::exp;
  return exp(x);
}

template <class T>
T log_t(const T& x) {
  using std::log;
  return log(x);
}

template <class T>
T max_t(const T& x, double floor) {
  return x > T(floor) ? x : T(floor);
}

template <class T>
T square_t(const T& x) {
  return x * x;
}

template <class T>
T logistic_selectivity_t(const T& age, const T& a50, const T& slope) {
  return T(1.0) / (T(1.0) + exp_t(-slope * (age - a50)));
}

class RedSnapperQuadraObjective {
 public:
  explicit RedSnapperQuadraObjective(std::vector<Observation> observations)
      : observations_(std::move(observations)) {}

  template <class T>
  T operator()(const std::vector<T>& par) const {
    if (par.size() < 3) {
      throw std::runtime_error(
          "RedSnapperQuadraObjective expected parameters: log_r0, log_fbar, log_q");
    }

    const T log_r0 = par[0];
    const T log_fbar = par[1];
    const T log_q = par[2];

    const T r0 = exp_t(log_r0);
    const T m = T(0.18);
    const T fbar = exp_t(log_fbar);
    const T q = exp_t(log_q);

    const T sigma_log_index = T(0.20);
    const T sigma_log_catch = T(0.15);
    const double min_positive = 1.0e-12;

    const auto weight = default_weight_at_age();
    const auto maturity = default_maturity_at_age();

    std::array<T, kAges> selectivity{};
    for (int a = 0; a < kAges; ++a) {
      selectivity[static_cast<std::size_t>(a)] =
          logistic_selectivity_t(T(a + 1), T(4.0), T(1.2));
    }

    std::array<T, kAges> n{};
    n[0] = r0;
    for (int a = 1; a < kAges; ++a) {
      n[static_cast<std::size_t>(a)] =
          n[static_cast<std::size_t>(a - 1)] * exp_t(-m);
    }
    n[static_cast<std::size_t>(kAges - 1)] =
        n[static_cast<std::size_t>(kAges - 1)] /
        (T(1.0) - exp_t(-m));

    T nll = T(0.0);

    for (const auto& obs : observations_) {
      T biomass = T(0.0);
      for (int a = 0; a < kAges; ++a) {
        biomass = biomass +
                  n[static_cast<std::size_t>(a)] *
                      T(weight[static_cast<std::size_t>(a)]);
      }

      T catch_hat = T(0.0);
      for (int a = 0; a < kAges; ++a) {
        const auto i = static_cast<std::size_t>(a);
        const T f_a = fbar * selectivity[i];
        const T z_a = m + f_a;
        const T harvest_rate =
            (f_a / z_a) * (T(1.0) - exp_t(-z_a));
        catch_hat = catch_hat + n[i] * T(weight[i]) * harvest_rate;
      }

      const T index_hat = q * biomass;

      if (obs.index > 0.0) {
        const T z = (log_t(T(obs.index)) -
                     log_t(max_t(index_hat, min_positive))) /
                    sigma_log_index;
        nll = nll + T(0.5) * square_t(z);
      }

      if (obs.catch_mt > 0.0) {
        const T z = (log_t(T(obs.catch_mt)) -
                     log_t(max_t(catch_hat, min_positive))) /
                    sigma_log_catch;
        nll = nll + T(0.5) * square_t(z);
      }

      std::array<T, kAges> next{};
      next[0] = r0;

      for (int a = 1; a < kAges; ++a) {
        const auto prev = static_cast<std::size_t>(a - 1);
        const T f_prev = fbar * selectivity[prev];
        const T z_prev = m + f_prev;
        next[static_cast<std::size_t>(a)] = n[prev] * exp_t(-z_prev);
      }

      const auto last = static_cast<std::size_t>(kAges - 1);
      const T f_last = fbar * selectivity[last];
      const T z_last = m + f_last;
      next[last] = next[last] + n[last] * exp_t(-z_last);

      n = next;
    }

    return nll;
  }

 private:
  std::vector<Observation> observations_;
};

void write_fit_summary(const std::string& path,
                       const quadra::OptResult& fit) {
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Could not open fit summary CSV: " + path);
  }

  out << "field,value\n";
  out << std::setprecision(12);
  out << "objective," << fit.value << "\n";
  out << "grad_norm," << fit.grad_norm << "\n";
  out << "iterations," << fit.iterations << "\n";
  out << "converged," << (fit.converged ? "yes" : "no") << "\n";
  out << "message," << fit.message << "\n";

  if (fit.par.size() >= 3) {
    out << "log_r0," << fit.par[0] << "\n";
    out << "r0," << std::exp(fit.par[0]) << "\n";
    out << "log_fbar," << fit.par[1] << "\n";
    out << "fbar," << std::exp(fit.par[1]) << "\n";
    out << "log_q," << fit.par[2] << "\n";
    out << "q," << std::exp(fit.par[2]) << "\n";
  }
}

}  // namespace sefsc_red_snapper

int main() {
  const std::string input_path =
      "examples/NMFS/sefsc_red_snapper/data/synthetic_red_snapper_observations.csv";
  const std::string summary_path =
      "examples/NMFS/sefsc_red_snapper/outputs/quadra_fit_summary.csv";

  const auto observations = sefsc_red_snapper::read_observations(input_path);

  sefsc_red_snapper::RedSnapperQuadraObjective objective(observations);

  quadra::ParameterVector params;
  params.add_fixed("log_r0", std::log(1200.0));
  params.add_fixed("log_fbar", std::log(0.025));
  params.add_fixed("log_q", std::log(0.00005));

  quadra::LaplaceOptions opts;
  opts.max_outer_iterations = 200;
  opts.fixed_grad_tol = 1.0e-8;

  auto fit = quadra::optimize_lbfgs(objective, params, opts);

  sefsc_red_snapper::write_fit_summary(summary_path, fit);

  std::cout << "SEFSC red-snapper-style Quadra fixed-effect fit\n";
  std::cout << "objective:  " << fit.value << "\n";
  std::cout << "grad_norm:  " << fit.grad_norm << "\n";
  std::cout << "converged:  " << (fit.converged ? "yes" : "no") << "\n";
  std::cout << "message:    " << fit.message << "\n";
  std::cout << "wrote:      " << summary_path << "\n";

  return 0;
}
CPP

cat > "$BASE/run_red_snapper_quadra_fit.sh" <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p examples/NMFS/sefsc_red_snapper/outputs

c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  -Icore \
  -Iexternal/LBFGSpp/include \
  -o examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit \
  examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit.cpp

./examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit
SH
chmod +x "$BASE/run_red_snapper_quadra_fit.sh"

cat > "$BASE/validation/quadra_fit_checklist.md" <<'MD'
# Quadra Fit Checklist

- [x] fixed-effect objective adapter added
- [x] Quadra optimizer path used
- [ ] fixed-effect fit compiles against current local API
- [ ] fit summary output validated
- [ ] derived trajectory generated at fitted parameters
- [ ] age-composition likelihood added
- [ ] recruitment deviations added as random effects
- [ ] Laplace uncertainty added
MD

echo
echo "Added Quadra optimizer adapter."
echo
echo "Run:"
echo "  ./examples/NMFS/sefsc_red_snapper/run_red_snapper_quadra_fit.sh"
echo
echo "If the compile fails, inspect the current optimizer API with:"
echo "  grep -R \"add_fixed\\|optimize_lbfgs\\|struct OptResult\\|struct LaplaceOptions\" -n core examples | head -80"
