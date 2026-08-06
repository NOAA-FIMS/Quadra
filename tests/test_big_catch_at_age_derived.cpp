#include "../core/laplace/laplace_profiled_ad_gradient.hpp"
#include "../examples/big/catch_at_age_derived.hpp"

#include <cmath>
#include <iostream>
#include <vector>

int main() {
  example::CatchAtAgeLaplaceModel model;

  std::vector<double> theta{
      std::log(900.0), std::log(0.25), std::log(0.15), std::log(0.18), 0.0,
      std::log(1.25),  std::log(0.20), std::log(0.18), std::log(0.35)};

  std::vector<double> u(static_cast<std::size_t>(model.data.n_years), 0.0);

  auto derived =
      example::evaluate_catch_at_age_derived_quantities(model, theta, u);

  if (!std::isfinite(derived.terminal_ssb_proxy_m) ||
      !std::isfinite(derived.terminal_depletion_m) ||
      !std::isfinite(derived.mean_f_m)) {
    std::cerr << "FAIL: non-finite derived quantity\n";
    return 1;
  }

  if (!(derived.terminal_ssb_proxy_m > 0.0)) {
    std::cerr << "FAIL: terminal SSB proxy is not positive\n";
    return 1;
  }

  if (!(derived.terminal_depletion_m > 0.0)) {
    std::cerr << "FAIL: terminal depletion is not positive\n";
    return 1;
  }

  if (!(derived.mean_f_m > 0.0)) {
    std::cerr << "FAIL: mean F is not positive\n";
    return 1;
  }

  const auto depletion_profiled = quadra::evaluate_profiled_ad_gradient_blocks(
      [&model](const auto &fixed, const auto &random) {
        return example::evaluate_terminal_depletion_ad(model, fixed, random);
      },
      theta, u);
  const auto ssb_profiled = quadra::evaluate_profiled_ad_gradient_blocks(
      [&model](const auto &fixed, const auto &random) {
        return example::evaluate_terminal_ssb_proxy_ad(model, fixed, random);
      },
      theta, u);
  if (!depletion_profiled.success_m || !ssb_profiled.success_m ||
      std::abs(depletion_profiled.estimate_m - derived.terminal_depletion_m) >
          1.0e-10 ||
      std::abs(ssb_profiled.estimate_m - derived.terminal_ssb_proxy_m) >
          1.0e-8) {
    std::cerr << "FAIL: profiled derived estimates do not match direct "
                 "derived estimates\n";
    return 1;
  }
  if (!(depletion_profiled.gradient_random_m.norm() > 0.0) ||
      !(ssb_profiled.gradient_random_m.norm() > 0.0)) {
    std::cerr << "FAIL: profiled derived quantities discarded random effects\n";
    return 1;
  }

  std::cout << "PASS: big catch-at-age derived quantities\n";
  std::cout << "  terminal_ssb_proxy: " << derived.terminal_ssb_proxy_m << "\n";
  std::cout << "  terminal_depletion: " << derived.terminal_depletion_m << "\n";
  std::cout << "  mean_f: " << derived.mean_f_m << "\n";
  std::cout << "  depletion random-gradient norm: "
            << depletion_profiled.gradient_random_m.norm() << "\n";

  return 0;
}
