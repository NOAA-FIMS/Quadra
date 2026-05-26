#include "../core/laplace/laplace_implicit_derivatives.hpp"
#include "../examples/big/catch_at_age_shared.hpp"

#include <cmath>
#include <iostream>
#include <vector>

int main() {
  example::CatchAtAgeLaplaceModel model;

  quadra::ParameterSet parameters;
  parameters.add("log_R0", std::log(900.0),
                 quadra::ParameterTransform::Identity, false);
  parameters.add("log_M", std::log(0.25), quadra::ParameterTransform::Identity,
                 false);
  parameters.add("log_q", std::log(0.15), quadra::ParameterTransform::Identity,
                 false);
  parameters.add("log_Fbar", std::log(0.18),
                 quadra::ParameterTransform::Identity, false);
  parameters.add("sel50_raw", 0.0, quadra::ParameterTransform::Identity, false);
  parameters.add("log_sel_slope", std::log(1.25),
                 quadra::ParameterTransform::Identity, false);
  parameters.add("log_sigma_index", std::log(0.20),
                 quadra::ParameterTransform::Identity, false);
  parameters.add("log_sigma_catch", std::log(0.18),
                 quadra::ParameterTransform::Identity, false);
  parameters.add("log_sigma_rec", std::log(0.35),
                 quadra::ParameterTransform::Identity, false);

  std::vector<double> fixed{
      std::log(900.0), std::log(0.25), std::log(0.15), std::log(0.18), 0.0,
      std::log(1.25),  std::log(0.20), std::log(0.18), std::log(0.35)};

  std::vector<double> random_initial;
  random_initial.reserve(static_cast<std::size_t>(model.data.n_years));

  for (int y = 0; y < model.data.n_years; ++y) {
    parameters.add("rec_dev_" + std::to_string(y + 1), 0.0,
                   quadra::ParameterTransform::Identity, true);

    random_initial.push_back(0.0);
  }

  const auto partition = quadra::partition_parameters(parameters);

  const auto implicit = quadra::evaluate_laplace_implicit_derivatives(
      model, fixed, random_initial, parameters);

  if (!implicit.success_m) {
    std::cerr << "FAIL: implicit derivative solve failed: "
              << implicit.message_m << "\n";
    return 1;
  }

  const auto H_ad =
      quadra::ad_mixed_hessian(model, fixed, implicit.u_hat_m, partition);

  const auto H_fd = quadra::finite_difference_mixed_hessian(
      model, fixed, implicit.u_hat_m, partition, 1.0e-5);

  if (H_ad.rows() != H_fd.rows() || H_ad.cols() != H_fd.cols()) {
    std::cerr << "FAIL: dimension mismatch\n";
    return 1;
  }

  const double max_abs_error = (H_ad - H_fd).cwiseAbs().maxCoeff();

  const double max_abs_ad = H_ad.cwiseAbs().maxCoeff();

  const double relative_error =
      max_abs_ad > 0.0 ? max_abs_error / max_abs_ad : max_abs_error;

  if (!std::isfinite(max_abs_error) || !std::isfinite(relative_error)) {
    std::cerr << "FAIL: non-finite mixed Hessian comparison\n";
    return 1;
  }

  if (max_abs_error > 1.0e-3 && relative_error > 1.0e-4) {
    std::cerr << "FAIL: AD mixed Hessian differs from FD\n";
    std::cerr << "max_abs_error: " << max_abs_error << "\n";
    std::cerr << "relative_error: " << relative_error << "\n";
    std::cerr << "max_abs_ad: " << max_abs_ad << "\n";
    return 1;
  }

  std::cout << "PASS: AD mixed Hessian matches finite difference\n";
  std::cout << "  dims: " << H_ad.rows() << " x " << H_ad.cols() << "\n";
  std::cout << "  max_abs_error: " << max_abs_error << "\n";
  std::cout << "  relative_error: " << relative_error << "\n";
  std::cout << "  max_abs_ad: " << max_abs_ad << "\n";

  return 0;
}
