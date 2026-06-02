#include "../core/laplace/laplace_implicit_derivatives.hpp"
#include "../examples/big/catch_at_age_shared.hpp"

#include <Eigen/Dense>

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

  const auto implicit = quadra::evaluate_laplace_implicit_derivatives(
      model, fixed, random_initial, parameters);

  if (!implicit.success_m) {
    std::cerr << "FAIL: implicit derivatives failed: " << implicit.message_m
              << "\n";
    return 1;
  }

  const int n_random = model.data.n_years;
  const int n_fixed = 9;

  if (implicit.du_dtheta_m.rows() != n_random ||
      implicit.du_dtheta_m.cols() != n_fixed) {
    std::cerr << "FAIL: unexpected du_dtheta dimensions: "
              << implicit.du_dtheta_m.rows() << " x "
              << implicit.du_dtheta_m.cols() << "\n";
    return 1;
  }

  if (!implicit.du_dtheta_m.allFinite()) {
    std::cerr << "FAIL: non-finite du_dtheta\n";
    return 1;
  }

  std::cout << "PASS: big Laplace implicit derivative smoke test\n";
  std::cout << "  H_uu dims: " << implicit.H_uu_m.rows() << " x "
            << implicit.H_uu_m.cols() << "\n";
  std::cout << "  H_u_theta dims: " << implicit.H_u_theta_m.rows() << " x "
            << implicit.H_u_theta_m.cols() << "\n";
  std::cout << "  du_dtheta dims: " << implicit.du_dtheta_m.rows() << " x "
            << implicit.du_dtheta_m.cols() << "\n";
  std::cout << "  du_dtheta max abs: "
            << implicit.du_dtheta_m.cwiseAbs().maxCoeff() << "\n";

  return 0;
}
