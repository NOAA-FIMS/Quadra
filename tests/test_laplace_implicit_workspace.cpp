#include "../core/laplace/laplace_implicit_workspace.hpp"
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

  for (int y = 0; y < model.data.n_years; ++y) {
    parameters.add("rec_dev_" + std::to_string(y + 1), 0.0,
                   quadra::ParameterTransform::Identity, true);

    random_initial.push_back(0.0);
  }

  const auto workspace = quadra::build_laplace_implicit_workspace(
      model, fixed, random_initial, parameters);

  if (!workspace.success_m) {
    std::cerr << "FAIL: workspace build failed: " << workspace.message_m
              << "\n";
    return 1;
  }

  if (!workspace.factorization_m->factorized()) {
    std::cerr << "FAIL: workspace factorization not ready\n";
    return 1;
  }

  if (workspace.du_dtheta_m.rows() == 0 || workspace.du_dtheta_m.cols() == 0) {
    std::cerr << "FAIL: du_dtheta empty\n";
    return 1;
  }

  std::cout << "PASS: Laplace implicit workspace\n";
  std::cout << "  H_uu dims: " << workspace.H_uu_m.rows() << " x "
            << workspace.H_uu_m.cols() << "\n";

  std::cout << "  du_dtheta dims: " << workspace.du_dtheta_m.rows() << " x "
            << workspace.du_dtheta_m.cols() << "\n";

  std::cout << "  factorization nnz: " << workspace.factorization_m->nonzeros()
            << "\n";

  return 0;
}
