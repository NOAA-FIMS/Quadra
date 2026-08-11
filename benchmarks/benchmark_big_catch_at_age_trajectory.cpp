#include <iomanip>
#include <iostream>
#include <vector>

#include "../core/optimizer.hpp"
#include "../examples/big/catch_at_age_shared.hpp"

int main() {
  example::CatchAtAgeLaplaceModel model;
  quadra::ParameterVector parameters;
  const std::vector<double> initial = {
      std::log(900.0), std::log(0.25), std::log(0.15), std::log(0.18), 0.0,
      std::log(1.25), std::log(0.20), std::log(0.15), std::log(0.35),
      std::log(40.0)};
  const std::vector<const char *> names = {
      "log_R0",          "log_M",           "log_q",
      "log_Fbar",        "sel50_raw",       "log_sel_slope",
      "log_sigma_index", "log_sigma_catch", "log_sigma_rec",
      "log_comp_concentration"};
  for (std::size_t i = 0; i < initial.size(); ++i)
    parameters.add({names[i], initial[i], quadra::ParameterTransform::Identity,
                    false});
  for (int year = 0; year < model.data.n_years; ++year)
    parameters.add({"rec_dev_" + std::to_string(year + 1), 0.0,
                    quadra::ParameterTransform::Identity, true});

  quadra::LaplaceOptions options = quadra::default_laplace_options();
  options.use_hutchinson_trace = false;
  options.hessian_drop_tol = 0.0;
  const auto fit = quadra::optimize_lbfgs(model, parameters, options);

  std::cout << "TRAJECTORY_CSV\n";
  std::cout << "evaluation,objective,gradient_norm";
  for (const char *name : names)
    std::cout << ',' << name;
  std::cout << '\n' << std::setprecision(17);
  for (const auto &point : fit.trajectory) {
    std::cout << point.evaluation << ',' << point.objective << ','
              << point.gradient_norm;
    for (double value : point.fixed)
      std::cout << ',' << value;
    std::cout << '\n';
  }
}
