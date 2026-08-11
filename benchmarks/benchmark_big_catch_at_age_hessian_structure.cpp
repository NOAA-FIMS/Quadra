#include <Eigen/Dense>
#include <Eigen/SVD>

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../examples/big/catch_at_age_shared.hpp"
#include "../include/quadra/stats.hpp"

int main() {
  using Clock = std::chrono::steady_clock;
  example::CatchAtAgeLaplaceModel model;
  quadra::ParameterSet parameters;
  const std::vector<double> fixed = {
      8.199973, -0.369740, 0.005970, -1.194726, 0.293207,
      0.350404, -1.709454, -2.741071, -3.407494, 7.595303};
  const std::vector<const char *> names = {
      "log_R0",          "log_M",           "log_q",
      "log_Fbar",        "sel50_raw",       "log_sel_slope",
      "log_sigma_index", "log_sigma_catch", "log_sigma_rec",
      "log_comp_concentration"};
  for (std::size_t i = 0; i < fixed.size(); ++i)
    parameters.add(names[i], fixed[i], quadra::ParameterTransform::Identity,
                   false);
  std::vector<double> random(static_cast<std::size_t>(model.data.n_years), 0.0);
  for (int year = 0; year < model.data.n_years; ++year)
    parameters.add("rec_dev_" + std::to_string(year + 1), 0.0,
                   quadra::ParameterTransform::Identity, true);

  quadra::LaplaceObjectiveOptions options;
  options.include_constant_m = true;
  options.compute_mixed_derivatives_m = true;
  options.newton_m.gradient_tolerance_m = 1.0e-6;
  quadra::stats::LaplaceEvaluator<example::CatchAtAgeLaplaceModel> evaluator(
      model, random, parameters, options);
  const auto start = Clock::now();
  const auto result = evaluator.evaluate(fixed);
  const double evaluation_ms =
      std::chrono::duration<double, std::milli>(Clock::now() - start).count();
  if (!result.converged_m || !result.logdet_ok_m) {
    std::cerr << "Laplace evaluation failed: " << result.message_m << '\n';
    return 1;
  }

  const Eigen::MatrixXd H(result.hessian_random_m);
  const double hnorm = H.norm();
  std::cout << "n,nnz,density,evaluation_ms,hessian_frobenius\n";
  std::cout << H.rows() << ',' << result.hessian_random_m.nonZeros() << ','
            << std::fixed << std::setprecision(6)
            << static_cast<double>(result.hessian_random_m.nonZeros()) /
                   static_cast<double>(H.size())
            << ',' << evaluation_ms << ',' << hnorm << "\n\n";
  std::cout << "bandwidth,residual_relative_frobenius,residual_rank_1e-6,"
               "residual_rank_1e-8,leading_residual_singular_value\n";

  for (int bandwidth : std::vector<int>{0, 1, 2, 3, 5, 10}) {
    Eigen::MatrixXd residual = H;
    for (int row = 0; row < H.rows(); ++row)
      for (int col = 0; col < H.cols(); ++col)
        if (std::abs(row - col) <= bandwidth)
          residual(row, col) = 0.0;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(residual);
    const Eigen::VectorXd singular = svd.singularValues();
    const double leading = singular.size() == 0 ? 0.0 : singular[0];
    int rank_1e6 = 0;
    int rank_1e8 = 0;
    for (Eigen::Index i = 0; i < singular.size(); ++i) {
      rank_1e6 += singular[i] > 1.0e-6 * std::max(1.0, leading);
      rank_1e8 += singular[i] > 1.0e-8 * std::max(1.0, leading);
    }
    std::cout << bandwidth << ',' << residual.norm() / hnorm << ','
              << rank_1e6 << ',' << rank_1e8 << ',' << leading << '\n';
  }
}
