#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/resource.h>
#include <vector>

#include "../core/laplace/sparse_huu_factorization.hpp"
#include "../core/laplace/third_order_dense_hdot.hpp"
#include "../examples/big/catch_at_age_shared.hpp"
#include "../include/quadra/stats.hpp"

namespace {
using Clock = std::chrono::steady_clock;

double peak_rss_mb() {
  rusage usage{};
  if (getrusage(RUSAGE_SELF, &usage) != 0)
    return std::nan("");
#ifdef __APPLE__
  return static_cast<double>(usage.ru_maxrss) / (1024.0 * 1024.0);
#else
  return static_cast<double>(usage.ru_maxrss) / 1024.0;
#endif
}

Eigen::SparseMatrix<double> retain_band(const Eigen::MatrixXd &H, int band) {
  std::vector<Eigen::Triplet<double>> entries;
  entries.reserve(static_cast<std::size_t>(H.rows() * (2 * band + 1)));
  for (int col = 0; col < H.cols(); ++col)
    for (int row = std::max(0, col - band);
         row <= std::min<int>(H.rows() - 1, col + band); ++row)
      entries.emplace_back(row, col, H(row, col));
  Eigen::SparseMatrix<double> out(H.rows(), H.cols());
  out.setFromTriplets(entries.begin(), entries.end());
  out.makeCompressed();
  return out;
}
} // namespace

int main(int argc, char **argv) {
  const int repetitions = argc > 1 ? std::max(1, std::atoi(argv[1])) : 5;
  const int band_filter = argc > 2 ? std::atoi(argv[2]) : -1;
  example::CatchAtAgeLaplaceModel model;
  quadra::ParameterSet parameters;
  std::vector<double> fixed = {
      8.199973, -0.369740, 0.005970, -1.194726, 0.293207,
      0.350404, -1.709454, -2.741071, -3.407494, 7.595303};
  if (argc == 13) {
    for (std::size_t i = 0; i < fixed.size(); ++i)
      fixed[i] = std::strtod(argv[3 + i], nullptr);
  } else if (argc != 1 && argc != 2 && argc != 3) {
    std::cerr << "usage: benchmark_big_catch_at_age_banded_approximation "
                 "[repetitions] [band|-1] [10 fixed parameters]\n";
    return 2;
  }
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
  const auto base = evaluator.evaluate(fixed);
  if (!base.converged_m || !base.logdet_ok_m) {
    std::cerr << "Laplace evaluation failed: " << base.message_m << '\n';
    return 1;
  }

  const int nf = static_cast<int>(fixed.size());
  const int nr = static_cast<int>(base.u_hat_m.size());
  const Eigen::MatrixXd H(base.hessian_random_m);
  std::vector<double> combined = fixed;
  combined.insert(combined.end(), base.u_hat_m.begin(), base.u_hat_m.end());
  std::vector<int> random_indices(static_cast<std::size_t>(nr));
  for (int i = 0; i < nr; ++i)
    random_indices[static_cast<std::size_t>(i)] = nf + i;
  auto objective = [&model](const auto &x) {
    quadra::ModelReportContext context;
    model.initialize(context);
    return model.evaluate(x, context);
  };

  const Eigen::VectorXd joint_gradient = Eigen::Map<const Eigen::VectorXd>(
      base.gradient_fixed_joint_m.data(), nf);
  auto evaluate_gradient = [&](quadra::laplace::SparseHuuFactorization &factor,
                               Eigen::MatrixXd &directions,
                               int hdot_bandwidth) {
    directions = -factor.solve(base.mixed_hessian_m);
    Eigen::VectorXd gradient = joint_gradient;
    for (int j = 0; j < nf; ++j) {
      std::vector<double> total_direction(static_cast<std::size_t>(nf + nr),
                                          0.0);
      total_direction[static_cast<std::size_t>(j)] = 1.0;
      for (int i = 0; i < nr; ++i)
        total_direction[static_cast<std::size_t>(nf + i)] = directions(i, j);
      const double trace =
          quadra::laplace::dense_hdot_trace_third_order_polarized(
              objective, combined, total_direction, random_indices,
              [&](int column) {
                Eigen::VectorXd rhs = Eigen::VectorXd::Zero(nr);
                rhs[column] = 1.0;
                return factor.solve(rhs);
              },
              hdot_bandwidth);
      gradient[j] += 0.5 * trace;
    }
    return gradient;
  };

  quadra::laplace::SparseHuuFactorization full_factor(H.sparseView());
  Eigen::MatrixXd full_directions;
  const auto reference_start = Clock::now();
  Eigen::VectorXd reference_gradient;
  for (int repetition = 0; repetition < repetitions; ++repetition)
    reference_gradient = evaluate_gradient(full_factor, full_directions, -1);
  const double reference_gradient_ms =
      std::chrono::duration<double, std::milli>(Clock::now() - reference_start)
          .count() /
      repetitions;
  const double reference_logdet = full_factor.logdet();

  std::cout << "bandwidth,nnz,logdet_abs_error,logdet_relative_error,"
               "sensitivity_relative_frobenius,reference_gradient_l2,"
               "approximate_gradient_l2,gradient_l2_error,gradient_cosine,"
               "approximate_is_exact_descent,"
               "gradient_max_abs_error,factorization_ms,gradient_ms,"
               "peak_rss_mb\n";
  std::cout << "full," << H.size() << ",0,0,0,"
            << reference_gradient.norm() << ',' << reference_gradient.norm()
            << ",0,1,1,0,0,"
            << std::fixed << std::setprecision(6) << reference_gradient_ms
            << ',' << peak_rss_mb() << '\n';

  for (int band : std::vector<int>{0, 1, 2, 3, 5, 10}) {
    if (band_filter >= 0 && band != band_filter)
      continue;
    try {
      const Eigen::SparseMatrix<double> Hband = retain_band(H, band);
      const auto factor_start = Clock::now();
      quadra::laplace::SparseHuuFactorization factor(Hband);
      const double factor_ms =
          std::chrono::duration<double, std::milli>(Clock::now() - factor_start)
              .count();
      Eigen::MatrixXd directions;
      const auto gradient_start = Clock::now();
      Eigen::VectorXd gradient;
      for (int repetition = 0; repetition < repetitions; ++repetition)
        gradient = evaluate_gradient(factor, directions, band);
      const double gradient_ms =
          std::chrono::duration<double, std::milli>(Clock::now() -
                                                     gradient_start)
              .count() /
          repetitions;
      const double logdet_error = std::abs(factor.logdet() - reference_logdet);
      const Eigen::VectorXd gradient_error = gradient - reference_gradient;
      const double cosine =
          gradient.dot(reference_gradient) /
          std::max(1.0e-30, gradient.norm() * reference_gradient.norm());
      std::cout << band << ',' << Hband.nonZeros() << ',' << logdet_error << ','
                << logdet_error / std::max(1.0, std::abs(reference_logdet))
                << ','
                << (directions - full_directions).norm() /
                       std::max(1.0, full_directions.norm())
                << ',' << reference_gradient.norm() << ','
                << gradient.norm() << ',' << gradient_error.norm() << ','
                << cosine << ','
                << (gradient.dot(reference_gradient) > 0.0 ? 1 : 0) << ','
                << gradient_error.cwiseAbs().maxCoeff() << ',' << factor_ms
                << ',' << gradient_ms << ',' << peak_rss_mb() << '\n';
    } catch (const std::exception &error) {
      std::cout << band << ",ERROR," << error.what() << '\n';
    }
  }
}
