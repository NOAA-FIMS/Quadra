#pragma once

#include "random_effect_hessian.hpp"
#include "third_order_dense_hdot.hpp"

#include <Eigen/Dense>
#include <Eigen/SparseCholesky>

#include <chrono>
#include <stdexcept>
#include <vector>

namespace quadra {
namespace laplace {

struct ExactLaplaceDirectionalCurvatureResult {
  double curvature_m = 0.0;
  double joint_profile_curvature_m = 0.0;
  double logdet_curvature_m = 0.0;
  Eigen::VectorXd du_m;
  Eigen::VectorXd d2u_m;
  double total_ms_m = 0.0;
};

struct ExactLaplaceHessianResult {
  Eigen::MatrixXd hessian_m;
  int directional_evaluations_m = 0;
  double total_ms_m = 0.0;
};

// Exact second directional derivative of
// f(theta,u_hat(theta)) + 0.5 log det f_uu(theta,u_hat(theta)).
// This reference implementation uses cubic and quartic polarization. It is
// intentionally matrix-free in fixed-effect space, but is expensive in the
// number of random effects; it provides the correctness oracle for a future
// reverse fourth-order contraction backend.
template <class Model>
ExactLaplaceDirectionalCurvatureResult exact_laplace_directional_curvature(
    Model &model, const std::vector<double> &fixed,
    const std::vector<double> &u_hat, const ParameterPartition &partition,
    const Eigen::VectorXd &fixed_direction) {
  using Clock = std::chrono::steady_clock;
  const auto start = Clock::now();
  const int nf = static_cast<int>(fixed.size());
  const int nr = static_cast<int>(u_hat.size());
  if (fixed_direction.size() != nf)
    throw std::invalid_argument("fixed curvature direction has wrong length");

  RandomEffectHessianWorkspace<Model> workspace(model, fixed, u_hat,
                                                 partition, false);
  const auto blocks = workspace.EvaluateJointHessianBlocks(fixed, u_hat);
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> factor;
  factor.compute(blocks.random_hessian_m);
  if (factor.info() != Eigen::Success)
    throw std::runtime_error("random Hessian factorization failed");
  const Eigen::VectorXd du =
      -factor.solve(blocks.mixed_hessian_m * fixed_direction);

  std::vector<double> full = merge_parameters(fixed, u_hat, partition);
  std::vector<double> velocity(full.size(), 0.0);
  for (int i = 0; i < nf; ++i)
    velocity[partition.fixed_indices_m[static_cast<size_t>(i)]] =
        fixed_direction[i];
  for (int i = 0; i < nr; ++i)
    velocity[partition.random_indices_m[static_cast<size_t>(i)]] = du[i];
  const auto objective = [&model](const auto &x) { return model.evaluate(x); };

  const Eigen::MatrixXd Huu = Eigen::MatrixXd(blocks.random_hessian_m);
  const double joint_curvature =
      fixed_direction.dot(blocks.fixed_hessian_m * fixed_direction) +
      2.0 * du.dot(blocks.mixed_hessian_m * fixed_direction) +
      du.dot(Huu * du);

  const auto cubic = [&](const std::vector<double> &direction) {
    return had::third_directional_derivative(objective, full, direction);
  };
  Eigen::VectorXd third_u_vv(nr);
  for (int i = 0; i < nr; ++i) {
    std::vector<double> plus = velocity, minus = velocity, unit(full.size(), 0.0);
    const size_t index = partition.random_indices_m[static_cast<size_t>(i)];
    plus[index] += 1.0;
    minus[index] -= 1.0;
    unit[index] = 1.0;
    third_u_vv[i] = (cubic(plus) - cubic(minus) - 2.0 * cubic(unit)) / 6.0;
  }
  const Eigen::VectorXd d2u = -factor.solve(third_u_vv);
  std::vector<double> acceleration(full.size(), 0.0);
  for (int i = 0; i < nr; ++i)
    acceleration[partition.random_indices_m[static_cast<size_t>(i)]] = d2u[i];

  const std::vector<int> random_indices = random_indices_as_ints(partition);
  const Eigen::MatrixXd H1 = dense_hdot_third_order_polarized(
      objective, full, velocity, random_indices);
  const Eigen::MatrixXd third_acceleration =
      dense_hdot_third_order_polarized(objective, full, acceleration,
                                        random_indices);
  Eigen::MatrixXd fourth_velocity = Eigen::MatrixXd::Zero(nr, nr);
  for (int a = 0; a < nr; ++a) {
    for (int b = 0; b <= a; ++b) {
      double polarized = 0.0;
      for (int mask = 0; mask < 16; ++mask) {
        const int signs[4] = {(mask & 1) ? 1 : -1, (mask & 2) ? 1 : -1,
                              (mask & 4) ? 1 : -1, (mask & 8) ? 1 : -1};
        std::vector<double> direction(full.size(), 0.0);
        for (size_t k = 0; k < full.size(); ++k)
          direction[k] = (signs[2] + signs[3]) * velocity[k];
        direction[partition.random_indices_m[static_cast<size_t>(a)]] += signs[0];
        direction[partition.random_indices_m[static_cast<size_t>(b)]] += signs[1];
        const double fourth =
            had::fourth_directional_derivative(objective, full, direction);
        polarized += signs[0] * signs[1] * signs[2] * signs[3] * fourth;
      }
      fourth_velocity(a, b) = fourth_velocity(b, a) = polarized / 384.0;
    }
  }
  const Eigen::MatrixXd H2 = fourth_velocity + third_acceleration;
  const Eigen::MatrixXd Hinv_H1 = factor.solve(H1);
  const double logdet_curvature =
      factor.solve(H2).trace() - (Hinv_H1 * Hinv_H1).trace();

  ExactLaplaceDirectionalCurvatureResult result;
  result.joint_profile_curvature_m = joint_curvature;
  result.logdet_curvature_m = logdet_curvature;
  result.curvature_m = joint_curvature + 0.5 * logdet_curvature;
  result.du_m = du;
  result.d2u_m = d2u;
  result.total_ms_m =
      std::chrono::duration<double, std::milli>(Clock::now() - start).count();
  return result;
}

template <class Model>
ExactLaplaceHessianResult exact_laplace_hessian_fourth_order(
    Model &model, const std::vector<double> &fixed,
    const std::vector<double> &u_hat, const ParameterPartition &partition) {
  using Clock = std::chrono::steady_clock;
  const auto start = Clock::now();
  const int n = static_cast<int>(fixed.size());
  ExactLaplaceHessianResult result;
  result.hessian_m = Eigen::MatrixXd::Zero(n, n);
  for (int i = 0; i < n; ++i) {
    Eigen::VectorXd direction = Eigen::VectorXd::Zero(n);
    direction[i] = 1.0;
    result.hessian_m(i, i) =
        exact_laplace_directional_curvature(model, fixed, u_hat, partition,
                                             direction)
            .curvature_m;
    ++result.directional_evaluations_m;
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      Eigen::VectorXd direction = Eigen::VectorXd::Zero(n);
      direction[i] = direction[j] = 1.0;
      const double pair =
          exact_laplace_directional_curvature(model, fixed, u_hat, partition,
                                               direction)
              .curvature_m;
      result.hessian_m(i, j) = result.hessian_m(j, i) =
          0.5 * (pair - result.hessian_m(i, i) - result.hessian_m(j, j));
      ++result.directional_evaluations_m;
    }
  }
  result.total_ms_m =
      std::chrono::duration<double, std::milli>(Clock::now() - start).count();
  return result;
}

} // namespace laplace
} // namespace quadra
