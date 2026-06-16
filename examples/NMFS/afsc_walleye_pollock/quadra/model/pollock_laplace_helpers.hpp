#pragma once

#include "pollock_model.hpp"

#include "../../../../../core/optimizer.hpp"
#include "../../../../../core/laplace/laplace_structure_report.hpp"
#include "../../../../../core/laplace/functional_analysis_report.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace pollock_example {

std::vector<int> pollock_fixed_indices(const quadra::ParameterVector &params) {
std::vector<int> out;
for (std::size_t i = 0; i < params.params.size(); ++i) {
if (!params.params[i].is_random) {
  out.push_back(static_cast<int>(i));
}
}
return out;
}


std::vector<int> pollock_random_indices(const quadra::ParameterVector &params) {
std::vector<int> out;
for (std::size_t i = 0; i < params.params.size(); ++i) {
if (params.params[i].is_random) {
  out.push_back(static_cast<int>(i));
}
}
return out;
}


double pollock_joint_objective_at_x_u(
PollockModel &model,
quadra::ParameterVector params,
const std::vector<int> &fixed_idx,
const std::vector<int> &random_idx,
const Eigen::VectorXd &x,
const std::vector<double> &u) {
for (std::size_t k = 0; k < fixed_idx.size(); ++k) {
params.params[static_cast<std::size_t>(fixed_idx[k])].value =
    x[static_cast<Eigen::Index>(k)];
}
for (std::size_t k = 0; k < random_idx.size(); ++k) {
params.params[static_cast<std::size_t>(random_idx[k])].value = u[k];
}

std::vector<double> p_double;
p_double.reserve(static_cast<std::size_t>(params.size()));
for (int i = 0; i < params.size(); ++i) {
p_double.emplace_back(params.params[static_cast<std::size_t>(i)].value);
}

return model(p_double);
}


Eigen::VectorXd pollock_fixed_values(const quadra::ParameterVector &params) {
const auto fixed_idx = pollock_fixed_indices(params);
Eigen::VectorXd x(static_cast<Eigen::Index>(fixed_idx.size()));
for (std::size_t i = 0; i < fixed_idx.size(); ++i) {
x(static_cast<Eigen::Index>(i)) =
    params.params[static_cast<std::size_t>(fixed_idx[i])].value;
}
return x;
}


std::vector<double> pollock_profile_gradient_fd_at_x(
PollockModel &model,
quadra::ParameterVector params,
const Eigen::VectorXd &x,
const std::vector<double> &u_hat,
const quadra::LaplaceOptions &opts,
double fd_step = 1.0e-5) {
const auto fixed_idx = pollock_fixed_indices(params);
const auto random_idx = pollock_random_indices(params);

std::vector<double> grad(static_cast<std::size_t>(x.size()), 0.0);

auto eval = [&](const Eigen::VectorXd &x_eval) -> double {
had::ADGraph graph;
auto res = quadra::laplace_eval_at_u_star(
    model, params, fixed_idx, random_idx, x_eval, u_hat, graph, opts);
return res.value;
};

for (Eigen::Index j = 0; j < x.size(); ++j) {
const double h = fd_step * std::max(1.0, std::abs(x(j)));
Eigen::VectorXd xp = x;
Eigen::VectorXd xm = x;
xp(j) += h;
xm(j) -= h;
grad[static_cast<std::size_t>(j)] = (eval(xp) - eval(xm)) / (2.0 * h);
}

return grad;
}


Eigen::MatrixXd pollock_fd_huu(
PollockModel &model,
quadra::ParameterVector params,
const quadra::OptResult &fit,
double eps = 1.0e-4) {
const auto fixed_idx = quadra::build_fixed_index(params);
const auto random_idx = quadra::build_random_index(params);

const Eigen::Index n = static_cast<Eigen::Index>(random_idx.size());
Eigen::MatrixXd H = Eigen::MatrixXd::Zero(n, n);

if (n == 0 || fit.par.size() != fixed_idx.size() ||
  fit.u_hat.size() != random_idx.size()) {
return H;
}

Eigen::VectorXd x(static_cast<Eigen::Index>(fixed_idx.size()));
for (std::size_t i = 0; i < fixed_idx.size(); ++i) {
x[static_cast<Eigen::Index>(i)] = fit.par[i];
}

std::vector<double> u = fit.u_hat;
const double f0 = pollock_joint_objective_at_x_u(
  model, params, fixed_idx, random_idx, x, u);

for (Eigen::Index i = 0; i < n; ++i) {
std::vector<double> up = u;
std::vector<double> um = u;
up[static_cast<std::size_t>(i)] += eps;
um[static_cast<std::size_t>(i)] -= eps;

const double fp = pollock_joint_objective_at_x_u(
    model, params, fixed_idx, random_idx, x, up);
const double fm = pollock_joint_objective_at_x_u(
    model, params, fixed_idx, random_idx, x, um);

H(i, i) = (fp - 2.0 * f0 + fm) / (eps * eps);

for (Eigen::Index j = i + 1; j < n; ++j) {
  std::vector<double> upp = u;
  std::vector<double> upm = u;
  std::vector<double> ump = u;
  std::vector<double> umm = u;

  upp[static_cast<std::size_t>(i)] += eps;
  upp[static_cast<std::size_t>(j)] += eps;

  upm[static_cast<std::size_t>(i)] += eps;
  upm[static_cast<std::size_t>(j)] -= eps;

  ump[static_cast<std::size_t>(i)] -= eps;
  ump[static_cast<std::size_t>(j)] += eps;

  umm[static_cast<std::size_t>(i)] -= eps;
  umm[static_cast<std::size_t>(j)] -= eps;

  const double fpp = pollock_joint_objective_at_x_u(
      model, params, fixed_idx, random_idx, x, upp);
  const double fpm = pollock_joint_objective_at_x_u(
      model, params, fixed_idx, random_idx, x, upm);
  const double fmp = pollock_joint_objective_at_x_u(
      model, params, fixed_idx, random_idx, x, ump);
  const double fmm = pollock_joint_objective_at_x_u(
      model, params, fixed_idx, random_idx, x, umm);

  const double hij = (fpp - fpm - fmp + fmm) / (4.0 * eps * eps);
  H(i, j) = hij;
  H(j, i) = hij;
}
}

return H;
}


}  // namespace pollock_example

// Compatibility aliases for existing Pollock diagnostics/driver call sites.
using pollock_example::pollock_fixed_indices;
using pollock_example::pollock_random_indices;
using pollock_example::pollock_joint_objective_at_x_u;
using pollock_example::pollock_fixed_values;
using pollock_example::pollock_profile_gradient_fd_at_x;
using pollock_example::pollock_fd_huu;
