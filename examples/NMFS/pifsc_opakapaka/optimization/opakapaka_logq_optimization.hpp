#pragma once

#include "../quadra/opakapaka_model.hpp"

#include "../../../../core/optimizer.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace opakapaka_example {

template <class Model>
void polish_single_logq_if_helpful(Model &model,
                                 quadra::ParameterVector &params,
                                 quadra::LaplaceOptions &opts,
                                 quadra::OptResult &fit)
{
constexpr double OPAKAPAKA_POLISH_MIN_MEANINGFUL_STEP = 1.0e-8;
constexpr double OPAKAPAKA_POLISH_MIN_MEANINGFUL_DECREASE = 1.0e-10;
if (fit.par.size() != 1)
{
  return;
}

const std::vector<int> fixed_idx = {0};
std::vector<int> random_idx;
for (std::size_t i = 1; i < params.size(); ++i)
{
  random_idx.push_back(static_cast<int>(i));
}

auto eval_at = [&](double theta,
                   std::vector<double> *out_u_hat = nullptr) -> double
{
  auto tmp = params;
  tmp.params.at(0).value = theta;

  Eigen::VectorXd x(1);
  x[0] = theta;

  had::ADGraph graph;
  auto u_hat = quadra::solve_random_effects_laplace(model, tmp, x, fixed_idx,
                                                    random_idx, graph);

  auto res = quadra::laplace_eval_at_u_star(model, tmp, fixed_idx, random_idx,
                                            x, u_hat, graph, opts);

  if (out_u_hat != nullptr)
  {
    *out_u_hat = u_hat;
  }

  return res.value;
};

const double theta0 = fit.par.at(0);
const double f0 = fit.value;
const double h = std::max(1.0e-5, 1.0e-4 * (1.0 + std::abs(theta0)));

const double fm = eval_at(theta0 - h);
const double fp = eval_at(theta0 + h);

if (!std::isfinite(fm) || !std::isfinite(fp) || !std::isfinite(f0))
{
  return;
}

const double g = (fp - fm) / (2.0 * h);
const double curv = (fp - 2.0 * f0 + fm) / (h * h);

if (!std::isfinite(g) || !std::isfinite(curv) || curv <= 0.0)
{
  return;
}

double step = -g / curv;
if (std::abs(step) < OPAKAPAKA_POLISH_MIN_MEANINGFUL_STEP)
{
  return;
}
const double max_step = 0.05;
if (step > max_step)
  step = max_step;
if (step < -max_step)
  step = -max_step;

if (!std::isfinite(step) || std::abs(step) < 1.0e-12)
{
  return;
}

std::vector<double> polished_u_hat;
const double theta1 = theta0 + step;
const double f1 = eval_at(theta1, &polished_u_hat);

if (!std::isfinite(f1) || f1 >= f0)
{
  std::cout << "Opakapaka log_q polish rejected: " << "step = " << step
            << ", f0 = " << f0 << ", f1 = " << f1 << ", fd_grad = " << g
            << ", fd_curvature = " << curv << "\n";
  return;
}

const double h2 = std::max(1.0e-5, 1.0e-4 * (1.0 + std::abs(theta1)));
const double fm2 = eval_at(theta1 - h2);
const double fp2 = eval_at(theta1 + h2);
double g2 = std::numeric_limits<double>::quiet_NaN();
if (std::isfinite(fm2) && std::isfinite(fp2))
{
  g2 = (fp2 - fm2) / (2.0 * h2);
}

fit.par.at(0) = theta1;
fit.u_hat = polished_u_hat;
fit.value = f1;
if (std::isfinite(g2))
{
  fit.grad_norm = std::abs(g2);
}
fit.converged = true;
fit.message = "accepted safeguarded one-dimensional log_q polish after "
              "line-search stall";

std::cout << "Opakapaka log_q polish accepted: " << "step = " << step
          << ", objective = " << fit.value << ", fd_grad_before = " << g
          << ", fd_curvature = " << curv << ", fd_grad_after = " << g2
          << "\n";
}


template <class Model>
quadra::OptResult fit_log_q_fd_newton_fallback(Model &model,
                                             quadra::ParameterVector &params,
                                             quadra::LaplaceOptions &opts,
                                             double initial_log_q)
{
const std::vector<int> fixed_idx = {0};
std::vector<int> random_idx;
for (std::size_t i = 1; i < params.size(); ++i)
{
  random_idx.push_back(static_cast<int>(i));
}

struct Eval
{
  double value = std::numeric_limits<double>::infinity();
  std::vector<double> u_hat;
};

auto eval_at = [&](double theta) -> Eval
{
  auto tmp = params;
  tmp.params.at(0).value = theta;

  Eigen::VectorXd x(1);
  x[0] = theta;

  had::ADGraph graph;
  Eval out;
  out.u_hat = quadra::solve_random_effects_laplace(model, tmp, x, fixed_idx,
                                                   random_idx, graph);

  auto res = quadra::laplace_eval_at_u_star(model, tmp, fixed_idx, random_idx,
                                            x, out.u_hat, graph, opts);

  out.value = res.value;
  return out;
};

double theta = initial_log_q;
Eval cur = eval_at(theta);
double grad = std::numeric_limits<double>::infinity();
double curv = std::numeric_limits<double>::quiet_NaN();
int iter = 0;

for (; iter < 25; ++iter)
{
  const double h = std::max(1.0e-5, 1.0e-4 * (1.0 + std::abs(theta)));
  const Eval left = eval_at(theta - h);
  const Eval right = eval_at(theta + h);

  if (!std::isfinite(left.value) || !std::isfinite(right.value) ||
      !std::isfinite(cur.value))
  {
    break;
  }

  grad = (right.value - left.value) / (2.0 * h);
  curv = (right.value - 2.0 * cur.value + left.value) / (h * h);

  if (std::abs(grad) < 1.0e-4)
  {
    break;
  }
  if (!std::isfinite(curv) || curv <= 0.0)
  {
    break;
  }

  double step = -grad / curv;
  step = std::max(-1.0, std::min(1.0, step));

  bool accepted = false;
  for (int bt = 0; bt < 20; ++bt)
  {
    const double trial_theta = theta + step;
    Eval trial = eval_at(trial_theta);
    if (std::isfinite(trial.value) && trial.value <= cur.value)
    {
      theta = trial_theta;
      cur = std::move(trial);
      accepted = true;
      break;
    }
    step *= 0.5;
  }

  if (!accepted || std::abs(step) < 1.0e-10)
  {
    break;
  }
}

// One final centered derivative at the returned point.
{
  const double h = std::max(1.0e-5, 1.0e-4 * (1.0 + std::abs(theta)));
  const Eval left = eval_at(theta - h);
  const Eval right = eval_at(theta + h);
  if (std::isfinite(left.value) && std::isfinite(right.value))
  {
    grad = (right.value - left.value) / (2.0 * h);
  }
}

params.params.at(0).value = theta;

quadra::OptResult out;
out.par = std::vector<double>{theta};
out.value = cur.value;
out.grad_norm = std::abs(grad);
out.converged = std::abs(grad) < 1.0e-4;
out.iterations = iter;
out.message = out.converged ? "accepted local safeguarded one-dimensional "
                              "log_q fallback after LBFGS line-search stall"
                            : "local safeguarded one-dimensional log_q "
                              "fallback stopped before requested tolerance";
out.u_hat = cur.u_hat;
return out;
}


}  // namespace opakapaka_example

using opakapaka_example::fit_log_q_fd_newton_fallback;
using opakapaka_example::polish_single_logq_if_helpful;
