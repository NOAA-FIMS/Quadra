#pragma once

#include "../../../../../core/optimizer.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace pifsc_bigeye_tuna {

template <class Objective>
double bigeye_wiggle_profiled_laplace_value_at_fixed(
    Objective &objective,
    quadra::ParameterVector params,
    const std::vector<double> &theta,
    quadra::LaplaceOptions opts)
{
  for (std::size_t i = 0; i < theta.size() && i < params.params.size(); ++i)
  {
    params.params[i].value = theta[i];
  }

  std::vector<int> fixed_idx;
  std::vector<int> random_idx;

  for (std::size_t i = 0; i < params.params.size(); ++i)
  {
    if (params.params[i].is_random)
      random_idx.push_back(static_cast<int>(i));
    else
      fixed_idx.push_back(static_cast<int>(i));
  }

  Eigen::VectorXd x(static_cast<Eigen::Index>(fixed_idx.size()));
  for (std::size_t i = 0; i < fixed_idx.size(); ++i)
  {
    x[static_cast<Eigen::Index>(i)] = params.params[fixed_idx[i]].value;
  }

  had::ADGraph graph;
  const auto u_hat = quadra::solve_random_effects_laplace(
      objective, params, x, fixed_idx, random_idx, graph);

  const auto res = quadra::laplace_eval_at_u_star(
      objective, params, fixed_idx, random_idx, x, u_hat, graph, opts);

  return res.value;
}

inline std::vector<std::string> bigeye_wiggle_fixed_effect_names(
    const quadra::ParameterVector &params,
    const quadra::OptResult &fit)
{
  if (fit.fixed_gradient_names.size() == fit.par.size())
    return fit.fixed_gradient_names;

  std::vector<std::string> names;
  for (const auto &p : params.params)
  {
    if (!p.is_random)
      names.push_back(p.name);
  }

  while (names.size() < fit.par.size())
    names.push_back("fixed_" + std::to_string(names.size()));

  return names;
}

struct BigeyeWiggleRow
{
  std::size_t index = 0;
  std::string name;
  double base_value = std::numeric_limits<double>::quiet_NaN();
  double step = std::numeric_limits<double>::quiet_NaN();
  double f_base = std::numeric_limits<double>::quiet_NaN();
  double f_minus = std::numeric_limits<double>::quiet_NaN();
  double f_plus = std::numeric_limits<double>::quiet_NaN();
  double delta_minus = std::numeric_limits<double>::quiet_NaN();
  double delta_plus = std::numeric_limits<double>::quiet_NaN();
  double central_gradient = std::numeric_limits<double>::quiet_NaN();
  double central_curvature = std::numeric_limits<double>::quiet_NaN();
  double abs_sensitivity = std::numeric_limits<double>::quiet_NaN();
};

struct BigeyeWiggleDiagnostics
{
  double relative_step = 0.05;
  double absolute_min_step = 1.0e-3;
  std::vector<BigeyeWiggleRow> rows;
};

template <class Objective>
BigeyeWiggleDiagnostics make_bigeye_fixed_effect_wiggle_diagnostics(
    Objective &objective,
    const quadra::ParameterVector &params,
    const quadra::OptResult &fit,
    quadra::LaplaceOptions opts,
    double relative_step = 0.05,
    double absolute_min_step = 1.0e-3)
{
  BigeyeWiggleDiagnostics out;
  out.relative_step = relative_step;
  out.absolute_min_step = absolute_min_step;

  const auto names = bigeye_wiggle_fixed_effect_names(params, fit);

  if (fit.par.empty())
    return out;

  const double f0 = bigeye_wiggle_profiled_laplace_value_at_fixed(
      objective, params, fit.par, opts);

  for (std::size_t i = 0; i < fit.par.size(); ++i)
  {
    BigeyeWiggleRow row;
    row.index = i;
    row.name = i < names.size() ? names[i] : ("fixed_" + std::to_string(i));
    row.base_value = fit.par[i];
    row.step = std::max(absolute_min_step,
                        relative_step * (1.0 + std::abs(fit.par[i])));
    row.f_base = f0;

    auto minus = fit.par;
    auto plus = fit.par;
    minus[i] -= row.step;
    plus[i] += row.step;

    row.f_minus = bigeye_wiggle_profiled_laplace_value_at_fixed(
        objective, params, minus, opts);
    row.f_plus = bigeye_wiggle_profiled_laplace_value_at_fixed(
        objective, params, plus, opts);

    row.delta_minus = row.f_minus - f0;
    row.delta_plus = row.f_plus - f0;
    row.central_gradient = (row.f_plus - row.f_minus) / (2.0 * row.step);
    row.central_curvature =
        (row.f_plus - 2.0 * f0 + row.f_minus) / (row.step * row.step);
    row.abs_sensitivity =
        std::max(std::abs(row.delta_minus), std::abs(row.delta_plus));

    out.rows.push_back(row);
  }

  std::sort(out.rows.begin(), out.rows.end(),
            [](const BigeyeWiggleRow &a, const BigeyeWiggleRow &b)
            {
              return a.abs_sensitivity < b.abs_sensitivity;
            });

  return out;
}

inline void write_bigeye_fixed_effect_wiggle_diagnostics_csv(
    const BigeyeWiggleDiagnostics &diagnostics,
    const std::string &path)
{
  std::ofstream out(path);
  if (!out)
    throw std::runtime_error("Could not open wiggle diagnostics CSV: " + path);

  out << "index,name,base_value,step,f_base,f_minus,f_plus,"
      << "delta_minus,delta_plus,central_gradient,central_curvature,"
      << "abs_sensitivity\n";

  out << std::setprecision(15);
  for (const auto &row : diagnostics.rows)
  {
    out << row.index << ","
        << row.name << ","
        << row.base_value << ","
        << row.step << ","
        << row.f_base << ","
        << row.f_minus << ","
        << row.f_plus << ","
        << row.delta_minus << ","
        << row.delta_plus << ","
        << row.central_gradient << ","
        << row.central_curvature << ","
        << row.abs_sensitivity << "\n";
  }
}

inline void write_bigeye_fixed_effect_wiggle_diagnostics_text(
    const BigeyeWiggleDiagnostics &diagnostics,
    const std::string &path)
{
  std::ofstream out(path);
  if (!out)
    throw std::runtime_error("Could not open wiggle diagnostics text: " + path);

  out << "Fixed Effect Wiggle Diagnostics\n";
  out << "===============================\n\n";
  out << std::setprecision(15);
  out << "relative_step:       " << diagnostics.relative_step << "\n";
  out << "absolute_min_step:   " << diagnostics.absolute_min_step << "\n\n";

  out << "Interpretation\n";
  out << "--------------\n";
  out << "Each fixed effect is perturbed plus/minus one step while random effects are\n";
  out << "re-solved through the Laplace machinery. Small objective changes indicate a\n";
  out << "locally weak or compensable parameter. Large asymmetric changes may indicate\n";
  out << "nonlinearity, boundary behavior, or a poor local quadratic approximation.\n\n";

  out << "Rows sorted by abs_sensitivity ascending\n";
  out << "----------------------------------------\n";
  out << "index,name,base_value,step,f_base,f_minus,f_plus,"
      << "delta_minus,delta_plus,central_gradient,central_curvature,"
      << "abs_sensitivity\n";

  for (const auto &row : diagnostics.rows)
  {
    out << row.index << ","
        << row.name << ","
        << row.base_value << ","
        << row.step << ","
        << row.f_base << ","
        << row.f_minus << ","
        << row.f_plus << ","
        << row.delta_minus << ","
        << row.delta_plus << ","
        << row.central_gradient << ","
        << row.central_curvature << ","
        << row.abs_sensitivity << "\n";
  }
}

template <class Objective>
void write_bigeye_fixed_effect_wiggle_diagnostics(
    const std::string &text_path,
    const std::string &csv_path,
    Objective &objective,
    const quadra::ParameterVector &params,
    const quadra::OptResult &fit,
    quadra::LaplaceOptions opts,
    double relative_step = 0.05,
    double absolute_min_step = 1.0e-3)
{
  const auto diagnostics =
      make_bigeye_fixed_effect_wiggle_diagnostics(
          objective, params, fit, opts, relative_step, absolute_min_step);

  write_bigeye_fixed_effect_wiggle_diagnostics_text(diagnostics, text_path);
  write_bigeye_fixed_effect_wiggle_diagnostics_csv(diagnostics, csv_path);
}

} // namespace pifsc_bigeye_tuna

using pifsc_bigeye_tuna::BigeyeWiggleDiagnostics;
using pifsc_bigeye_tuna::write_bigeye_fixed_effect_wiggle_diagnostics;
