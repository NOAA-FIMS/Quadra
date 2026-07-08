#pragma once

#include "../../../../../core/optimizer.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace pifsc_bigeye_tuna {

struct SafeLaplaceValue
{
  bool ok = false;
  double value = std::numeric_limits<double>::quiet_NaN();
  std::string message;
};

template <class Objective>
SafeLaplaceValue safe_profiled_laplace_value_at_fixed(
    Objective &objective,
    quadra::ParameterVector params,
    const std::vector<double> &theta,
    quadra::LaplaceOptions opts)
{
  SafeLaplaceValue out;

  try
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

    out.ok = std::isfinite(res.value);
    out.value = res.value;
    out.message = out.ok ? "ok" : "nonfinite value";
  }
  catch (const std::exception &e)
  {
    out.ok = false;
    out.message = e.what();
  }
  catch (...)
  {
    out.ok = false;
    out.message = "unknown exception";
  }

  return out;
}

inline std::vector<std::string> safe_fixed_effect_names(
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

struct SafeWiggleRow
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
  std::string status = "not_run";
  std::string message;
};

struct SafeWiggleDiagnostics
{
  double relative_step = 0.02;
  double absolute_min_step = 1.0e-3;
  std::vector<SafeWiggleRow> rows;
};

template <class Objective>
SafeWiggleDiagnostics make_safe_fixed_effect_wiggle_diagnostics(
    Objective &objective,
    const quadra::ParameterVector &params,
    const quadra::OptResult &fit,
    quadra::LaplaceOptions opts,
    double relative_step = 0.02,
    double absolute_min_step = 1.0e-3)
{
  SafeWiggleDiagnostics out;
  out.relative_step = relative_step;
  out.absolute_min_step = absolute_min_step;

  const auto names = safe_fixed_effect_names(params, fit);
  const auto base = safe_profiled_laplace_value_at_fixed(
      objective, params, fit.par, opts);

  for (std::size_t i = 0; i < fit.par.size(); ++i)
  {
    SafeWiggleRow row;
    row.index = i;
    row.name = i < names.size() ? names[i] : ("fixed_" + std::to_string(i));
    row.base_value = fit.par[i];
    row.step = std::max(absolute_min_step,
                        relative_step * (1.0 + std::abs(fit.par[i])));
    row.f_base = base.value;

    if (!base.ok)
    {
      row.status = "base_failed";
      row.message = base.message;
      out.rows.push_back(row);
      continue;
    }

    auto minus = fit.par;
    auto plus = fit.par;
    minus[i] -= row.step;
    plus[i] += row.step;

    const auto fm = safe_profiled_laplace_value_at_fixed(
        objective, params, minus, opts);
    const auto fp = safe_profiled_laplace_value_at_fixed(
        objective, params, plus, opts);

    row.f_minus = fm.value;
    row.f_plus = fp.value;

    if (fm.ok && fp.ok)
    {
      row.delta_minus = row.f_minus - row.f_base;
      row.delta_plus = row.f_plus - row.f_base;
      row.central_gradient = (row.f_plus - row.f_minus) / (2.0 * row.step);
      row.central_curvature =
          (row.f_plus - 2.0 * row.f_base + row.f_minus) / (row.step * row.step);
      row.abs_sensitivity =
          std::max(std::abs(row.delta_minus), std::abs(row.delta_plus));
      row.status = "ok";
      row.message = "ok";
    }
    else
    {
      row.status = "perturbation_failed";
      row.message = "minus=" + fm.message + "; plus=" + fp.message;
    }

    out.rows.push_back(row);
  }

  std::sort(out.rows.begin(), out.rows.end(),
            [](const SafeWiggleRow &a, const SafeWiggleRow &b)
            {
              const bool a_finite = std::isfinite(a.abs_sensitivity);
              const bool b_finite = std::isfinite(b.abs_sensitivity);
              if (a_finite != b_finite)
                return a_finite;
              return a.abs_sensitivity < b.abs_sensitivity;
            });

  return out;
}

inline void write_safe_fixed_effect_wiggle_diagnostics_csv(
    const SafeWiggleDiagnostics &diagnostics,
    const std::string &path)
{
  std::ofstream out(path);
  if (!out)
    throw std::runtime_error("Could not open safe wiggle diagnostics CSV: " +
                             path);

  out << "index,name,base_value,step,f_base,f_minus,f_plus,"
      << "delta_minus,delta_plus,central_gradient,central_curvature,"
      << "abs_sensitivity,status,message\n";
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
        << row.abs_sensitivity << ","
        << row.status << ","
        << '"' << row.message << '"' << "\n";
  }
}

inline void write_safe_fixed_effect_wiggle_diagnostics_text(
    const SafeWiggleDiagnostics &diagnostics,
    const std::string &path)
{
  std::ofstream out(path);
  if (!out)
    throw std::runtime_error("Could not open safe wiggle diagnostics text: " +
                             path);

  out << "Safe Fixed Effect Wiggle Diagnostics\n";
  out << "====================================\n\n";
  out << std::setprecision(15);
  out << "relative_step:       " << diagnostics.relative_step << "\n";
  out << "absolute_min_step:   " << diagnostics.absolute_min_step << "\n\n";

  out << "Interpretation\n";
  out << "--------------\n";
  out << "Perturbations that fail are retained as diagnostic rows instead of aborting\n";
  out << "the model run. Failed perturbations identify regions where the profiled\n";
  out << "Laplace solve or Huu factorization is not stable.\n\n";

  out << "Rows sorted by finite abs_sensitivity, failed rows last\n";
  out << "------------------------------------------------------\n";
  out << "index,name,base_value,step,f_base,f_minus,f_plus,"
      << "delta_minus,delta_plus,central_gradient,central_curvature,"
      << "abs_sensitivity,status,message\n";

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
        << row.abs_sensitivity << ","
        << row.status << ","
        << row.message << "\n";
  }
}

template <class Objective>
void write_safe_fixed_effect_wiggle_diagnostics(
    const std::string &text_path,
    const std::string &csv_path,
    Objective &objective,
    const quadra::ParameterVector &params,
    const quadra::OptResult &fit,
    quadra::LaplaceOptions opts,
    double relative_step = 0.02,
    double absolute_min_step = 1.0e-3)
{
  const auto diagnostics = make_safe_fixed_effect_wiggle_diagnostics(
      objective, params, fit, opts, relative_step, absolute_min_step);

  write_safe_fixed_effect_wiggle_diagnostics_text(diagnostics, text_path);
  write_safe_fixed_effect_wiggle_diagnostics_csv(diagnostics, csv_path);
}

} // namespace pifsc_bigeye_tuna

using pifsc_bigeye_tuna::SafeWiggleDiagnostics;
using pifsc_bigeye_tuna::write_safe_fixed_effect_wiggle_diagnostics;
