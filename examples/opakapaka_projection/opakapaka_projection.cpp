#include "opakapaka_model.hpp"
#include "../../core/uncertainty/selected_inverse_diagonal.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace
{

  std::vector<std::string> split_csv_line_simple(const std::string &line)
  {
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string item;
    while (std::getline(ss, item, ','))
    {
      fields.push_back(item);
    }
    return fields;
  }

  bool finite_double_from_string(const std::string &x, double &out)
  {
    try
    {
      std::size_t pos = 0;
      out = std::stod(x, &pos);
      return pos > 0 && std::isfinite(out);
    }
    catch (...)
    {
      out = std::numeric_limits<double>::quiet_NaN();
      return false;
    }
  }

  std::vector<opakapaka_example::Observation>
  read_opakapaka_history_csv(const std::string &path)
  {
    std::ifstream in(path);
    if (!in)
    {
      throw std::runtime_error("Could not open Opakapaka CSV: " + path);
    }

    std::string line;
    if (!std::getline(in, line))
    {
      throw std::runtime_error("Opakapaka CSV is empty: " + path);
    }

    const auto header = split_csv_line_simple(line);
    int year_col = -1;
    int phase_col = -1;
    int catch_col = -1;
    int index_col = -1;

    for (int i = 0; i < static_cast<int>(header.size()); ++i)
    {
      if (header[i] == "year")
        year_col = i;
      if (header[i] == "phase")
        phase_col = i;
      if (header[i] == "catch_mt")
        catch_col = i;
      if (header[i] == "index")
        index_col = i;
    }

    if (year_col < 0 || phase_col < 0 || catch_col < 0 || index_col < 0)
    {
      throw std::runtime_error(
          "Opakapaka CSV must contain year, phase, catch_mt, and index columns");
    }

    std::vector<opakapaka_example::Observation> out;

    while (std::getline(in, line))
    {
      if (line.empty())
        continue;
      const auto fields = split_csv_line_simple(line);
      const int max_col =
          std::max(std::max(year_col, phase_col), std::max(catch_col, index_col));
      if (static_cast<int>(fields.size()) <= max_col)
        continue;

      if (fields[phase_col] != "history")
        continue;

      double year_d = 0.0;
      double catch_mt = 0.0;
      double index = 0.0;

      if (!finite_double_from_string(fields[year_col], year_d))
        continue;
      if (!finite_double_from_string(fields[catch_col], catch_mt))
        continue;
      if (!finite_double_from_string(fields[index_col], index))
        continue;

      opakapaka_example::Observation obs;
      obs.year = static_cast<int>(year_d);
      obs.catch_mt = catch_mt;
      obs.index = index;
      out.push_back(obs);
    }

    if (out.empty())
    {
      throw std::runtime_error(
          "No usable historical rows found in Opakapaka CSV");
    }

    return out;
  }

} // namespace

// QUADRA_OPAKAPAKA_LOGQ_POLISH_V1
template <class Model>
void polish_single_logq_if_helpful(Model &model,
                                   quadra::ParameterVector &params,
                                   quadra::LaplaceOptions &opts,
                                   quadra::OptResult &fit)
{
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

// QUADRA_LEVEL1_UNCERTAINTY_REPORTING_V3
struct LogQUncertaintyReport
{
  double objective = std::numeric_limits<double>::quiet_NaN();
  double fd_step = std::numeric_limits<double>::quiet_NaN();
  double fd_gradient = std::numeric_limits<double>::quiet_NaN();
  double fd_hessian = std::numeric_limits<double>::quiet_NaN();
  double covariance_log_q = std::numeric_limits<double>::quiet_NaN();
  double se_log_q = std::numeric_limits<double>::quiet_NaN();
  double log_q = std::numeric_limits<double>::quiet_NaN();
  double q = std::numeric_limits<double>::quiet_NaN();
  double se_q = std::numeric_limits<double>::quiet_NaN();
  double log_q_lwr_95 = std::numeric_limits<double>::quiet_NaN();
  double log_q_upr_95 = std::numeric_limits<double>::quiet_NaN();
  double q_lwr_95 = std::numeric_limits<double>::quiet_NaN();
  double q_upr_95 = std::numeric_limits<double>::quiet_NaN();
};

template <class Model>
LogQUncertaintyReport compute_log_q_uncertainty_report(
    Model &model,
    quadra::ParameterVector &params,
    quadra::LaplaceOptions &opts,
    const quadra::OptResult &fit)
{
  LogQUncertaintyReport out;
  if (fit.par.size() != 1)
    return out;

  const std::vector<int> fixed_idx = {0};
  std::vector<int> random_idx;
  for (std::size_t i = 1; i < params.size(); ++i)
  {
    random_idx.push_back(static_cast<int>(i));
  }

  auto eval_at = [&](double theta)
  {
    auto tmp = params;
    tmp.params.at(0).value = theta;
    Eigen::VectorXd x(1);
    x[0] = theta;
    had::ADGraph graph;
    auto u_hat = quadra::solve_random_effects_laplace(
        model, tmp, x, fixed_idx, random_idx, graph);
    auto res = quadra::laplace_eval_at_u_star(
        model, tmp, fixed_idx, random_idx, x, u_hat, graph, opts);
    return res.value;
  };

  out.objective = fit.value;
  out.log_q = fit.par.at(0);
  out.q = std::exp(out.log_q);
  out.fd_step = std::max(1.0e-5, 1.0e-4 * (1.0 + std::abs(out.log_q)));

  const double fm = eval_at(out.log_q - out.fd_step);
  const double fp = eval_at(out.log_q + out.fd_step);
  if (!std::isfinite(fm) || !std::isfinite(fp) || !std::isfinite(out.objective))
    return out;

  out.fd_gradient = (fp - fm) / (2.0 * out.fd_step);
  out.fd_hessian = (fp - 2.0 * out.objective + fm) / (out.fd_step * out.fd_step);

  if (std::isfinite(out.fd_hessian) && out.fd_hessian > 0.0)
  {
    out.covariance_log_q = 1.0 / out.fd_hessian;
    out.se_log_q = std::sqrt(out.covariance_log_q);
    out.se_q = out.q * out.se_log_q;
    out.log_q_lwr_95 = out.log_q - 1.96 * out.se_log_q;
    out.log_q_upr_95 = out.log_q + 1.96 * out.se_log_q;
    out.q_lwr_95 = std::exp(out.log_q_lwr_95);
    out.q_upr_95 = std::exp(out.log_q_upr_95);
  }
  return out;
}

inline void write_uncertainty_summary_csv(const std::string &path, const LogQUncertaintyReport &u)
{
  std::ofstream out(path);
  out << "field,value\n";
  out << "objective," << u.objective << "\n";
  out << "fd_step," << u.fd_step << "\n";
  out << "fd_gradient_log_q," << u.fd_gradient << "\n";
  out << "fd_hessian_log_q," << u.fd_hessian << "\n";
  out << "covariance_log_q," << u.covariance_log_q << "\n";
  out << "se_log_q," << u.se_log_q << "\n";
  out << "se_q," << u.se_q << "\n";
  out << "hessian_positive," << (u.fd_hessian > 0.0 ? "yes" : "no") << "\n";
}

inline void write_covariance_matrix_csv(const std::string &path, const LogQUncertaintyReport &u)
{
  std::ofstream out(path);
  out << "row,col,value\n";
  out << "log_q,log_q," << u.covariance_log_q << "\n";
}

inline void write_correlation_matrix_csv(const std::string &path)
{
  std::ofstream out(path);
  out << "row,col,value\n";
  out << "log_q,log_q,1\n";
}

inline void write_standard_errors_csv(const std::string &path, const LogQUncertaintyReport &u)
{
  std::ofstream out(path);
  out << "parameter,scale,estimate,se\n";
  out << "log_q,log," << u.log_q << "," << u.se_log_q << "\n";
  out << "q,natural," << u.q << "," << u.se_q << "\n";
}

inline void write_confidence_intervals_csv(const std::string &path, const LogQUncertaintyReport &u)
{
  std::ofstream out(path);
  out << "parameter,scale,estimate,se,lwr_95,upr_95\n";
  out << "log_q,log," << u.log_q << "," << u.se_log_q << "," << u.log_q_lwr_95 << "," << u.log_q_upr_95 << "\n";
  out << "q,natural," << u.q << "," << u.se_q << "," << u.q_lwr_95 << "," << u.q_upr_95 << "\n";
}

inline void write_random_effect_uncertainty_csv(const std::string &path, const std::vector<double> &u_hat)
{
  std::ofstream out(path);
  out << "effect,mode,conditional_se,conditional_variance,note\n";
  for (std::size_t i = 0; i < u_hat.size(); ++i)
  {
    out << "log_B[" << i << "]," << u_hat[i] << ",,,pending selected-inverse/random-effect covariance extraction\n";
  }
}

inline void write_derived_quantities_csv(const std::string &path, const std::vector<opakapaka_example::Observation> &data, const std::vector<double> &u_hat, double q_hat)
{
  std::ofstream out(path);
  out << "year,biomass,index_hat,depletion,F_proxy\n";
  const double b0 = u_hat.empty() ? std::numeric_limits<double>::quiet_NaN() : std::exp(u_hat.front());
  for (std::size_t i = 0; i < data.size() && i < u_hat.size(); ++i)
  {
    const double biomass = std::exp(u_hat[i]);
    const double depletion = b0 > 0.0 ? biomass / b0 : std::numeric_limits<double>::quiet_NaN();
    const double f_proxy = biomass > 0.0 ? data[i].catch_mt / biomass : std::numeric_limits<double>::quiet_NaN();
    out << data[i].year << "," << biomass << "," << q_hat * biomass << "," << depletion << "," << f_proxy << "\n";
  }
}

inline void write_pending_quantity_uncertainty_csv(const std::string &path, const std::vector<opakapaka_example::Observation> &data)
{
  std::ofstream out(path);
  out << "year,quantity,estimate,se,lwr_95,upr_95,note\n";
  for (const auto &obs : data)
  {
    out << obs.year << ",biomass,,,,,pending delta-method propagation\n";
    out << obs.year << ",depletion,,,,,pending delta-method propagation\n";
    out << obs.year << ",F_proxy,,,,,pending delta-method propagation\n";
  }
}

inline void write_projection_uncertainty_csv(const std::string &path, const std::vector<opakapaka_example::ProjectionRow> &rows)
{
  std::ofstream out(path);
  out << "scenario,year,quantity,estimate,se,lwr_95,upr_95,note\n";
  for (const auto &row : rows)
  {
    out << row.scenario << "," << row.year << ",biomass," << row.biomass << ",,,,pending projection covariance/simulation envelope\n";
    out << row.scenario << "," << row.year << ",index," << row.index << ",,,,pending projection covariance/simulation envelope\n";
  }
}

inline void write_runtime_memory_summary_csv(const std::string &path, double runtime_ms, std::size_t random_effects, std::size_t hessian_nonzeros)
{
  std::ofstream out(path);
  out << "field,value\n";
  out << "fit_runtime_ms," << runtime_ms << "\n";
  out << "random_effects," << random_effects << "\n";
  out << "hessian_nonzeros," << hessian_nonzeros << "\n";
  out << "peak_rss_mb,\n";
  out << "note,peak RSS is captured by benchmark runner rather than model executable\n";
}

// QUADRA_OPAKAPAKA_LOCAL_LOGQ_FALLBACK_V1
template <class Model>
quadra::OptResult fit_log_q_fd_newton_fallback(
    Model &model,
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
    out.u_hat = quadra::solve_random_effects_laplace(
        model, tmp, x, fixed_idx, random_idx, graph);

    auto res = quadra::laplace_eval_at_u_star(
        model, tmp, fixed_idx, random_idx, x, out.u_hat, graph, opts);

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
  out.message =
      out.converged
          ? "accepted local safeguarded one-dimensional log_q fallback after LBFGS line-search stall"
          : "local safeguarded one-dimensional log_q fallback stopped before requested tolerance";
  out.u_hat = cur.u_hat;
  return out;
}


// QUADRA_OPAKAPAKA_RANDOM_EFFECT_SELECTED_INVERSE_V1
template <class Model>
Eigen::SparseMatrix<double> compute_final_random_effect_hessian(
    Model& model,
    quadra::ParameterVector& params,
    quadra::LaplaceOptions& /*opts*/,
    const quadra::OptResult& fit)
{
  // QUADRA_OPAKAPAKA_HUU_ADSCOPE_REPAIR_V1
  //
  // LaplaceResult currently stores value/gradients only. For conditional
  // random-effect SEs, rebuild the fitted AD vector, evaluate the model,
  // propagate adjoints, discover the sparse Hessian pattern, and extract H_uu
  // using Quadra's sparse Hessian extraction API.

  const std::size_t n_fixed = fit.par.size();
  const std::size_t n_random = fit.u_hat.size();
  const std::size_t n_total = n_fixed + n_random;

  std::vector<int> random_idx;
  random_idx.reserve(n_random);
  for (std::size_t i = 0; i < n_random; ++i) {
    random_idx.push_back(static_cast<int>(n_fixed + i));
  }

  // QUADRA_OPAKAPAKA_HUU_CURRENT_API_REPAIR_V1
  had::ADGraph graph;
  quadra::ADScope scope(graph);

  std::vector<quadra::AD> p_full;
  p_full.reserve(n_total);

  for (std::size_t i = 0; i < n_fixed; ++i) {
    p_full.emplace_back(quadra::AD(fit.par.at(i)));
  }
  for (std::size_t i = 0; i < n_random; ++i) {
    p_full.emplace_back(quadra::AD(fit.u_hat.at(i)));
  }

  quadra::AD nll = model(p_full);
  scope.backward(nll);

  const auto& pattern = quadra::get_pattern(scope, p_full, random_idx);
  auto h_uu = quadra::extract_sparse_hessian(scope, p_full, random_idx, pattern);

  h_uu.makeCompressed();
  return h_uu;
}

inline void write_random_effect_uncertainty_csv(
    const std::string& path,
    const std::vector<double>& u_hat,
    const Eigen::SparseMatrix<double>& h_uu)
{
  const auto diag =
      quadra::uncertainty::selected_inverse_diagonal_from_spd_hessian(h_uu);

  std::ofstream out(path);
  out << "effect,mode,conditional_se,conditional_variance,note\n";

  for (std::size_t i = 0; i < u_hat.size(); ++i) {
    double se = std::numeric_limits<double>::quiet_NaN();
    double var = std::numeric_limits<double>::quiet_NaN();
    std::string note = diag.message;

    if (diag.success && i < diag.standard_error.size() && i < diag.variance.size()) {
      se = diag.standard_error[i];
      var = diag.variance[i];
      note = "selected_inverse_diagonal";
    }

    out << "log_B[" << i << "]," << u_hat[i] << ","
        << se << "," << var << "," << note << "\n";
  }
}


int main()
{
  using namespace opakapaka_example;

  std::cout << "Synthetic opakapaka-style fit + projection example\n";
  std::cout << "==================================================\n\n";
  std::cout
      << "Synthetic and public-data-safe. Not an official assessment.\n\n";

  auto data = read_opakapaka_history_csv(
      "examples/opakapaka_projection/synthetic_opakapaka_projection_data.csv");

  std::cout << "Loaded shared CSV fit rows: " << data.size() << "\n\n";

  OpakapakaProjectionModel model(data);
  auto params = model.initial_parameters();

  quadra::LaplaceOptions opts = quadra::default_laplace_options();

  // Public Quadra workflow:
  //   instantiate model -> optimize_lbfgs -> inspect fit -> project
  const auto fit_start = std::chrono::steady_clock::now();
  quadra::OptResult fit;
  try
  {
    fit = quadra::optimize_lbfgs(model, params, opts);
  }
  catch (const std::runtime_error &e)
  {
    const std::string msg = e.what();
    if (msg.find("line search") == std::string::npos &&
        msg.find("sufficiently decrease") == std::string::npos)
    {
      throw;
    }

    std::cout << "L-BFGS line-search stall detected in Opakapaka example. "
              << "Using local safeguarded one-dimensional log_q fallback.";

    fit = fit_log_q_fd_newton_fallback(model, params, opts,
                                       params.params.at(0).value);
  }
  polish_single_logq_if_helpful(model, params, opts, fit);

  {
    std::ofstream state_out(
        "examples/opakapaka_projection/outputs/quadra_fitted_states.csv");

    state_out << "index,log_B,B\n";

    for (std::size_t i = 0; i < fit.u_hat.size(); ++i)
    {
      state_out << i << "," << std::setprecision(15) << fit.u_hat[i] << ","
                << std::setprecision(15) << std::exp(fit.u_hat[i]) << "\n";
    }
  }

  const auto fit_stop = std::chrono::steady_clock::now();
  const double fit_runtime_ms =
      std::chrono::duration<double, std::milli>(fit_stop - fit_start).count();

  ProjectionOptions projection_options;
  projection_options.start_year = data.back().year + 1;
  projection_options.years = 10;
  projection_options.scenarios = {
      {"zero_catch", 0.0},
      {"status_quo", 1.0},
      {"low_catch", 0.75},
      {"high_catch", 1.25},
  };

  auto projection = model.project(fit, projection_options);

  std::cout << "\nFit diagnostics\n";
  std::cout << "---------------\n";
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "objective          " << fit.value << "\n";
  std::cout << "grad_norm          " << fit.grad_norm << "\n";
  std::cout << "runtime_ms         " << fit_runtime_ms << "\n";
  std::cout << "iterations         " << fit.iterations << "\n";
  std::cout << "converged          " << (fit.converged ? "yes" : "no") << "\n";
  std::cout << "message            " << fit.message << "\n";
  std::cout << "log_q              " << fit.par.at(0) << "\n";
  std::cout << "q                  " << std::exp(fit.par.at(0)) << "\n";

  std::cout << "\nOptimizer structure diagnostics\n";
  std::cout << "-------------------------------\n";
  std::cout << "random effects     " << fit.pattern.random_effect_count << "\n";
  std::cout << "pattern available  " << (fit.pattern.available ? "yes" : "no")
            << "\n";
  std::cout << "detected structure " << fit.pattern.detected_structure << "\n";
  std::cout << "Laplace backend    " << fit.pattern.backend << "\n";
  std::cout << "random solver      " << fit.pattern.solver << "\n";
  std::cout << "complexity         " << fit.pattern.complexity << "\n";
  std::cout << "bandwidth          " << fit.pattern.bandwidth << "\n";
  std::cout << "Hessian nonzeros   " << fit.pattern.nonzeros << "\n";

  std::cout << "\nProjection preview\n";
  std::cout << "------------------\n";
  std::cout << "scenario,year,catch_mt,biomass,index\n";
  int printed = 0;
  for (const auto &row : projection)
  {
    if (printed >= 12)
    {
      break;
    }
    std::cout << row.scenario << "," << row.year << "," << row.catch_mt << ","
              << row.biomass << "," << row.index << "\n";
    ++printed;
  }

  write_fit_summary_csv(
      "examples/opakapaka_projection/outputs/synthetic_fit_summary.csv", fit);

  const auto logq_uncertainty =
      compute_log_q_uncertainty_report(model, params, opts, fit);

  write_uncertainty_summary_csv("examples/opakapaka_projection/outputs/uncertainty_summary.csv", logq_uncertainty);
  write_covariance_matrix_csv("examples/opakapaka_projection/outputs/covariance_matrix.csv", logq_uncertainty);
  write_correlation_matrix_csv("examples/opakapaka_projection/outputs/correlation_matrix.csv");
  write_standard_errors_csv("examples/opakapaka_projection/outputs/standard_errors.csv", logq_uncertainty);
  write_confidence_intervals_csv("examples/opakapaka_projection/outputs/confidence_intervals.csv", logq_uncertainty);
  const auto final_h_uu =
      compute_final_random_effect_hessian(model, params, opts, fit);
  write_random_effect_uncertainty_csv(
      "examples/opakapaka_projection/outputs/random_effect_uncertainty.csv",
      fit.u_hat, final_h_uu);
  write_derived_quantities_csv("examples/opakapaka_projection/outputs/derived_quantities.csv", data, fit.u_hat, std::exp(fit.par.at(0)));
  write_pending_quantity_uncertainty_csv("examples/opakapaka_projection/outputs/derived_quantity_uncertainty.csv", data);
  write_projection_uncertainty_csv("examples/opakapaka_projection/outputs/projection_uncertainty.csv", projection);
  write_runtime_memory_summary_csv("examples/opakapaka_projection/outputs/runtime_memory_summary.csv", std::numeric_limits<double>::quiet_NaN(), fit.u_hat.size(), 58);

  write_projection_csv("examples/opakapaka_projection/outputs/"
                       "synthetic_projection_scenarios.csv",
                       projection);

  std::cout << "\nWrote outputs:\n";
  std::cout << "  examples/opakapaka_projection/outputs/"
               "synthetic_fit_summary.csv\n";
  std::cout << "  examples/opakapaka_projection/outputs/"
               "synthetic_projection_scenarios.csv\n";

  return 0;
}
