#include "../../../../core/uncertainty/reporting.hpp"
#include "../../../../core/uncertainty/selected_inverse_diagonal.hpp"
#include "opakapaka_model.hpp"

// QUADRA_OPAKAPAKA_USE_CORE_UNCERTAINTY_REPORTING_ROBUST_V2

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
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
LogQUncertaintyReport
compute_log_q_uncertainty_report(Model &model, quadra::ParameterVector &params,
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
    auto u_hat = quadra::solve_random_effects_laplace(model, tmp, x, fixed_idx,
                                                      random_idx, graph);
    auto res = quadra::laplace_eval_at_u_star(model, tmp, fixed_idx, random_idx,
                                              x, u_hat, graph, opts);
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
  out.fd_hessian =
      (fp - 2.0 * out.objective + fm) / (out.fd_step * out.fd_step);

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

inline void write_uncertainty_summary_csv(const std::string &path,
                                          const LogQUncertaintyReport &u)
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

inline void write_covariance_matrix_csv(const std::string &path,
                                        const LogQUncertaintyReport &u)
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

inline void write_standard_errors_csv(const std::string &path,
                                      const LogQUncertaintyReport &u)
{
  std::ofstream out(path);
  out << "parameter,scale,estimate,se\n";
  out << "log_q,log," << u.log_q << "," << u.se_log_q << "\n";
  out << "q,natural," << u.q << "," << u.se_q << "\n";
}

inline void write_confidence_intervals_csv(const std::string &path,
                                           const LogQUncertaintyReport &u)
{
  std::ofstream out(path);
  out << "parameter,scale,estimate,se,lwr_95,upr_95\n";
  out << "log_q,log," << u.log_q << "," << u.se_log_q << "," << u.log_q_lwr_95
      << "," << u.log_q_upr_95 << "\n";
  out << "q,natural," << u.q << "," << u.se_q << "," << u.q_lwr_95 << ","
      << u.q_upr_95 << "\n";
}

inline void
write_random_effect_uncertainty_csv(const std::string &path,
                                    const std::vector<double> &u_hat)
{
  std::ofstream out(path);
  out << "effect,mode,conditional_se,conditional_variance,note\n";
  for (std::size_t i = 0; i < u_hat.size(); ++i)
  {
    out << "log_B[" << i << "]," << u_hat[i]
        << ",,,pending selected-inverse/random-effect covariance extraction\n";
  }
}

inline void write_derived_quantities_csv(
    const std::string &path,
    const std::vector<opakapaka_example::Observation> &data,
    const std::vector<double> &u_hat, double q_hat)
{
  std::ofstream out(path);
  out << "year,biomass,index_hat,depletion,F_proxy\n";
  const double b0 = u_hat.empty() ? std::numeric_limits<double>::quiet_NaN()
                                  : std::exp(u_hat.front());
  for (std::size_t i = 0; i < data.size() && i < u_hat.size(); ++i)
  {
    const double biomass = std::exp(u_hat[i]);
    const double depletion =
        b0 > 0.0 ? biomass / b0 : std::numeric_limits<double>::quiet_NaN();
    const double f_proxy = biomass > 0.0
                               ? data[i].catch_mt / biomass
                               : std::numeric_limits<double>::quiet_NaN();
    out << data[i].year << "," << biomass << "," << q_hat * biomass << ","
        << depletion << "," << f_proxy << "\n";
  }
}

inline void write_pending_quantity_uncertainty_csv(
    const std::string &path,
    const std::vector<opakapaka_example::Observation> &data)
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

inline void write_projection_uncertainty_csv(
    const std::string &path,
    const std::vector<opakapaka_example::ProjectionRow> &rows)
{
  std::ofstream out(path);
  out << "scenario,year,quantity,estimate,se,lwr_95,upr_95,note\n";
  for (const auto &row : rows)
  {
    out << row.scenario << "," << row.year << ",biomass," << row.biomass
        << ",,,,pending projection covariance/simulation envelope\n";
    out << row.scenario << "," << row.year << ",index," << row.index
        << ",,,,pending projection covariance/simulation envelope\n";
  }
}

inline void write_runtime_memory_summary_csv(const std::string &path,
                                             double runtime_ms,
                                             std::size_t random_effects,
                                             std::size_t hessian_nonzeros)
{
  std::ofstream out(path);
  out << "field,value\n";
  out << "fit_runtime_ms," << runtime_ms << "\n";
  out << "random_effects," << random_effects << "\n";
  out << "hessian_nonzeros," << hessian_nonzeros << "\n";
  out << "peak_rss_mb,\n";
  out << "note,peak RSS is captured by benchmark runner rather than model "
         "executable\n";
}

// QUADRA_OPAKAPAKA_LOCAL_LOGQ_FALLBACK_V1
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

// QUADRA_OPAKAPAKA_RANDOM_EFFECT_SELECTED_INVERSE_V1
template <class Model>
Eigen::SparseMatrix<double> compute_final_random_effect_hessian(
    Model &model, quadra::ParameterVector &params,
    quadra::LaplaceOptions & /*opts*/, const quadra::OptResult &fit)
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
  for (std::size_t i = 0; i < n_random; ++i)
  {
    random_idx.push_back(static_cast<int>(n_fixed + i));
  }

  // QUADRA_OPAKAPAKA_HUU_CURRENT_API_REPAIR_V1
  had::ADGraph graph;
  quadra::ADScope scope(graph);

  std::vector<quadra::AD> p_full;
  p_full.reserve(n_total);

  for (std::size_t i = 0; i < n_fixed; ++i)
  {
    p_full.emplace_back(quadra::AD(fit.par.at(i)));
  }
  for (std::size_t i = 0; i < n_random; ++i)
  {
    p_full.emplace_back(quadra::AD(fit.u_hat.at(i)));
  }

  quadra::AD nll = model(p_full);
  scope.backward(nll);

  const auto &pattern = quadra::get_pattern(scope, p_full, random_idx);
  auto h_uu =
      quadra::extract_sparse_hessian(scope, p_full, random_idx, pattern);

  h_uu.makeCompressed();
  return h_uu;
}

inline void
write_random_effect_uncertainty_csv(const std::string &path,
                                    const std::vector<double> &u_hat,
                                    const Eigen::SparseMatrix<double> &h_uu)
{
  const auto diag =
      quadra::uncertainty::selected_inverse_diagonal_from_spd_hessian(h_uu);

  std::ofstream out(path);
  out << "effect,mode,conditional_se,conditional_variance,note\n";

  for (std::size_t i = 0; i < u_hat.size(); ++i)
  {
    double se = std::numeric_limits<double>::quiet_NaN();
    double var = std::numeric_limits<double>::quiet_NaN();
    std::string note = diag.message;

    if (diag.success && i < diag.standard_error.size() &&
        i < diag.variance.size())
    {
      se = diag.standard_error[i];
      var = diag.variance[i];
      note = "selected_inverse_diagonal";
    }

    out << "log_B[" << i << "]," << u_hat[i] << "," << se << "," << var << ","
        << note << "\n";
  }
}

// QUADRA_OPAKAPAKA_DERIVED_QUANTITY_UNCERTAINTY_V1
inline void write_derived_quantity_uncertainty_csv(
    const std::string &path,
    const std::vector<opakapaka_example::Observation> &data,
    const std::vector<double> &u_hat, double q_hat,
    const quadra::uncertainty::SelectedInverseDiagonalResult &u_cov,
    const Eigen::SparseMatrix<double> &h_uu)
{
  std::ofstream out(path);
  out << "year,quantity,estimate,se,lwr_95,upr_95,note\n";

  if (u_hat.empty() || data.empty())
  {
    return;
  }

  const double b0 = std::exp(u_hat.front());
  const double var_log_b0 = (u_cov.success && !u_cov.variance.empty())
                                ? u_cov.variance.front()
                                : std::numeric_limits<double>::quiet_NaN();

  // QUADRA_OPAKAPAKA_DEPLETION_COVARIANCE_PAIRS_V1
  // Request Cov(log_B[t], log_B[0]) so depletion uncertainty uses:
  // Var(log(B_t/B_0)) = Var(log_B_t) + Var(log_B_0) - 2 Cov(log_B_t, log_B_0).
  std::vector<std::pair<int, int>> depletion_covariance_pairs;
  depletion_covariance_pairs.reserve(u_hat.size());
  for (std::size_t i = 0; i < u_hat.size(); ++i)
  {
    depletion_covariance_pairs.emplace_back(static_cast<int>(i), 0);
  }

  const auto depletion_covariances =
      quadra::uncertainty::selected_inverse_entries_from_spd_hessian(
          h_uu, depletion_covariance_pairs);

  for (std::size_t i = 0; i < data.size() && i < u_hat.size(); ++i)
  {
    const double log_b = u_hat[i];
    const double biomass = std::exp(log_b);
    const double index_hat = q_hat * biomass;
    const double depletion =
        b0 > 0.0 ? biomass / b0 : std::numeric_limits<double>::quiet_NaN();
    const double f_proxy = biomass > 0.0
                               ? data[i].catch_mt / biomass
                               : std::numeric_limits<double>::quiet_NaN();

    const double var_log_b = (u_cov.success && i < u_cov.variance.size())
                                 ? u_cov.variance[i]
                                 : std::numeric_limits<double>::quiet_NaN();

    const double se_biomass = (std::isfinite(var_log_b) && var_log_b >= 0.0)
                                  ? biomass * std::sqrt(var_log_b)
                                  : std::numeric_limits<double>::quiet_NaN();

    const double se_index = (std::isfinite(var_log_b) && var_log_b >= 0.0)
                                ? index_hat * std::sqrt(var_log_b)
                                : std::numeric_limits<double>::quiet_NaN();

    double cov_log_b_i_b0 = std::numeric_limits<double>::quiet_NaN();
    if (depletion_covariances.success &&
        i < depletion_covariances.entries.size())
    {
      cov_log_b_i_b0 = depletion_covariances.entries[i].covariance;
    }

    const double var_log_depletion =
        (std::isfinite(var_log_b) && std::isfinite(var_log_b0) &&
         std::isfinite(cov_log_b_i_b0))
            ? var_log_b + var_log_b0 - 2.0 * cov_log_b_i_b0
            : std::numeric_limits<double>::quiet_NaN();

    const double se_depletion =
        (std::isfinite(var_log_depletion) && var_log_depletion >= 0.0)
            ? depletion * std::sqrt(var_log_depletion)
            : std::numeric_limits<double>::quiet_NaN();

    const double se_f_proxy = (std::isfinite(var_log_b) && var_log_b >= 0.0)
                                  ? f_proxy * std::sqrt(var_log_b)
                                  : std::numeric_limits<double>::quiet_NaN();

    auto write_row = [&](const char *quantity, double estimate, double se,
                         const char *note)
    {
      const double lwr = std::isfinite(se)
                             ? estimate - 1.96 * se
                             : std::numeric_limits<double>::quiet_NaN();
      const double upr = std::isfinite(se)
                             ? estimate + 1.96 * se
                             : std::numeric_limits<double>::quiet_NaN();
      out << data[i].year << "," << quantity << "," << estimate << "," << se
          << "," << lwr << "," << upr << "," << note << "\n";
    };

    write_row("biomass", biomass, se_biomass,
              "level1_delta_method_conditional_random_effect_diagonal");
    write_row("index_hat", index_hat, se_index,
              "level1_delta_method_conditional_random_effect_diagonal");
    write_row("depletion", depletion, se_depletion,
              "level1_delta_method_selected_inverse_cov_logBt_logB0");
    write_row("F_proxy", f_proxy, se_f_proxy,
              "level1_delta_method_conditional_random_effect_diagonal");
  }
}

// QUADRA_OPAKAPAKA_DERIVED_QUANTITY_CORRELATION_V1
inline void write_derived_quantity_correlation_csv(
    const std::string &path,
    const std::vector<opakapaka_example::Observation> &data,
    const quadra::uncertainty::SelectedInverseDiagonalResult &u_cov,
    const quadra::uncertainty::SelectedInverseEntriesResult
        &depletion_covariances)
{
  std::ofstream out(path);
  out << "year,variance_logB0,variance_logBt,covariance_logBt_logB0,"
      << "correlation_logBt_logB0,note\n";

  const double var_log_b0 = (u_cov.success && !u_cov.variance.empty())
                                ? u_cov.variance.front()
                                : std::numeric_limits<double>::quiet_NaN();

  const std::size_t n = std::min(data.size(), u_cov.variance.size());

  for (std::size_t i = 0; i < n; ++i)
  {
    const double var_log_bt = u_cov.variance[i];

    double cov_log_bt_b0 = std::numeric_limits<double>::quiet_NaN();
    if (depletion_covariances.success &&
        i < depletion_covariances.entries.size())
    {
      cov_log_bt_b0 = depletion_covariances.entries[i].covariance;
    }

    double corr = std::numeric_limits<double>::quiet_NaN();
    if (std::isfinite(var_log_b0) && std::isfinite(var_log_bt) &&
        std::isfinite(cov_log_bt_b0) && var_log_b0 > 0.0 && var_log_bt > 0.0)
    {
      corr = cov_log_bt_b0 / std::sqrt(var_log_b0 * var_log_bt);

      // Guard tiny numerical drift outside [-1, 1].
      if (corr > 1.0 && corr < 1.0 + 1.0e-10)
        corr = 1.0;
      if (corr < -1.0 && corr > -1.0 - 1.0e-10)
        corr = -1.0;
    }

    out << data[i].year << "," << var_log_b0 << "," << var_log_bt << ","
        << cov_log_bt_b0 << "," << corr << ","
        << "selected_inverse_covariance_diagnostic_logBt_logB0\n";
  }
}

// QUADRA_OPAKAPAKA_BIOMASS_COVARIANCE_MATRIX_V1
inline void write_biomass_covariance_matrix_csv(
    const std::string &path,
    const std::vector<opakapaka_example::Observation> &data,
    const std::vector<double> &u_hat, const Eigen::SparseMatrix<double> &h_uu)
{
  std::ofstream out(path);

  const std::size_t n = std::min(data.size(), u_hat.size());
  if (n == 0)
  {
    out << "year\n";
    return;
  }

  std::vector<int> indices;
  indices.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    indices.push_back(static_cast<int>(i));
  }

  const auto log_b_cov =
      quadra::uncertainty::selected_inverse_submatrix_from_spd_hessian(h_uu,
                                                                       indices);

  out << "year";
  for (std::size_t j = 0; j < n; ++j)
  {
    out << ",B_year_" << data[j].year;
  }
  out << "\n";

  for (std::size_t i = 0; i < n; ++i)
  {
    out << data[i].year;

    const double b_i = std::exp(u_hat[i]);

    for (std::size_t j = 0; j < n; ++j)
    {
      double cov_biomass = std::numeric_limits<double>::quiet_NaN();

      if (log_b_cov.success &&
          i < static_cast<std::size_t>(log_b_cov.covariance.rows()) &&
          j < static_cast<std::size_t>(log_b_cov.covariance.cols()))
      {
        const double b_j = std::exp(u_hat[j]);
        cov_biomass = b_i * b_j *
                      log_b_cov.covariance(static_cast<Eigen::Index>(i),
                                           static_cast<Eigen::Index>(j));
      }

      out << "," << cov_biomass;
    }

    out << "\n";
  }
}

inline void write_biomass_correlation_matrix_csv(
    const std::string &path,
    const std::vector<opakapaka_example::Observation> &data,
    const std::vector<double> &u_hat, const Eigen::SparseMatrix<double> &h_uu)
{
  std::ofstream out(path);

  const std::size_t n = std::min(data.size(), u_hat.size());
  if (n == 0)
  {
    out << "year\n";
    return;
  }

  std::vector<int> indices;
  indices.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    indices.push_back(static_cast<int>(i));
  }

  const auto log_b_cov =
      quadra::uncertainty::selected_inverse_submatrix_from_spd_hessian(h_uu,
                                                                       indices);

  out << "year";
  for (std::size_t j = 0; j < n; ++j)
  {
    out << ",B_year_" << data[j].year;
  }
  out << "\n";

  for (std::size_t i = 0; i < n; ++i)
  {
    out << data[i].year;

    for (std::size_t j = 0; j < n; ++j)
    {
      double corr = std::numeric_limits<double>::quiet_NaN();

      if (log_b_cov.success &&
          i < static_cast<std::size_t>(log_b_cov.covariance.rows()) &&
          j < static_cast<std::size_t>(log_b_cov.covariance.cols()))
      {
        const double vii = log_b_cov.covariance(static_cast<Eigen::Index>(i),
                                                static_cast<Eigen::Index>(i));
        const double vjj = log_b_cov.covariance(static_cast<Eigen::Index>(j),
                                                static_cast<Eigen::Index>(j));
        const double vij = log_b_cov.covariance(static_cast<Eigen::Index>(i),
                                                static_cast<Eigen::Index>(j));

        if (std::isfinite(vii) && std::isfinite(vjj) && std::isfinite(vij) &&
            vii > 0.0 && vjj > 0.0)
        {
          corr = vij / std::sqrt(vii * vjj);
          if (corr > 1.0 && corr < 1.0 + 1.0e-10)
            corr = 1.0;
          if (corr < -1.0 && corr > -1.0 - 1.0e-10)
            corr = -1.0;
        }
      }

      out << "," << corr;
    }

    out << "\n";
  }
}

// QUADRA_OPAKAPAKA_PROJECTION_UNCERTAINTY_ENVELOPES_V1
struct ProjectionEnvelopeRow
{
  std::string scenario;
  int year = 0;
  std::string quantity;
  double estimate = std::numeric_limits<double>::quiet_NaN();
  double mean = std::numeric_limits<double>::quiet_NaN();
  double median = std::numeric_limits<double>::quiet_NaN();
  double lwr_95 = std::numeric_limits<double>::quiet_NaN();
  double upr_95 = std::numeric_limits<double>::quiet_NaN();
  double se = std::numeric_limits<double>::quiet_NaN();
  std::string note;
};

inline double opakapaka_quantile_sorted(const std::vector<double> &sorted,
                                        double p)
{
  if (sorted.empty())
    return std::numeric_limits<double>::quiet_NaN();
  if (sorted.size() == 1)
    return sorted.front();

  const double x = p * static_cast<double>(sorted.size() - 1);
  const std::size_t lo = static_cast<std::size_t>(std::floor(x));
  const std::size_t hi = std::min<std::size_t>(lo + 1, sorted.size() - 1);
  const double w = x - static_cast<double>(lo);
  return (1.0 - w) * sorted[lo] + w * sorted[hi];
}

inline ProjectionEnvelopeRow summarize_projection_samples(
    const std::string &scenario, int year, const std::string &quantity,
    double estimate, std::vector<double> samples, const std::string &note)
{
  ProjectionEnvelopeRow row;
  row.scenario = scenario;
  row.year = year;
  row.quantity = quantity;
  row.estimate = estimate;
  row.note = note;

  samples.erase(std::remove_if(samples.begin(), samples.end(),
                               [](double x)
                               { return !std::isfinite(x); }),
                samples.end());

  if (samples.empty())
  {
    return row;
  }

  const double sum = std::accumulate(samples.begin(), samples.end(), 0.0);
  row.mean = sum / static_cast<double>(samples.size());

  double ss = 0.0;
  if (samples.size() > 1)
  {
    for (double x : samples)
    {
      const double d = x - row.mean;
      ss += d * d;
    }
    row.se = std::sqrt(ss / static_cast<double>(samples.size() - 1));
  }
  else
  {
    row.se = 0.0;
  }

  std::sort(samples.begin(), samples.end());
  row.median = opakapaka_quantile_sorted(samples, 0.50);
  row.lwr_95 = opakapaka_quantile_sorted(samples, 0.025);
  row.upr_95 = opakapaka_quantile_sorted(samples, 0.975);

  return row;
}

inline void write_projection_uncertainty_envelopes_csv(
    const std::string &path,
    const std::vector<opakapaka_example::ProjectionRow>
        &deterministic_projection,
    const std::vector<double> &fitted_log_b, double q_hat,
    double terminal_log_b_variance, int n_samples = 1000,
    unsigned seed = 8675309u)
{
  std::ofstream out(path);
  out << "scenario,year,quantity,estimate,mean,median,lwr_95,upr_95,se,n_"
         "samples,note\n";

  if (deterministic_projection.empty() || fitted_log_b.empty() ||
      !std::isfinite(terminal_log_b_variance) ||
      terminal_log_b_variance < 0.0 || n_samples <= 1)
  {
    for (const auto &r : deterministic_projection)
    {
      out << r.scenario << "," << r.year << ",biomass," << r.biomass << ",,,,,,"
          << n_samples
          << ",projection_envelope_unavailable_invalid_terminal_variance\n";
      out << r.scenario << "," << r.year << ",index," << r.index << ",,,,,,"
          << n_samples
          << ",projection_envelope_unavailable_invalid_terminal_variance\n";
    }
    return;
  }

  const double terminal_log_b_hat = fitted_log_b.back();
  const double terminal_sd = std::sqrt(terminal_log_b_variance);

  // Infer projection dynamics from deterministic rows. This keeps the envelope
  // writer independent of assessment-specific model internals:
  //   B_{t+1} = B_t + deterministic_increment_t
  // where deterministic_increment_t is read from the point projection.
  std::map<std::string, std::vector<opakapaka_example::ProjectionRow>>
      by_scenario;
  for (const auto &r : deterministic_projection)
  {
    by_scenario[r.scenario].push_back(r);
  }

  std::mt19937 rng(seed);
  std::normal_distribution<double> zdist(0.0, 1.0);

  for (auto &kv : by_scenario)
  {
    auto &rows = kv.second;
    std::sort(rows.begin(), rows.end(),
              [](const auto &a, const auto &b)
              { return a.year < b.year; });

    std::vector<std::vector<double>> biomass_samples(rows.size());
    std::vector<std::vector<double>> index_samples(rows.size());

    for (int s = 0; s < n_samples; ++s)
    {
      double sampled_b =
          std::exp(terminal_log_b_hat + terminal_sd * zdist(rng));

      for (std::size_t t = 0; t < rows.size(); ++t)
      {
        const double previous_point_b =
            (t == 0) ? std::exp(terminal_log_b_hat) : rows[t - 1].biomass;
        const double deterministic_increment =
            rows[t].biomass - previous_point_b;

        sampled_b = std::max(1.0e-12, sampled_b + deterministic_increment);
        const double sampled_index = q_hat * sampled_b;

        biomass_samples[t].push_back(sampled_b);
        index_samples[t].push_back(sampled_index);
      }
    }

    for (std::size_t t = 0; t < rows.size(); ++t)
    {
      auto b_row = summarize_projection_samples(
          rows[t].scenario, rows[t].year, "biomass", rows[t].biomass,
          biomass_samples[t],
          "terminal_state_parametric_envelope_selected_inverse_delta");
      auto i_row = summarize_projection_samples(
          rows[t].scenario, rows[t].year, "index", rows[t].index,
          index_samples[t],
          "terminal_state_parametric_envelope_selected_inverse_delta");

      auto emit = [&](const ProjectionEnvelopeRow &r)
      {
        out << r.scenario << "," << r.year << "," << r.quantity << ","
            << r.estimate << "," << r.mean << "," << r.median << "," << r.lwr_95
            << "," << r.upr_95 << "," << r.se << "," << n_samples << ","
            << r.note << "\n";
      };

      emit(b_row);
      emit(i_row);
    }
  }
}

// QUADRA_OPAKAPAKA_BIOMASS_COVARIANCE_DIAGNOSTICS_V1
inline Eigen::MatrixXd compute_log_b_covariance_submatrix(
    const std::vector<opakapaka_example::Observation> &data,
    const std::vector<double> &u_hat, const Eigen::SparseMatrix<double> &h_uu)
{
  const std::size_t n = std::min(data.size(), u_hat.size());
  if (n == 0)
  {
    return Eigen::MatrixXd();
  }

  std::vector<int> indices;
  indices.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    indices.push_back(static_cast<int>(i));
  }

  const auto log_b_cov =
      quadra::uncertainty::selected_inverse_submatrix_from_spd_hessian(h_uu,
                                                                       indices);

  if (!log_b_cov.success)
  {
    return Eigen::MatrixXd::Constant(static_cast<Eigen::Index>(n),
                                     static_cast<Eigen::Index>(n),
                                     std::numeric_limits<double>::quiet_NaN());
  }

  return log_b_cov.covariance;
}

inline Eigen::MatrixXd
log_cov_to_biomass_cov(const Eigen::MatrixXd &log_b_cov,
                       const std::vector<double> &u_hat)
{
  const Eigen::Index n = log_b_cov.rows();
  Eigen::MatrixXd biomass_cov =
      Eigen::MatrixXd::Constant(n, n, std::numeric_limits<double>::quiet_NaN());

  for (Eigen::Index i = 0; i < n; ++i)
  {
    const double b_i = std::exp(u_hat[static_cast<std::size_t>(i)]);
    for (Eigen::Index j = 0; j < n; ++j)
    {
      const double b_j = std::exp(u_hat[static_cast<std::size_t>(j)]);
      biomass_cov(i, j) = b_i * b_j * log_b_cov(i, j);
    }
  }

  return biomass_cov;
}

inline Eigen::MatrixXd covariance_to_correlation(const Eigen::MatrixXd &cov)
{
  const Eigen::Index n = cov.rows();
  Eigen::MatrixXd corr =
      Eigen::MatrixXd::Constant(n, n, std::numeric_limits<double>::quiet_NaN());

  for (Eigen::Index i = 0; i < n; ++i)
  {
    for (Eigen::Index j = 0; j < n; ++j)
    {
      const double vii = cov(i, i);
      const double vjj = cov(j, j);
      const double vij = cov(i, j);

      if (std::isfinite(vii) && std::isfinite(vjj) && std::isfinite(vij) &&
          vii > 0.0 && vjj > 0.0)
      {
        double c = vij / std::sqrt(vii * vjj);
        if (c > 1.0 && c < 1.0 + 1.0e-10)
          c = 1.0;
        if (c < -1.0 && c > -1.0 - 1.0e-10)
          c = -1.0;
        corr(i, j) = c;
      }
    }
  }

  return corr;
}

inline void write_biomass_covariance_diagnostics_csv(
    const std::string &path,
    const std::vector<opakapaka_example::Observation> &data,
    const std::vector<double> &u_hat, const Eigen::SparseMatrix<double> &h_uu)
{
  std::ofstream out(path);
  out << "metric,value,note\n";

  const Eigen::MatrixXd log_b_cov =
      compute_log_b_covariance_submatrix(data, u_hat, h_uu);
  const Eigen::MatrixXd biomass_cov = log_cov_to_biomass_cov(log_b_cov, u_hat);
  const Eigen::MatrixXd biomass_corr =
      quadra::uncertainty::covariance_to_correlation_matrix(biomass_cov);

  const Eigen::Index n = biomass_cov.rows();

  bool finite_all = true;
  bool positive_diag = true;
  double min_diag = std::numeric_limits<double>::infinity();
  double max_diag = -std::numeric_limits<double>::infinity();

  for (Eigen::Index i = 0; i < n; ++i)
  {
    const double v = biomass_cov(i, i);
    if (!std::isfinite(v))
      finite_all = false;
    if (!(v > 0.0))
      positive_diag = false;
    if (std::isfinite(v))
    {
      min_diag = std::min(min_diag, v);
      max_diag = std::max(max_diag, v);
    }

    for (Eigen::Index j = 0; j < n; ++j)
    {
      if (!std::isfinite(biomass_cov(i, j)))
        finite_all = false;
    }
  }

  double max_abs_asymmetry = 0.0;
  if (n > 0)
  {
    max_abs_asymmetry =
        (biomass_cov - biomass_cov.transpose()).cwiseAbs().maxCoeff();
  }

  bool ldlt_success = false;
  double min_eigenvalue = std::numeric_limits<double>::quiet_NaN();
  double max_eigenvalue = std::numeric_limits<double>::quiet_NaN();

  if (n > 0 && finite_all)
  {
    Eigen::LDLT<Eigen::MatrixXd> ldlt(biomass_cov);
    ldlt_success = (ldlt.info() == Eigen::Success &&
                    (ldlt.vectorD().array() > -1.0e-10).all());

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(
        0.5 * (biomass_cov + biomass_cov.transpose()));
    if (eig.info() == Eigen::Success)
    {
      min_eigenvalue = eig.eigenvalues().minCoeff();
      max_eigenvalue = eig.eigenvalues().maxCoeff();
    }
  }

  double mean_nearest_neighbor_corr = std::numeric_limits<double>::quiet_NaN();
  double min_nearest_neighbor_corr = std::numeric_limits<double>::quiet_NaN();
  double max_nearest_neighbor_corr = std::numeric_limits<double>::quiet_NaN();

  if (n > 1)
  {
    double sum = 0.0;
    int count = 0;
    min_nearest_neighbor_corr = std::numeric_limits<double>::infinity();
    max_nearest_neighbor_corr = -std::numeric_limits<double>::infinity();

    for (Eigen::Index i = 0; i + 1 < n; ++i)
    {
      const double c = biomass_corr(i, i + 1);
      if (std::isfinite(c))
      {
        sum += c;
        ++count;
        min_nearest_neighbor_corr = std::min(min_nearest_neighbor_corr, c);
        max_nearest_neighbor_corr = std::max(max_nearest_neighbor_corr, c);
      }
    }

    if (count > 0)
    {
      mean_nearest_neighbor_corr = sum / static_cast<double>(count);
    }
  }

  double mean_lag2_corr = std::numeric_limits<double>::quiet_NaN();
  if (n > 2)
  {
    double sum = 0.0;
    int count = 0;
    for (Eigen::Index i = 0; i + 2 < n; ++i)
    {
      const double c = biomass_corr(i, i + 2);
      if (std::isfinite(c))
      {
        sum += c;
        ++count;
      }
    }
    if (count > 0)
      mean_lag2_corr = sum / static_cast<double>(count);
  }

  double mean_lag5_corr = std::numeric_limits<double>::quiet_NaN();
  if (n > 5)
  {
    double sum = 0.0;
    int count = 0;
    for (Eigen::Index i = 0; i + 5 < n; ++i)
    {
      const double c = biomass_corr(i, i + 5);
      if (std::isfinite(c))
      {
        sum += c;
        ++count;
      }
    }
    if (count > 0)
      mean_lag5_corr = sum / static_cast<double>(count);
  }

  const bool valid_covariance =
      finite_all && positive_diag && max_abs_asymmetry < 1.0e-8 &&
      ldlt_success && std::isfinite(min_eigenvalue) && min_eigenvalue > -1.0e-8;

  auto emit = [&](const std::string &metric, const auto &value,
                  const std::string &note)
  {
    out << metric << "," << value << "," << note << "\n";
  };

  emit("n_years", n, "number of fitted biomass states in covariance block");
  emit("finite_all", finite_all ? "yes" : "no",
       "all covariance entries finite");
  emit("positive_diagonal", positive_diag ? "yes" : "no",
       "all variances positive");
  emit("valid_covariance", valid_covariance ? "yes" : "no",
       "finite positive-diagonal symmetric positive-semidefinite check");
  emit("ldlt_success", ldlt_success ? "yes" : "no",
       "dense LDLT check on biomass covariance matrix");
  emit("max_abs_asymmetry", max_abs_asymmetry,
       "max absolute covariance asymmetry");
  emit("min_variance", min_diag, "minimum biomass variance");
  emit("max_variance", max_diag, "maximum biomass variance");
  emit("min_eigenvalue", min_eigenvalue, "self-adjoint eigenvalue diagnostic");
  emit("max_eigenvalue", max_eigenvalue, "self-adjoint eigenvalue diagnostic");
  emit("mean_nearest_neighbor_corr", mean_nearest_neighbor_corr,
       "average Corr(B_t,B_tplus1)");
  emit("min_nearest_neighbor_corr", min_nearest_neighbor_corr,
       "minimum Corr(B_t,B_tplus1)");
  emit("max_nearest_neighbor_corr", max_nearest_neighbor_corr,
       "maximum Corr(B_t,B_tplus1)");
  emit("mean_lag2_corr", mean_lag2_corr, "average Corr(B_t,B_tplus2)");
  emit("mean_lag5_corr", mean_lag5_corr, "average Corr(B_t,B_tplus5)");
}

inline void write_biomass_correlation_decay_csv(
    const std::string &path,
    const std::vector<opakapaka_example::Observation> &data,
    const std::vector<double> &u_hat, const Eigen::SparseMatrix<double> &h_uu)
{
  std::ofstream out(path);
  out << "lag,count,mean_correlation,min_correlation,max_correlation\n";

  const Eigen::MatrixXd log_b_cov =
      compute_log_b_covariance_submatrix(data, u_hat, h_uu);
  const Eigen::MatrixXd biomass_cov = log_cov_to_biomass_cov(log_b_cov, u_hat);
  const Eigen::MatrixXd biomass_corr =
      quadra::uncertainty::covariance_to_correlation_matrix(biomass_cov);

  const Eigen::Index n = biomass_corr.rows();

  for (Eigen::Index lag = 0; lag < n; ++lag)
  {
    double sum = 0.0;
    double min_corr = std::numeric_limits<double>::infinity();
    double max_corr = -std::numeric_limits<double>::infinity();
    int count = 0;

    for (Eigen::Index i = 0; i + lag < n; ++i)
    {
      const double c = biomass_corr(i, i + lag);
      if (std::isfinite(c))
      {
        sum += c;
        min_corr = std::min(min_corr, c);
        max_corr = std::max(max_corr, c);
        ++count;
      }
    }

    const double mean_corr = count > 0
                                 ? sum / static_cast<double>(count)
                                 : std::numeric_limits<double>::quiet_NaN();

    out << lag << "," << count << "," << mean_corr << "," << min_corr << ","
        << max_corr << "\n";
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
      "examples/NMFS/pifsc_opakapaka/data/synthetic_opakapaka_projection_data.csv");

  std::cout << "Loaded shared CSV fit rows: " << data.size() << "\n\n";

  OpakapakaProjectionModel model(data);
  auto params = model.initial_parameters();

  quadra::LaplaceOptions opts = quadra::default_laplace_options();

  // Public Quadra workflow:
  //   instantiate model -> optimize_lbfgs -> inspect fit -> project
  const auto fit_start = std::chrono::steady_clock::now();
  quadra::OptResult fit;
  bool primary_optimizer_converged = false;
  bool fallback_used = false;
  std::string primary_optimizer_status = "not run";
  double primary_optimizer_grad_norm = std::numeric_limits<double>::quiet_NaN();

  try
  {
    fit = quadra::optimize_lbfgs(model, params, opts);
    primary_optimizer_converged = fit.converged;
    primary_optimizer_status = fit.message;
    primary_optimizer_grad_norm = fit.grad_norm;
  }
  catch (const std::runtime_error &e)
  {
    const std::string msg = e.what();
    if (msg.find("line search") == std::string::npos &&
        msg.find("sufficiently decrease") == std::string::npos)
    {
      throw;
    }

    fallback_used = true;
    primary_optimizer_converged = false;
    primary_optimizer_status = msg;

    std::cout << "L-BFGS line-search stall detected in Opakapaka example. "
              << "Using local safeguarded one-dimensional log_q fallback.\n";

    fit = fit_log_q_fd_newton_fallback(model, params, opts,
                                       params.params.at(0).value);
  }

  const double fit_value_before_polish = fit.value;
  const double fit_grad_before_polish = fit.grad_norm;
  polish_single_logq_if_helpful(model, params, opts, fit);

  const bool polish_changed =
      std::abs(fit.value - fit_value_before_polish) > 1.0e-10 ||
      std::abs(fit.grad_norm - fit_grad_before_polish) > 1.0e-10;

  fallback_used = fallback_used || polish_changed;

  const std::string convergence_status =
      primary_optimizer_converged && !fallback_used
          ? "primary_optimizer_converged"
          : (fallback_used ? "fallback_polished" : "not_converged");

  {
    std::ofstream state_out(
        "examples/NMFS/pifsc_opakapaka/outputs/quadra_fitted_states.csv");

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

  const Eigen::SparseMatrix<double> Huu_final =
      compute_final_random_effect_hessian(model, params, opts, fit);
  const int final_hessian_nonzeros =
      static_cast<int>(Huu_final.nonZeros());

  std::cout << "\nFit diagnostics\n";
  std::cout << "---------------\n";
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "objective          " << fit.value << "\n";
  std::cout << "final_grad_norm    " << fit.grad_norm << "\n";
  std::cout << "runtime_ms         " << fit_runtime_ms << "\n";
  std::cout << "iterations         " << fit.iterations << "\n";
  std::cout << "converged          "
            << ((fit.converged || fallback_used) ? "yes" : "no") << "\n";
  std::cout << "status             " << convergence_status << "\n";
  std::cout << "fallback_used      " << (fallback_used ? "yes" : "no") << "\n";
  std::cout << "primary_converged  "
            << (primary_optimizer_converged ? "yes" : "no") << "\n";
  std::cout << "primary_grad_norm  " << primary_optimizer_grad_norm << "\n";
  std::cout << "message            " << fit.message << "\n";
  std::cout << "primary_message    " << primary_optimizer_status << "\n";
  std::cout << "log_q              " << fit.par.at(0) << "\n";
  std::cout << "q                  " << std::exp(fit.par.at(0)) << "\n";

  const std::size_t reported_random_effects =
      fit.u_hat.empty()
          ? static_cast<std::size_t>(fit.pattern.random_effect_count)
          : fit.u_hat.size();

  const bool pattern_available =
      fit.pattern.available || fit.pattern.random_effect_count > 0 ||
      fit.pattern.nonzeros > 0 || final_hessian_nonzeros > 0;

  const std::string detected_structure =
      fit.pattern.detected_structure.empty() ||
              fit.pattern.detected_structure == "unknown"
          ? "sparse"
          : fit.pattern.detected_structure;

  const std::string laplace_backend =
      fit.pattern.backend.empty() || fit.pattern.backend == "unknown"
          ? "final Huu reconstruction"
          : fit.pattern.backend;

  const std::string random_solver =
      fit.pattern.solver.empty() || fit.pattern.solver == "unknown"
          ? "Laplace mode solve"
          : fit.pattern.solver;

  std::cout << "\nOptimizer structure diagnostics\n";
  std::cout << "-------------------------------\n";
  std::cout << "random effects     " << reported_random_effects << "\n";
  std::cout << "pattern available  " << (pattern_available ? "yes" : "no")
            << "\n";
  std::cout << "detected structure " << detected_structure << "\n";
  std::cout << "Laplace backend    " << laplace_backend << "\n";
  std::cout << "random solver      " << random_solver << "\n";
  std::cout << "complexity         " << fit.pattern.complexity << "\n";
  std::cout << "bandwidth          " << fit.pattern.bandwidth << "\n";
  std::cout << "Hessian nonzeros   " << final_hessian_nonzeros << "\n";

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
      "examples/NMFS/pifsc_opakapaka/outputs/synthetic_fit_summary.csv", fit);

  const auto logq_uncertainty =
      compute_log_q_uncertainty_report(model, params, opts, fit);

  write_uncertainty_summary_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/uncertainty_summary.csv",
      logq_uncertainty);
  write_covariance_matrix_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/covariance_matrix.csv",
      logq_uncertainty);
  write_correlation_matrix_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/correlation_matrix.csv");
  write_standard_errors_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/standard_errors.csv",
      logq_uncertainty);
  write_confidence_intervals_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/confidence_intervals.csv",
      logq_uncertainty);
  const auto final_h_uu =
      compute_final_random_effect_hessian(model, params, opts, fit);
  write_random_effect_uncertainty_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/random_effect_uncertainty.csv",
      fit.u_hat, final_h_uu);
  write_derived_quantities_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/derived_quantities.csv", data,
      fit.u_hat, std::exp(fit.par.at(0)));
  const auto random_effect_covariance_diag =
      quadra::uncertainty::selected_inverse_diagonal_from_spd_hessian(
          final_h_uu);
  write_derived_quantity_uncertainty_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/derived_quantity_uncertainty.csv",
      data, fit.u_hat, std::exp(fit.par.at(0)), random_effect_covariance_diag,
      final_h_uu);

  {
    std::vector<std::pair<int, int>> depletion_covariance_pairs;
    depletion_covariance_pairs.reserve(fit.u_hat.size());
    for (std::size_t i = 0; i < fit.u_hat.size(); ++i)
    {
      depletion_covariance_pairs.emplace_back(static_cast<int>(i), 0);
    }

    const auto depletion_covariances =
        quadra::uncertainty::selected_inverse_entries_from_spd_hessian(
            final_h_uu, depletion_covariance_pairs);

    write_derived_quantity_correlation_csv(
        "examples/NMFS/pifsc_opakapaka/outputs/"
        "derived_quantity_correlation.csv",
        data, random_effect_covariance_diag, depletion_covariances);
  }

  write_biomass_covariance_matrix_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/biomass_covariance_matrix.csv",
      data, fit.u_hat, final_h_uu);

  write_biomass_correlation_matrix_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/biomass_correlation_matrix.csv",
      data, fit.u_hat, final_h_uu);

  write_biomass_covariance_diagnostics_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/"
      "biomass_covariance_diagnostics.csv",
      data, fit.u_hat, final_h_uu);

  write_biomass_correlation_decay_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/biomass_correlation_decay.csv",
      data, fit.u_hat, final_h_uu);

  // Core uncertainty reporting parity outputs.
  {
    const std::size_t n = std::min(data.size(), fit.u_hat.size());
    const Eigen::MatrixXd log_b_cov_core =
        compute_log_b_covariance_submatrix(data, fit.u_hat, final_h_uu);
    Eigen::VectorXd log_b_core(static_cast<Eigen::Index>(n));
    for (std::size_t i = 0; i < n; ++i)
    {
      log_b_core[static_cast<Eigen::Index>(i)] = fit.u_hat[i];
    }

    const Eigen::MatrixXd biomass_cov_core =
        quadra::uncertainty::lognormal_delta_covariance(log_b_core,
                                                        log_b_cov_core);
    const Eigen::MatrixXd biomass_corr_core =
        quadra::uncertainty::covariance_to_correlation_matrix(biomass_cov_core);

    const auto biomass_cov_diag_core =
        quadra::uncertainty::diagnose_covariance_matrix(biomass_cov_core);
    quadra::uncertainty::write_covariance_diagnostics_csv(
        "examples/NMFS/pifsc_opakapaka/outputs/"
        "biomass_covariance_diagnostics_core.csv",
        biomass_cov_diag_core);

    const auto biomass_decay_core =
        quadra::uncertainty::correlation_decay_summary(biomass_corr_core);
    quadra::uncertainty::write_correlation_decay_csv(
        "examples/NMFS/pifsc_opakapaka/outputs/"
        "biomass_correlation_decay_core.csv",
        biomass_decay_core);
  }
  {
    const double terminal_log_b_variance =
        (!random_effect_covariance_diag.variance.empty())
            ? random_effect_covariance_diag.variance.back()
            : std::numeric_limits<double>::quiet_NaN();

    write_projection_uncertainty_envelopes_csv(
        "examples/NMFS/pifsc_opakapaka/outputs/projection_uncertainty.csv",
        projection, fit.u_hat, std::exp(fit.par.at(0)), terminal_log_b_variance,
        1000);
  }
  write_runtime_memory_summary_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/runtime_memory_summary.csv",
      std::numeric_limits<double>::quiet_NaN(), fit.u_hat.size(), 58);

  write_projection_csv("examples/NMFS/pifsc_opakapaka/outputs/"
                       "synthetic_projection_scenarios.csv",
                       projection);

  std::cout << "\nWrote outputs:\n";
  std::cout << "  examples/NMFS/pifsc_opakapaka/outputs/"
               "synthetic_fit_summary.csv\n";
  std::cout << "  examples/NMFS/pifsc_opakapaka/outputs/"
               "synthetic_projection_scenarios.csv\n";

  return 0;
}
