#ifndef OPTIMIZER_HPP
#define OPTIMIZER_HPP
#pragma once

#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "../external/LBFGSpp/include/LBFGS.h"
#include "../external/eigen/Eigen/Dense"

#include "autodiff.hpp"
#include "laplace.hpp"
#include "laplace/model_analysis_report.hpp"
#include "laplace/persistent_structured_runtime.hpp"

namespace quadra {

struct OptPatternInfo {
  bool available = false;

  std::string detected_structure = "unknown";
  std::string backend = "unknown";
  std::string solver = "unknown";
  std::string complexity = "unknown";

  int bandwidth = -1;
  std::size_t rows = 0;
  std::size_t cols = 0;
  std::size_t nonzeros = 0;
  std::size_t random_effect_count = 0;
};

struct OptResult {
  // Backward-compatible fixed-effect estimate.
  std::vector<double> par;

  // Random-effect mode at the final fixed-effect estimate.
  std::vector<double> u_hat;

  // Parameter indices used to construct par and u_hat.
  std::vector<int> fixed_index;
  std::vector<int> random_index;

  // Objective and outer-gradient diagnostics.
  double value = std::numeric_limits<double>::quiet_NaN();
  int iterations = 0;
  double grad_norm = std::numeric_limits<double>::quiet_NaN();

  bool converged = false;
  std::string message;

  // Random-effect Hessian / backend diagnostic payload.
  //
  // v1 fills the random-effect count and leaves detailed structure as unknown.
  // The next patch should wire this to the structure detector / backend
  // factory.
  OptPatternInfo pattern;
};

inline Eigen::VectorXd to_eigen(const std::vector<double> &x) {
  Eigen::VectorXd out(static_cast<Eigen::Index>(x.size()));
  for (Eigen::Index i = 0; i < out.size(); ++i) {
    out[i] = x[static_cast<size_t>(i)];
  }
  return out;
}

inline bool all_finite_eigen(const Eigen::VectorXd &v) {
  for (Eigen::Index i = 0; i < v.size(); ++i) {
    if (!std::isfinite(v[i])) {
      return false;
    }
  }
  return true;
}

inline double safe_eigen_norm(const Eigen::VectorXd &v) {
  if (!all_finite_eigen(v)) {
    return std::numeric_limits<double>::infinity();
  }
  return v.norm();
}

inline OptPatternInfo
make_opt_pattern_info_from_report(const laplace::ModelAnalysisReport &report) {
  OptPatternInfo info;
  info.available = true;
  info.detected_structure = laplace::ToString(report.structure);
  info.backend = laplace::ToString(report.backend);
  info.solver = laplace::ToString(report.solver);
  info.complexity = report.complexity;
  info.bandwidth = report.bandwidth;
  info.rows = static_cast<std::size_t>(report.rows);
  info.cols = static_cast<std::size_t>(report.cols);
  info.nonzeros = static_cast<std::size_t>(report.nnz);
  info.random_effect_count = static_cast<std::size_t>(report.random_effects);
  return info;
}

template <typename Model>
OptPatternInfo analyze_final_random_effect_pattern(
    Model &model, ParameterVector &params, const Eigen::VectorXd &x,
    const std::vector<double> &u_hat, const std::vector<int> &fixed_idx,
    const std::vector<int> &random_idx,
    const LaplaceOptions & /*options*/ = default_laplace_options()) {
  OptPatternInfo info;
  info.random_effect_count = random_idx.size();

  if (random_idx.empty() || u_hat.empty()) {
    info.available = false;
    info.detected_structure = "none";
    info.backend = "none";
    info.solver = "none";
    info.complexity = "none";
    return info;
  }

  if (u_hat.size() != random_idx.size()) {
    info.available = false;
    info.detected_structure = "unavailable";
    info.backend = "unavailable";
    info.solver = "unavailable";
    info.complexity = "random-effect mode size mismatch";
    return info;
  }

  try {
    had::ADGraph graph;
    ADScope scope(graph);

    std::vector<AD> p_full;
    p_full.reserve(params.size());

    for (int i = 0; i < params.size(); ++i) {
      p_full.emplace_back(AD(0.0));
    }

    inject_fixed_params(x, p_full, fixed_idx);
    inject_random_params(u_hat, p_full, random_idx);

    AD nll = model(p_full);
    scope.backward(nll);

    const auto &pattern = get_pattern(scope, p_full, random_idx);

    Eigen::SparseMatrix<double> H =
        extract_sparse_hessian(scope, p_full, random_idx, pattern, 0.0);

    laplace::StructureDetectorOptions detector_options;
    detector_options.prefer_dense_for_small_matrices = false;
    detector_options.dense_size_cutoff = 0;

    const laplace::ModelAnalysisReport report =
        laplace::analyze_hessian_structure(H, detector_options);

    return make_opt_pattern_info_from_report(report);
  } catch (const std::exception &e) {
    info.available = false;
    info.detected_structure = "unavailable";
    info.backend = "unavailable";
    info.solver = "unavailable";
    info.complexity = e.what();
    return info;
  } catch (...) {
    info.available = false;
    info.detected_structure = "unavailable";
    info.backend = "unavailable";
    info.solver = "unavailable";
    info.complexity = "unknown pattern-analysis failure";
    return info;
  }
}

template <typename Model>
LaplaceResult<Model> laplace_eval_at_u_star_persistent_structured(
    Model &model, ParameterVector &params, const std::vector<int> &fixed_idx,
    const std::vector<int> &random_idx, const Eigen::VectorXd &x,
    const std::vector<double> &u_star, had::ADGraph &graph,
    laplace::PersistentStructuredRuntimeState &structured_runtime,
    const LaplaceOptions &options = default_laplace_options()) {
  ADScope scope(graph);

  using Result = LaplaceResult<Model>;
  Result res;

  std::vector<AD> p_full;
  p_full.reserve(params.size());

  for (int i = 0; i < params.size(); ++i) {
    p_full.emplace_back(AD(0.0));
  }

  inject_fixed_params(x, p_full, fixed_idx);
  inject_random_params(u_star, p_full, random_idx);

  AD nll = model(p_full);

  scope.backward(nll);

  res.grad_x.resize(fixed_idx.size());
  for (size_t k = 0; k < fixed_idx.size(); ++k) {
    res.grad_x[k] = scope.grad(p_full[fixed_idx[k]]);
  }

  // Keep the existing exact logdet-gradient path for fixed effects.
  {
    Eigen::Map<const Eigen::VectorXd> u_star_eigen(
        u_star.data(), static_cast<Eigen::Index>(u_star.size()));

    Eigen::VectorXd g_logdet =
        laplace_logdet_gradient_exact(model, params, x, u_star_eigen, options);

    // laplace_logdet_gradient_exact builds temporary AD graphs.
    // Restore the graph for this outer evaluation before any further
    // grad/hess access through scope.
    had::g_ADGraph = &scope.graph;

    for (size_t k = 0; k < fixed_idx.size(); ++k) {
      res.grad_x[k] += g_logdet[static_cast<Eigen::Index>(k)];
    }
  }

  res.grad_u.resize(random_idx.size());
  for (size_t k = 0; k < random_idx.size(); ++k) {
    res.grad_u[k] = scope.grad(p_full[random_idx[k]]);
  }

  const auto &pattern = get_pattern(scope, p_full, random_idx);

  Eigen::SparseMatrix<double> H = extract_sparse_hessian(
      scope, p_full, random_idx, pattern, options.hessian_drop_tol);

  // Persistent structured bridge:
  //   First call: detect structure and choose backend.
  //   Later calls: update structured values only and reuse recommendation.
  laplace::StructureDetectorOptions detector_options;
  detector_options.prefer_dense_for_small_matrices = false;
  detector_options.dense_size_cutoff = 0;

  if (!structured_runtime.initialized) {
    structured_runtime.update_from_hessian(H, detector_options);
  } else {
    structured_runtime.update_values_only(H);
  }

  const double logdet = structured_runtime.logdet();

  res.value = value_of(nll) + 0.5 * logdet;

  return res;
}

struct LBFGSConvergedByGradient : public std::runtime_error {
  LBFGSConvergedByGradient()
      : std::runtime_error(
            "Quadra LBFGS reached requested gradient tolerance") {}
};

template <typename Model> class LBFGSObjective {
  void print(int iter, double fx, double gnorm) {
    const bool converged = std::isfinite(gnorm) && gnorm <= epsilon;
    std::cout << "L-BFGS: " << "outer eval = " << std::setw(3) << iter
              << ", fx = " << std::setw(14) << std::fixed
              << std::setprecision(6) << fx << ", |grad| = ";

    if (converged) {
      std::cout << "\033[1;32m";
    } else {
      std::cout << "\033[1;31m";
    }

    std::cout << std::setw(12) << std::fixed << std::setprecision(6) << gnorm
              << "\033[0m" << std::endl;
  }

public:
  double epsilon = 1e-6;
  Model &model;
  ParameterVector &params;
  std::vector<int> fixed_idx;
  std::vector<int> random_idx;
  LaplaceOptions options;

  int iter = 0;
  int print_every = 10;

  double last_fx = std::numeric_limits<double>::quiet_NaN();
  Eigen::VectorXd last_grad;
  Eigen::VectorXd last_x;
  std::vector<double> last_u_star;
  Eigen::VectorXd best_converged_x;
  Eigen::VectorXd best_converged_grad;
  std::vector<double> best_converged_u_star;
  double best_converged_fx = std::numeric_limits<double>::quiet_NaN();
  int best_converged_iter = 0;
  bool has_best_converged = false;

  // Best finite point seen by the optimizer, independent of the final
  // line-search bookkeeping state.
  double best_fx = std::numeric_limits<double>::infinity();
  double best_grad_norm = std::numeric_limits<double>::infinity();
  Eigen::VectorXd best_x;
  Eigen::VectorXd best_grad;
  std::vector<double> best_u_star;
  bool best_available = false;

  // Best point satisfying the configured fixed-effect gradient tolerance.
  double best_converged_grad_norm = std::numeric_limits<double>::infinity();
  laplace::PersistentStructuredRuntimeState structured_runtime;

  LBFGSObjective(Model &m, ParameterVector &p, std::vector<int> fixed,
                 std::vector<int> random,
                 const LaplaceOptions &opts = default_laplace_options())
      : model(m), params(p), fixed_idx(std::move(fixed)),
        random_idx(std::move(random)), options(opts) {
    laplace_pattern_cache().clear();
  }

  double operator()(const VectorXd &x, VectorXd &grad) {
    TapeContext tape;
    had::ADGraph &graph = tape.graph;

    ++iter;

    std::vector<double> u_star;
    const bool verbose_inner = ((iter % print_every) == 0) || iter == 1;

    try {
      u_star = solve_random_effects_laplace(model, params, x, fixed_idx,
                                            random_idx, graph);
      last_u_star = u_star;
    } catch (const std::exception &e) {
      std::cerr << "L-BFGS: random-effect mode solve failed; returning "
                   "penalty. reason="
                << e.what() << std::endl;
      const double penalty_gradient_scale = 1.0e3;

      for (int i = 0; i < grad.size(); ++i) {
        const double xi = (i < x.size() && std::isfinite(x[i])) ? x[i] : 1.0;
        grad[i] = penalty_gradient_scale * ((xi == 0.0) ? 1.0 : xi);
      }

      return std::numeric_limits<double>::max() / 100.0;
    }

    using Result = LaplaceResult<Model>;
    Result res;

    try {
      res = laplace_eval_at_u_star_persistent_structured(
          model, params, fixed_idx, random_idx, x, u_star, graph,
          structured_runtime, options);
    } catch (const std::exception &e) {
      std::cerr
          << "L-BFGS: Laplace evaluation failed; returning penalty. reason="
          << e.what() << std::endl;

      grad.resize(x.size());
      grad.setZero();
      last_grad = grad;
      last_x = x;
      last_fx = std::numeric_limits<double>::max() / 1.0e100;
      return last_fx;
    } catch (...) {
      std::cerr << "L-BFGS: Laplace evaluation failed with unknown exception; "
                   "returning penalty."
                << std::endl;

      grad.resize(x.size());
      grad.setZero();
      last_grad = grad;
      last_x = x;
      last_fx = std::numeric_limits<double>::max() / 1.0e100;
      return last_fx;
    }

    grad = to_eigen(res.grad_x);

    last_fx = res.value;
    last_grad = grad;
    last_x = x;

    const double gnorm = safe_eigen_norm(grad);

    if (std::isfinite(gnorm) && gnorm <= epsilon) {
      if (!has_best_converged || !std::isfinite(best_converged_fx) ||
          res.value < best_converged_fx) {
        best_converged_x = x;
        best_converged_grad = grad;
        best_converged_u_star = u_star;
        best_converged_fx = res.value;
        best_converged_iter = iter;
        has_best_converged = true;
      }

      print(iter, res.value, gnorm);
      throw LBFGSConvergedByGradient();
    }

    if ((iter % print_every) == 0 || iter == 1) {
      print(iter, res.value, gnorm);
    }

    return res.value;
  }
};

template <typename Model>
OptResult
optimize_lbfgs(Model &model, ParameterVector &params,
               const LaplaceOptions &options = default_laplace_options()) {
  using namespace LBFGSpp;
  using namespace Eigen;

  const auto fixed_idx = build_fixed_index(params);
  const auto random_idx = build_random_index(params);

  if (fixed_idx.empty()) {
    throw std::runtime_error(
        "No fixed parameters found — optimizer has zero dimension");
  }

  VectorXd x(static_cast<Eigen::Index>(fixed_idx.size()));
  for (size_t k = 0; k < fixed_idx.size(); ++k) {
    x[static_cast<Eigen::Index>(k)] =
        params.params[static_cast<size_t>(fixed_idx[k])].value;
  }

  LBFGSObjective<Model> fun(model, params, fixed_idx, random_idx, options);
  fun.print_every = 10;

  LBFGSParam<double> param;
  param.max_iterations = 400;
  // param.max_linesearch = 20;
  param.epsilon = 1.0e-4;
  fun.epsilon = param.epsilon;

  LBFGSSolver<double> solver(param);
  double fx = std::numeric_limits<double>::quiet_NaN();
  int niter = 0;

  try {
    niter = solver.minimize(fun, x, fx);

    // quadra_lbfgs_honest_convergence_report_v1
    double quadra_final_fixed_grad_norm =
        std::numeric_limits<double>::quiet_NaN();
    if (fun.last_grad.size() > 0) {
      quadra_final_fixed_grad_norm = 0.0;
      for (int quadra_i = 0; quadra_i < fun.last_grad.size(); ++quadra_i) {
        quadra_final_fixed_grad_norm +=
            fun.last_grad[quadra_i] * fun.last_grad[quadra_i];
      }
      quadra_final_fixed_grad_norm = std::sqrt(quadra_final_fixed_grad_norm);
    }

    const bool quadra_requested_tol_met =
        std::isfinite(quadra_final_fixed_grad_norm) &&
        quadra_final_fixed_grad_norm <= 1.0e-4;

    std::cout << "L-BFGS minimize status report" << std::endl;
    std::cout << "  iterations returned by solver: " << niter << std::endl;
    std::cout << "  final objective returned by solver: " << fx << std::endl;
    std::cout << "  final fixed-gradient norm: " << quadra_final_fixed_grad_norm
              << std::endl;
    std::cout << "  requested gradient tolerance: " << std::scientific << 1.0e-4
              << std::defaultfloat << std::endl;
    std::cout << "  configured max-iteration field: " << 400
              << " (LBFGSpp max_iterations)" << std::endl;
    std::cout << "  requested tolerance met: "
              << (quadra_requested_tol_met ? "yes" : "no") << std::endl;
    std::cout << "  outer convergence interpretation: "
              << (quadra_requested_tol_met
                      ? "converged to requested gradient tolerance"
                      : "stopped before requested gradient tolerance; inspect "
                        "LBFGS status/max iterations/line search")
              << std::endl;
  } catch (const LBFGSConvergedByGradient &) {
    if (fun.has_best_converged) {

      std::cout << "L-BFGS: stopped at first iterate satisfying requested "
                   "fixed-effect gradient tolerance."
                << std::endl;
    } else {
      throw;
    }
  } catch (const std::runtime_error &e) {
    const double gnorm = safe_eigen_norm(fun.last_grad);
    const double max_grad = (fun.last_grad.size() > 0)
                                ? fun.last_grad.cwiseAbs().maxCoeff()
                                : std::numeric_limits<double>::infinity();

    const std::string msg = e.what();

    const bool line_search_failed =
        msg.find("line search") != std::string::npos ||
        msg.find("Line search") != std::string::npos;

    // LBFGSpp may throw a line-search failure after the objective has
    // effectively plateaued. For public examples and diagnostic workflows,
    // return the best finite iterate instead of aborting, while keeping
    // result.converged honest via the stricter param.epsilon check below.
    const double convergence_like_grad = 2e-2;

    if (gnorm <= param.epsilon) {
      std::cout << "L-BFGS: optimization reached convergence criterion "
                << "(|grad| <= epsilon). max|grad| = " << max_grad << std::endl;

      if (fun.last_x.size() == x.size()) {
        x = fun.last_x;
      }

      fx = fun.last_fx;
      niter = fun.iter;
    } else if (line_search_failed && max_grad < convergence_like_grad) {
      std::cout
          << "L-BFGS: line search failed after a small fixed-effect gradient. "
          << "Returning the last finite iterate as a non-converged result. "
          << "max|grad| = " << max_grad << std::endl;

      if (fun.last_x.size() == x.size()) {
        x = fun.last_x;
      }

      fx = fun.last_fx;
      niter = fun.iter;
    } else {
      throw;
    }
  }

  OptResult result;

  Eigen::VectorXd selected_x;
  std::vector<double> selected_u_hat;
  double selected_fx = std::numeric_limits<double>::quiet_NaN();
  double selected_grad_norm = std::numeric_limits<double>::infinity();

  if (fun.has_best_converged) {
    selected_grad_norm = fun.best_converged_grad_norm;
  } else if (fun.best_available) {
    selected_x = fun.best_x;
    selected_u_hat = fun.best_u_star;
    selected_fx = fun.best_fx;
    selected_grad_norm = fun.best_grad_norm;
  } else if (fun.last_x.size() == x.size()) {
    selected_x = fun.last_x;
    selected_u_hat = fun.last_u_star;
    selected_fx = std::isfinite(fun.last_fx) ? fun.last_fx : fx;
    selected_grad_norm = safe_eigen_norm(fun.last_grad);
  } else {
    selected_x = x;
    selected_u_hat = fun.last_u_star;
    selected_fx = fx;
    selected_grad_norm = safe_eigen_norm(fun.last_grad);
  }

  for (size_t k = 0; k < fixed_idx.size(); ++k) {
    params.params[static_cast<size_t>(fixed_idx[k])].value =
        selected_x[static_cast<Eigen::Index>(k)];
  }

  result.par.assign(selected_x.data(), selected_x.data() + selected_x.size());
  result.u_hat = selected_u_hat;
  result.fixed_index = fixed_idx;
  result.random_index = random_idx;

  result.value = selected_fx;
  result.iterations = niter;
  result.grad_norm = std::isfinite(selected_grad_norm)
                         ? selected_grad_norm
                         : std::numeric_limits<double>::infinity();

  result.converged =
      std::isfinite(result.grad_norm) && result.grad_norm <= param.epsilon;

  result.message =
      result.converged
          ? "converged to requested fixed-effect gradient tolerance"
          : "stopped before requested fixed-effect gradient tolerance";

  Eigen::VectorXd pattern_x = selected_x;

  result.pattern = analyze_final_random_effect_pattern(
      model, params, pattern_x, result.u_hat, fixed_idx, random_idx, options);

  return result;
}

} // namespace quadra

#endif
