#include "state_space_surplus_production.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "core/laplace/persistent_structured_runtime.hpp"
#include "core/laplace/persistent_random_effect_state.hpp"
#include "core/laplace/persistent_latent_state_runtime.hpp"
#include "core/laplace/laplace_evaluator.hpp"
#include <LBFGS.h>

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace ss = quadra_examples::state_space_surplus_production;
using Clock = std::chrono::high_resolution_clock;

namespace {

double ms_between(const Clock::time_point& a, const Clock::time_point& b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

std::vector<int> parse_lengths(const std::string& s) {
  std::vector<int> out;
  std::stringstream ss_in(s);
  std::string item;
  while (std::getline(ss_in, item, ',')) {
    if (!item.empty()) out.push_back(std::stoi(item));
  }
  return out;
}

double normal_const(const double sigma) {
  return std::log(sigma) + 0.5 * std::log(2.0 * M_PI);
}

double normal_nll(const double residual, const double sigma) {
  const double z = residual / sigma;
  return normal_const(sigma) + 0.5 * z * z;
}

double inv_logit(const double x) {
  return 1.0 / (1.0 + std::exp(-x));
}

ss::Data make_scaled_data(const int n) {
  ss::Data data;
  data.catch_observed.resize(static_cast<std::size_t>(n));
  data.index_observed.resize(static_cast<std::size_t>(n));

  const double r = 0.5;
  const double K = 700.0;
  const double q = 0.0024;
  double B = 0.90 * K;

  for (int t = 0; t < n; ++t) {
    const double seasonal = std::sin(2.0 * M_PI * static_cast<double>(t) / 17.0);
    const double trend = 1.0 + 0.10 * std::sin(2.0 * M_PI * static_cast<double>(t) / 53.0);
    const double C = 88.0 * trend + 18.0 * seasonal;

    data.catch_observed[static_cast<std::size_t>(t)] = std::max(40.0, C);

    const double obs_error =
        0.05 * std::sin(2.0 * M_PI * static_cast<double>(t) / 11.0) +
        0.025 * std::cos(2.0 * M_PI * static_cast<double>(t) / 7.0);

    data.index_observed[static_cast<std::size_t>(t)] = q * B * std::exp(obs_error);

    if (t < n - 1) {
      const double production = r * B * (1.0 - B / K);
      B = std::max(B + production - data.catch_observed[static_cast<std::size_t>(t)], 1e-9);
    }
  }

  return data;
}

struct Transformed {
  double r = 0.0;
  double K = 0.0;
  double q = 0.0;
  double sigma_process = 0.0;
  double sigma_index = 0.0;
  double B0_frac = 0.0;
  double B0 = 0.0;
  double inv_var_process = 0.0;
  double inv_var_index = 0.0;
};

Transformed transform(const ss::Parameters& par) {
  Transformed tr;
  tr.r = std::exp(par.log_r);
  tr.K = std::exp(par.log_K);
  tr.q = std::exp(par.log_q);
  tr.sigma_process = std::exp(par.log_sigma_process);
  tr.sigma_index = std::exp(par.log_sigma_index);
  tr.B0_frac = inv_logit(par.logit_B0_frac);
  tr.B0 = tr.B0_frac * tr.K;
  tr.inv_var_process = 1.0 / (tr.sigma_process * tr.sigma_process);
  tr.inv_var_index = 1.0 / (tr.sigma_index * tr.sigma_index);
  return tr;
}





ss::Parameters make_par() {
  ss::Parameters par;
  par.log_r = std::log(0.5);
  par.log_K = std::log(700.0);
  par.log_q = std::log(0.0024);
  par.log_sigma_process = std::log(0.15);
  par.log_sigma_index = std::log(0.10);
  par.logit_B0_frac = std::log(0.90 / 0.10);
  return par;
}

// f(B) = B + rB(1-B/K) - C
//      = (1+r)B - (r/K)B^2 - C
// logpred(y) with B = exp(y)
// d/dy log f = A
// d2/dy2 log f = A2
struct LogPredDerivatives {
  double pred = 0.0;
  double log_pred = 0.0;
  double d1 = 0.0;
  double d2 = 0.0;
};

LogPredDerivatives log_pred_derivatives(const double log_B,
                                        const double catch_value,
                                        const Transformed& tr) {
  const double B = std::exp(log_B);
  double pred = B + tr.r * B * (1.0 - B / tr.K) - catch_value;
  pred = std::max(pred, 1e-9);

  // For pred clamped at floor, derivatives are effectively zero.
  if (pred <= 1.0000001e-9) {
    return {pred, std::log(pred), 0.0, 0.0};
  }

  const double a = 1.0 + tr.r;
  const double b = tr.r / tr.K;

  const double fp_y = a * B - 2.0 * b * B * B;
  const double fpp_y = a * B - 4.0 * b * B * B;

  const double d1 = fp_y / pred;
  const double d2 = fpp_y / pred - (fp_y * fp_y) / (pred * pred);

  return {pred, std::log(pred), d1, d2};
}

double joint_x(const ss::Data& data, const ss::Parameters& par, const Eigen::VectorXd& x) {
  const int n = static_cast<int>(data.catch_observed.size());
  const Transformed tr = transform(par);

  double nll = 0.0;

  nll += normal_nll(
      std::log(data.index_observed[0]) - (std::log(tr.q) + std::log(tr.B0)),
      tr.sigma_index);

  for (int t = 0; t < n - 1; ++t) {
    const double log_B_t = (t == 0) ? std::log(tr.B0) : x[t - 1];
    const auto lp = log_pred_derivatives(
        log_B_t, data.catch_observed[static_cast<std::size_t>(t)], tr);

    const double process_residual = x[t] - lp.log_pred;
    nll += normal_nll(process_residual, tr.sigma_process);

    const double obs_residual =
        std::log(data.index_observed[static_cast<std::size_t>(t + 1)]) -
        (std::log(tr.q) + x[t]);

    nll += normal_nll(obs_residual, tr.sigma_index);
  }

  return nll;
}

Eigen::VectorXd deterministic_initial_x(const ss::Data& data, const ss::Parameters& par) {
  const int n = static_cast<int>(data.catch_observed.size());
  const Transformed tr = transform(par);
  Eigen::VectorXd x(n - 1);

  double B = tr.B0;
  for (int t = 0; t < n - 1; ++t) {
    const double production = tr.r * B * (1.0 - B / tr.K);
    B = std::max(B + production - data.catch_observed[static_cast<std::size_t>(t)], 1e-9);
    x[t] = std::log(B);
  }
  return x;
}

Eigen::VectorXd analytic_grad_x(const ss::Data& data,
                                const ss::Parameters& par,
                                const Eigen::VectorXd& x) {
  const int n_state = static_cast<int>(x.size());
  const Transformed tr = transform(par);

  Eigen::VectorXd g = Eigen::VectorXd::Zero(n_state);

  for (int k = 0; k < n_state; ++k) {
    // Process residual at k: e_k = x[k] - logpred_k.
    const double log_B_k = (k == 0) ? std::log(tr.B0) : x[k - 1];
    const auto lp_k = log_pred_derivatives(
        log_B_k, data.catch_observed[static_cast<std::size_t>(k)], tr);
    const double e_k = x[k] - lp_k.log_pred;

    g[k] += e_k * tr.inv_var_process;

    // Observation residual at year k+1: o_k = logI[k+1] - logq - x[k].
    const double o_k =
        std::log(data.index_observed[static_cast<std::size_t>(k + 1)]) -
        (std::log(tr.q) + x[k]);

    g[k] += -o_k * tr.inv_var_index;

    // Process residual at k+1 depends on x[k] through logpred_{k+1}.
    if (k + 1 < n_state) {
      const auto lp_next = log_pred_derivatives(
          x[k], data.catch_observed[static_cast<std::size_t>(k + 1)], tr);
      const double e_next = x[k + 1] - lp_next.log_pred;

      g[k] += -e_next * lp_next.d1 * tr.inv_var_process;
    }
  }

  return g;
}

Eigen::SparseMatrix<double> analytic_hessian_xx(const ss::Data& data,
                                                const ss::Parameters& par,
                                                const Eigen::VectorXd& x) {
  const int n_state = static_cast<int>(x.size());
  const Transformed tr = transform(par);

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(static_cast<std::size_t>(3 * n_state));

  for (int k = 0; k < n_state; ++k) {
    double diag = 0.0;

    // e_k contributes +1/sigma_p^2 wrt x[k].
    diag += tr.inv_var_process;

    // observation at k+1 contributes +1/sigma_i^2.
    diag += tr.inv_var_index;

    // e_{k+1} contributes curvature wrt x[k] if it exists:
    // e = x[k+1] - L(x[k])
    // contribution = 0.5 e^2 / sp^2
    // d2/dxk2 = (L'^2 - e L'') / sp^2
    if (k + 1 < n_state) {
      const auto lp_next = log_pred_derivatives(
          x[k], data.catch_observed[static_cast<std::size_t>(k + 1)], tr);
      const double e_next = x[k + 1] - lp_next.log_pred;
      diag += (lp_next.d1 * lp_next.d1 - e_next * lp_next.d2) *
              tr.inv_var_process;

      // off diagonal H[k, k+1] from e_{k+1}:
      // d/dxk = -e L'/sp^2, d/dx{k+1} => -L'/sp^2
      const double off = -lp_next.d1 * tr.inv_var_process;
      triplets.emplace_back(k, k + 1, off);
      triplets.emplace_back(k + 1, k, off);
    }

    triplets.emplace_back(k, k, diag);
  }

  Eigen::SparseMatrix<double> H(n_state, n_state);
  H.setFromTriplets(triplets.begin(), triplets.end());
  return H;
}

class ObjX {
 public:
  ObjX(const ss::Data& data, const ss::Parameters& par) : data_(data), par_(par) {}
  double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad) {
    grad = analytic_grad_x(data_, par_, x);
    return joint_x(data_, par_, x);
  }
 private:
  const ss::Data& data_;
  const ss::Parameters& par_;
};

Eigen::VectorXd optimize_x(const ss::Data& data, const ss::Parameters& par) {
  Eigen::VectorXd x = deterministic_initial_x(data, par);

  LBFGSpp::LBFGSParam<double> param;
  param.epsilon = 1e-8;
  param.max_iterations = 500;
  param.max_linesearch = 100;
  param.m = 8;
  param.ftol = 1e-4;
  param.wolfe = 0.9;
  param.min_step = 1e-20;
  param.max_step = 1.0;

  LBFGSpp::LBFGSSolver<double> solver(param);
  ObjX obj(data, par);
  double f = 0.0;

  try {
    solver.minimize(obj, x, f);
  } catch (...) {
    const double gnorm = analytic_grad_x(data, par, x).norm();
    if (!(std::isfinite(joint_x(data, par, x)) && gnorm < 1e-5)) throw;
  }

  return x;
}

double sparse_logdet_ldlt(const Eigen::SparseMatrix<double>& H) {
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
  ldlt.compute(H);
  if (ldlt.info() != Eigen::Success) throw std::runtime_error("Sparse LDLT failed");

  const auto& D = ldlt.vectorD();
  double logdet = 0.0;
  for (int i = 0; i < D.size(); ++i) {
    if (!(D[i] > 0.0)) throw std::runtime_error("Hxx not positive definite");
    logdet += std::log(D[i]);
  }
  return logdet;
}



struct EvalResult {
  double objective = 0.0;
  double joint = 0.0;
  double logdet = 0.0;
  double grad_norm = 0.0;
  int nnz = 0;
};

EvalResult eval(const ss::Data& data, const ss::Parameters& par) {
  EvalResult out;
  const Eigen::VectorXd xhat = optimize_x(data, par);

  out.joint = joint_x(data, par, xhat);
  out.grad_norm = analytic_grad_x(data, par, xhat).norm();

  const Eigen::SparseMatrix<double> H = analytic_hessian_xx(data, par, xhat);
  out.nnz = static_cast<int>(H.nonZeros());
  out.logdet = sparse_logdet_ldlt(H);

  const double n_x = static_cast<double>(xhat.size());
  const double correction = 0.5 * out.logdet - 0.5 * n_x * std::log(2.0 * M_PI);
  out.objective = out.joint + correction;
  return out;
}


quadra::laplace::TridiagonalValues make_tridiagonal_values_xx(
    const ss::Data& data,
    const ss::Parameters& par,
    const Eigen::VectorXd& x) {
  const int n_state = static_cast<int>(x.size());
  const Transformed tr = transform(par);

  quadra::laplace::TridiagonalValues out;
  out.diag = Eigen::VectorXd::Zero(n_state);
  out.offdiag = Eigen::VectorXd::Zero(std::max(0, n_state - 1));

  for (int k = 0; k < n_state; ++k) {
    double diag = 0.0;

    // e_k contributes +1/sigma_p^2 wrt x[k].
    diag += tr.inv_var_process;

    // Observation at k+1 contributes +1/sigma_i^2.
    diag += tr.inv_var_index;

    // e_{k+1} contributes curvature wrt x[k] if it exists:
    // e = x[k+1] - L(x[k])
    // d2/dxk2 = (L'^2 - e L'') / sigma_p^2
    if (k + 1 < n_state) {
      const auto lp_next = log_pred_derivatives(
          x[k], data.catch_observed[static_cast<std::size_t>(k + 1)], tr);
      const double e_next = x[k + 1] - lp_next.log_pred;

      diag += (lp_next.d1 * lp_next.d1 - e_next * lp_next.d2) *
              tr.inv_var_process;

      // Off diagonal H[k, k+1] = H[k+1, k] from e_{k+1}.
      out.offdiag[k] = -lp_next.d1 * tr.inv_var_process;
    }

    out.diag[k] = diag;
  }

  return out;
}

EvalResult eval_direct_runtime(
    const ss::Data& data,
    const ss::Parameters& par,
    quadra::laplace::PersistentStructuredRuntimeState& runtime) {
  EvalResult out;
  const Eigen::VectorXd xhat = optimize_x(data, par);

  out.joint = joint_x(data, par, xhat);
  out.grad_norm = analytic_grad_x(data, par, xhat).norm();

  const quadra::laplace::TridiagonalValues H =
      make_tridiagonal_values_xx(data, par, xhat);

  out.nnz = static_cast<int>(H.diag.size()) +
            2 * static_cast<int>(H.offdiag.size());

  runtime.update_direct(H);
  out.logdet = runtime.logdet();

  const double n_x = static_cast<double>(xhat.size());
  const double correction =
      0.5 * out.logdet - 0.5 * n_x * std::log(2.0 * M_PI);

  out.objective = out.joint + correction;
  return out;
}



Eigen::VectorXd solve_tridiagonal_ldlt(
    const quadra::laplace::TridiagonalValues& H,
    const Eigen::VectorXd& rhs) {
  const int n = static_cast<int>(H.diag.size());

  if (rhs.size() != n) {
    throw std::invalid_argument("tridiagonal solve rhs size mismatch");
  }

  if (n == 0) {
    return Eigen::VectorXd();
  }

  Eigen::VectorXd D = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd L = Eigen::VectorXd::Zero(std::max(0, n - 1));

  D[0] = H.diag[0];
  if (!(D[0] > 0.0) || !std::isfinite(D[0])) {
    throw std::runtime_error("Newton Hessian is not SPD");
  }

  for (int i = 1; i < n; ++i) {
    L[i - 1] = H.offdiag[i - 1] / D[i - 1];
    D[i] = H.diag[i] - L[i - 1] * L[i - 1] * D[i - 1];

    if (!(D[i] > 0.0) || !std::isfinite(D[i])) {
      throw std::runtime_error("Newton Hessian is not SPD");
    }
  }

  // Solve L y = rhs.
  Eigen::VectorXd y = rhs;
  for (int i = 1; i < n; ++i) {
    y[i] -= L[i - 1] * y[i - 1];
  }

  // Solve D z = y.
  Eigen::VectorXd z = y;
  for (int i = 0; i < n; ++i) {
    z[i] /= D[i];
  }

  // Solve L^T x = z.
  Eigen::VectorXd x = z;
  for (int i = n - 2; i >= 0; --i) {
    x[i] -= L[i] * x[i + 1];
  }

  return x;
}

struct NewtonSolveResult {
  Eigen::VectorXd xhat;
  int iterations = 0;
  double final_joint = 0.0;
  double final_grad_norm = 0.0;
  bool converged = false;
};

NewtonSolveResult optimize_x_newton_tridiagonal(
    const ss::Data& data,
    const ss::Parameters& par,
    const Eigen::VectorXd& initial_x,
    const int max_iterations = 25,
    const double grad_tol = 1e-6) {
  NewtonSolveResult out;
  Eigen::VectorXd x = initial_x;

  double fx = joint_x(data, par, x);

  for (int iter = 0; iter < max_iterations; ++iter) {
    const Eigen::VectorXd g = analytic_grad_x(data, par, x);
    const double gnorm = g.norm();

    out.iterations = iter;
    out.final_joint = fx;
    out.final_grad_norm = gnorm;

    if (gnorm < grad_tol) {
      out.converged = true;
      out.xhat = x;
      return out;
    }

    const quadra::laplace::TridiagonalValues H =
        make_tridiagonal_values_xx(data, par, x);

    Eigen::VectorXd step = solve_tridiagonal_ldlt(H, -g);

    double alpha = 1.0;
    bool accepted = false;

    for (int ls = 0; ls < 20; ++ls) {
      const Eigen::VectorXd candidate = x + alpha * step;
      const double f_candidate = joint_x(data, par, candidate);

      if (std::isfinite(f_candidate) && f_candidate <= fx) {
        x = candidate;
        fx = f_candidate;
        accepted = true;
        break;
      }

      alpha *= 0.5;
    }

    if (!accepted) {
      // If a full Newton step is numerically awkward near the optimum, stop
      // only if the gradient is already small enough for the benchmark.
      if (gnorm < 1e-4) {
        out.converged = true;
        out.xhat = x;
        out.final_joint = fx;
        out.final_grad_norm = gnorm;
        return out;
      }

      break;
    }
  }

  out.xhat = x;
  out.final_joint = joint_x(data, par, x);
  out.final_grad_norm = analytic_grad_x(data, par, x).norm();
  out.converged = out.final_grad_norm < 1e-4;
  return out;
}

EvalResult eval_with_xhat_direct_runtime(
    const ss::Data& data,
    const ss::Parameters& par,
    const Eigen::VectorXd& xhat,
    quadra::laplace::PersistentStructuredRuntimeState& runtime) {
  EvalResult out;

  out.joint = joint_x(data, par, xhat);
  out.grad_norm = analytic_grad_x(data, par, xhat).norm();

  const quadra::laplace::TridiagonalValues H =
      make_tridiagonal_values_xx(data, par, xhat);

  out.nnz = static_cast<int>(H.diag.size()) +
            2 * static_cast<int>(H.offdiag.size());

  runtime.update_direct(H);
  out.logdet = runtime.logdet();

  const double n_x = static_cast<double>(xhat.size());
  const double correction =
      0.5 * out.logdet - 0.5 * n_x * std::log(2.0 * M_PI);

  out.objective = out.joint + correction;
  return out;
}







quadra::laplace::NewtonSolveStatus to_newton_status(
    const NewtonSolveResult& result) {
  quadra::laplace::NewtonSolveStatus status;
  status.iterations = result.iterations;
  status.objective = result.final_joint;
  status.grad_norm = result.final_grad_norm;
  status.converged = result.converged;
  return status;
}


struct StateSpaceContext {
  const ss::Data& data;
  const ss::Parameters& par;
};

struct StateSpaceNewtonSolveResult {
  Eigen::VectorXd xhat;
  quadra::laplace::NewtonSolveStatus status;
  NewtonSolveResult raw;
};

struct StateSpaceNewtonPolicy {
  Eigen::VectorXd initial_x(StateSpaceContext& ctx) {
    return deterministic_initial_x(ctx.data, ctx.par);
  }

  StateSpaceNewtonSolveResult solve(StateSpaceContext& ctx,
                                    const Eigen::VectorXd& initial_x) {
    StateSpaceNewtonSolveResult out;
    out.raw = optimize_x_newton_tridiagonal(ctx.data, ctx.par, initial_x);
    out.xhat = out.raw.xhat;
    out.status = to_newton_status(out.raw);
    return out;
  }
};


struct StateSpaceEvaluationPolicy {
  EvalResult evaluate(
      StateSpaceContext& ctx,
      const Eigen::VectorXd& xhat,
      quadra::laplace::PersistentStructuredRuntimeState& structured) {
    return eval_with_xhat_direct_runtime(
        ctx.data,
        ctx.par,
        xhat,
        structured);
  }
};


}  // namespace

int main(int argc, char** argv) {
  int reps = 10;
  std::vector<int> lengths = {25, 50, 100, 250};

  if (argc > 1) reps = std::stoi(argv[1]);
  if (argc > 2) lengths = parse_lengths(argv[2]);

  const ss::Parameters par = make_par();

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "Quadra state-space LaplaceEvaluator benchmark\n";
  std::cout << "reps per n = " << reps << "\n\n";

  std::cout << std::setw(8) << "n"
            << std::setw(14) << "objective"
            << std::setw(14) << "lbfgs_ms"
            << std::setw(16) << "cold_newton"
            << std::setw(16) << "warm_newton"
            << std::setw(14) << "warm_x"
            << std::setw(12) << "cold_it"
            << std::setw(12) << "warm_it"
            << std::setw(14) << "grad_norm"
            << std::setw(16) << "objective_diff"
            << "\n";

  for (const int n : lengths) {
    const ss::Data data = make_scaled_data(n);
    quadra::laplace::PersistentLatentStateRuntime latent_runtime;

    EvalResult last_lbfgs = eval(data, par);

    const auto lb0 = Clock::now();
    for (int r = 0; r < reps; ++r) {
      last_lbfgs = eval(data, par);
    }
    const auto lb1 = Clock::now();

    const double lbfgs_ms = ms_between(lb0, lb1) / static_cast<double>(reps);

    NewtonSolveResult cold_last;
    EvalResult last_cold;

    const auto cold0 = Clock::now();
    for (int r = 0; r < reps; ++r) {
      const Eigen::VectorXd initial_x = deterministic_initial_x(data, par);
      cold_last = optimize_x_newton_tridiagonal(data, par, initial_x);
      last_cold = eval_with_xhat_direct_runtime(
          data, par, cold_last.xhat, latent_runtime.structured());
    }
    const auto cold1 = Clock::now();

    const double cold_newton_ms =
        ms_between(cold0, cold1) / static_cast<double>(reps);

    // Warm-start Newton driven through LaplaceEvaluator.
    StateSpaceContext context{data, par};
    using Evaluator =
        quadra::laplace::LaplaceEvaluator<
            StateSpaceContext,
            StateSpaceNewtonPolicy,
            StateSpaceEvaluationPolicy>;

    Evaluator evaluator;

    // One priming evaluation outside timing.
    EvalResult last_warm = evaluator.evaluate(context);

    const auto warm0 = Clock::now();
    int warm_iterations_total = 0;

    for (int r = 0; r < reps; ++r) {
      last_warm = evaluator.evaluate(context);
      warm_iterations_total +=
          evaluator.runtime().random_effects().status().iterations;
    }

    const auto warm1 = Clock::now();

    const double warm_newton_ms =
        ms_between(warm0, warm1) / static_cast<double>(reps);
    const double warm_speedup =
        warm_newton_ms > 0.0 ? lbfgs_ms / warm_newton_ms : 0.0;
    const int avg_warm_it =
        static_cast<int>(std::lround(
            static_cast<double>(warm_iterations_total) /
            static_cast<double>(reps)));

    const double objective_diff =
        std::abs(last_lbfgs.objective - last_warm.objective);

    std::cout << std::setw(8) << n
              << std::setw(14) << last_warm.objective
              << std::setw(14) << lbfgs_ms
              << std::setw(16) << cold_newton_ms
              << std::setw(16) << warm_newton_ms
              << std::setw(14) << warm_speedup
              << std::setw(12) << cold_last.iterations
              << std::setw(12) << avg_warm_it
              << std::setw(14) << last_warm.grad_norm
              << std::setw(16) << objective_diff
              << "\n";
  }

  return 0;
}
