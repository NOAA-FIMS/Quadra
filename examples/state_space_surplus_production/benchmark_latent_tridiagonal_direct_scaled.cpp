#include "state_space_surplus_production.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <LBFGS.h>

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
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

double inv_logit(const double x) {
  return 1.0 / (1.0 + std::exp(-x));
}

double normal_nll(const double residual, const double sigma) {
  const double z = residual / sigma;
  return std::log(sigma) + 0.5 * std::log(2.0 * M_PI) + 0.5 * z * z;
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
  double r;
  double K;
  double q;
  double sigma_process;
  double sigma_index;
  double B0_frac;
  double B0;
};

Transformed transform(const ss::Parameters& par) {
  Transformed out;
  out.r = std::exp(par.log_r);
  out.K = std::exp(par.log_K);
  out.q = std::exp(par.log_q);
  out.sigma_process = std::exp(par.log_sigma_process);
  out.sigma_index = std::exp(par.log_sigma_index);
  out.B0_frac = inv_logit(par.logit_B0_frac);
  out.B0 = out.B0_frac * out.K;
  return out;
}

double predicted_next_biomass(const double B, const double catch_value, const Transformed& tr) {
  const double production = tr.r * B * (1.0 - B / tr.K);
  return std::max(B + production - catch_value, 1e-9);
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
    const double B_t = std::exp(log_B_t);
    const double pred_next =
        predicted_next_biomass(B_t, data.catch_observed[static_cast<std::size_t>(t)], tr);

    nll += normal_nll(x[t] - std::log(pred_next), tr.sigma_process);

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
    B = predicted_next_biomass(B, data.catch_observed[static_cast<std::size_t>(t)], tr);
    x[t] = std::log(B);
  }
  return x;
}

Eigen::VectorXd fd_grad_x(const ss::Data& data, const ss::Parameters& par, const Eigen::VectorXd& x) {
  Eigen::VectorXd grad(x.size());
  for (int i = 0; i < x.size(); ++i) {
    const double h = 1e-5 * (1.0 + std::abs(x[i]));
    Eigen::VectorXd xp = x;
    Eigen::VectorXd xm = x;
    xp[i] += h;
    xm[i] -= h;
    grad[i] = (joint_x(data, par, xp) - joint_x(data, par, xm)) / (2.0 * h);
  }
  return grad;
}

class ObjX {
 public:
  ObjX(const ss::Data& data, const ss::Parameters& par) : data_(data), par_(par) {}
  double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad) {
    const double f = joint_x(data_, par_, x);
    grad = fd_grad_x(data_, par_, x);
    return f;
  }
 private:
  const ss::Data& data_;
  const ss::Parameters& par_;
};

Eigen::VectorXd optimize_x(const ss::Data& data, const ss::Parameters& par) {
  Eigen::VectorXd x = deterministic_initial_x(data, par);

  LBFGSpp::LBFGSParam<double> param;
  param.epsilon = 1e-7;
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
    const double gnorm = fd_grad_x(data, par, x).norm();
    if (!(std::isfinite(joint_x(data, par, x)) && gnorm < 1e-3)) throw;
  }
  return x;
}

Eigen::SparseMatrix<double> fd_tridiagonal_hessian_xx(const ss::Data& data,
                                                      const ss::Parameters& par,
                                                      const Eigen::VectorXd& x) {
  const int n = static_cast<int>(x.size());
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(static_cast<std::size_t>(3 * n));

  for (int j = 0; j < n; ++j) {
    const double h = 1e-5 * (1.0 + std::abs(x[j]));
    Eigen::VectorXd xp = x;
    Eigen::VectorXd xm = x;
    xp[j] += h;
    xm[j] -= h;

    const Eigen::VectorXd gp = fd_grad_x(data, par, xp);
    const Eigen::VectorXd gm = fd_grad_x(data, par, xm);

    for (int i = std::max(0, j - 1); i <= std::min(n - 1, j + 1); ++i) {
      const double hij = (gp[i] - gm[i]) / (2.0 * h);
      if (std::abs(hij) > 1e-12) triplets.emplace_back(i, j, hij);
    }
  }

  Eigen::SparseMatrix<double> H(n, n);
  H.setFromTriplets(triplets.begin(), triplets.end());
  Eigen::SparseMatrix<double> Hsym = 0.5 * (H + Eigen::SparseMatrix<double>(H.transpose()));
  return Hsym;
}


struct TridiagonalValues {
  Eigen::VectorXd diag;
  Eigen::VectorXd offdiag;
};

TridiagonalValues fd_tridiagonal_values_xx(const ss::Data& data,
                                           const ss::Parameters& par,
                                           const Eigen::VectorXd& x) {
  const int n = static_cast<int>(x.size());
  TridiagonalValues out;
  out.diag = Eigen::VectorXd::Zero(n);
  out.offdiag = Eigen::VectorXd::Zero(std::max(0, n - 1));

  for (int j = 0; j < n; ++j) {
    const double h = 1e-5 * (1.0 + std::abs(x[j]));
    Eigen::VectorXd xp = x;
    Eigen::VectorXd xm = x;
    xp[j] += h;
    xm[j] -= h;

    const Eigen::VectorXd gp = fd_grad_x(data, par, xp);
    const Eigen::VectorXd gm = fd_grad_x(data, par, xm);

    out.diag[j] = (gp[j] - gm[j]) / (2.0 * h);
    if (j + 1 < n) {
      out.offdiag[j] = (gp[j + 1] - gm[j + 1]) / (2.0 * h);
    }
  }

  return out;
}

double logdet_tridiagonal_values_ldlt(const TridiagonalValues& H) {
  const int n = static_cast<int>(H.diag.size());
  if (n == 0) return 0.0;
  double d_prev = H.diag[0];
  if (!(d_prev > 0.0)) {
    throw std::runtime_error("Tridiagonal value Hessian is not positive definite");
  }
  double logdet = std::log(d_prev);
  for (int i = 1; i < n; ++i) {
    const double e = H.offdiag[i - 1];
    const double d = H.diag[i] - (e * e) / d_prev;
    if (!(d > 0.0)) {
      throw std::runtime_error("Tridiagonal value Hessian is not positive definite");
    }
    logdet += std::log(d);
    d_prev = d;
  }
  return logdet;
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
  out.grad_norm = fd_grad_x(data, par, xhat).norm();

  const TridiagonalValues H = fd_tridiagonal_values_xx(data, par, xhat);
  out.nnz = static_cast<int>(H.diag.size()) +
            2 * static_cast<int>(H.offdiag.size());
  out.logdet = logdet_tridiagonal_values_ldlt(H);

  const double n_x = static_cast<double>(xhat.size());
  const double correction = 0.5 * out.logdet - 0.5 * n_x * std::log(2.0 * M_PI);
  out.objective = out.joint + correction;
  return out;
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

}  // namespace

int main(int argc, char** argv) {
  int reps = 10;
  std::vector<int> lengths = {25, 50, 100, 250};

  if (argc > 1) reps = std::stoi(argv[1]);
  if (argc > 2) lengths = parse_lengths(argv[2]);

  const ss::Parameters par = make_par();

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "Quadra scaled direct latent-state tridiagonal Laplace benchmark\n";
  std::cout << "reps per n = " << reps << "\n\n";

  std::cout << std::setw(8) << "n"
            << std::setw(14) << "objective"
            << std::setw(14) << "joint"
            << std::setw(14) << "logdet"
            << std::setw(14) << "nnz"
            << std::setw(14) << "grad_norm"
            << std::setw(14) << "avg_ms"
            << "\n";

  for (const int n : lengths) {
    const ss::Data data = make_scaled_data(n);
    EvalResult last = eval(data, par);

    const auto t0 = Clock::now();
    for (int r = 0; r < reps; ++r) last = eval(data, par);
    const auto t1 = Clock::now();

    const double avg_ms = ms_between(t0, t1) / static_cast<double>(reps);

    std::cout << std::setw(8) << n
              << std::setw(14) << last.objective
              << std::setw(14) << last.joint
              << std::setw(14) << last.logdet
              << std::setw(14) << last.nnz
              << std::setw(14) << last.grad_norm
              << std::setw(14) << avg_ms
              << "\n";
  }

  return 0;
}
