#include <Eigen/Dense>

#include "LBFGS.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Data {
  std::vector<double> catch_mt;
  std::vector<double> index;
};

struct FixedParameters {
  double r = 0.30;
  double K = 900.0;
  double q = 0.00115;
  double sigma_process = 0.10;
  double sigma_index = 0.08;
  double B0_frac = 0.82;
};

struct LogPredDerivatives {
  double log_pred = 0.0;
  double d1 = 0.0;
  double d2 = 0.0;
};

struct NewtonResult {
  Eigen::VectorXd xhat;
  double joint = std::numeric_limits<double>::quiet_NaN();
  double grad_norm = std::numeric_limits<double>::quiet_NaN();
  int iterations = 0;
  bool converged = false;
};

struct FitResult {
  FixedParameters par;
  Eigen::VectorXd theta;
  Eigen::VectorXd xhat;
  double objective = std::numeric_limits<double>::quiet_NaN();
  double joint = std::numeric_limits<double>::quiet_NaN();
  double logdet = std::numeric_limits<double>::quiet_NaN();
  double grad_norm = std::numeric_limits<double>::quiet_NaN();
  int outer_iterations = 0;
  int inner_iterations = 0;
  int objective_evaluations = 0;
};

double inv_logit(const double x) {
  if (x >= 0.0) {
    const double z = std::exp(-x);
    return 1.0 / (1.0 + z);
  }
  const double z = std::exp(x);
  return z / (1.0 + z);
}

double bounded(const double z, const double lo, const double hi) {
  return lo + (hi - lo) * inv_logit(z);
}

double unbounded_from_value(const double x, const double lo, const double hi) {
  const double p = std::min(0.999999, std::max(0.000001, (x - lo) / (hi - lo)));
  return std::log(p / (1.0 - p));
}

FixedParameters theta_to_par(const Eigen::VectorXd &theta) {
  if (theta.size() != 3) {
    throw std::runtime_error("theta must have length 3");
  }

  FixedParameters par;
  par.r = bounded(theta[0], 0.05, 0.80);
  par.K = bounded(theta[1], 400.0, 2000.0);
  par.q = bounded(theta[2], 0.0002, 0.0050);

  // Fixed nuisance quantities for a stable public example.
  par.sigma_process = 0.10;
  par.sigma_index = 0.08;
  par.B0_frac = 0.82;

  return par;
}

Eigen::VectorXd par_to_theta(const FixedParameters &par) {
  Eigen::VectorXd theta(3);
  theta << unbounded_from_value(par.r, 0.05, 0.80),
      unbounded_from_value(par.K, 400.0, 2000.0),
      unbounded_from_value(par.q, 0.0002, 0.0050);
  return theta;
}

double B0(const FixedParameters &par) { return par.B0_frac * par.K; }

double inv_var_process(const FixedParameters &par) {
  return 1.0 / (par.sigma_process * par.sigma_process);
}

double inv_var_index(const FixedParameters &par) {
  return 1.0 / (par.sigma_index * par.sigma_index);
}

LogPredDerivatives log_pred_derivatives(const double log_B,
                                        const double catch_t,
                                        const FixedParameters &par) {
  const double B = std::exp(log_B);
  const double production = par.r * B * (1.0 - B / par.K);
  const double pred = std::max(B + production - catch_t, 1e-9);

  const double dB_dx = B;
  const double dprod_dB = par.r * (1.0 - 2.0 * B / par.K);
  const double dA_dB = 1.0 + dprod_dB;
  const double dA_dx = dA_dB * dB_dx;

  const double d2prod_dB2 = -2.0 * par.r / par.K;
  const double d2A_dx2 = d2prod_dB2 * B * B + dA_dB * B;

  LogPredDerivatives out;
  out.log_pred = std::log(pred);
  out.d1 = dA_dx / pred;
  out.d2 = d2A_dx2 / pred - (dA_dx * dA_dx) / (pred * pred);
  return out;
}

Data make_synthetic_data() {
  Data data;
  const int n = 30;
  data.catch_mt.resize(n);
  data.index.resize(n);

  for (int t = 0; t < n; ++t) {
    const double trend = 92.0 - 0.45 * static_cast<double>(t);
    const double cycle = 4.0 * std::sin(0.38 * static_cast<double>(t));
    data.catch_mt[t] = std::max(45.0, trend + cycle);
  }

  FixedParameters true_par;
  true_par.r = 0.34;
  true_par.K = 950.0;
  true_par.q = 0.00115;
  true_par.sigma_process = 0.10;
  true_par.sigma_index = 0.08;
  true_par.B0_frac = 0.82;

  double B = B0(true_par);

  for (int t = 0; t < n; ++t) {
    const double index_noise = 0.030 * std::sin(0.7 * static_cast<double>(t)) +
                               0.018 * std::cos(0.23 * static_cast<double>(t));
    data.index[t] = true_par.q * B * std::exp(index_noise);

    const double process_noise =
        0.025 * std::sin(0.49 * static_cast<double>(t + 1));
    const double production = true_par.r * B * (1.0 - B / true_par.K);
    B = std::max(B + production - data.catch_mt[t], 1e-9);
    B *= std::exp(process_noise);
  }

  return data;
}

Eigen::VectorXd deterministic_initial_x(const Data &data,
                                        const FixedParameters &par) {
  const int n_state = static_cast<int>(data.catch_mt.size()) - 1;
  Eigen::VectorXd x(n_state);
  double B = B0(par);

  for (int t = 0; t < n_state; ++t) {
    const double production = par.r * B * (1.0 - B / par.K);
    B = std::max(B + production - data.catch_mt[t], 1e-9);
    x[t] = std::log(B);
  }

  return x;
}

double joint_x(const Data &data, const FixedParameters &par,
               const Eigen::VectorXd &x) {
  const int n_state = static_cast<int>(x.size());
  const double ivp = inv_var_process(par);
  const double ivi = inv_var_index(par);

  double nll = 0.0;

  for (int k = 0; k < n_state; ++k) {
    const double log_B_prev = (k == 0) ? std::log(B0(par)) : x[k - 1];
    const auto lp = log_pred_derivatives(log_B_prev, data.catch_mt[k], par);
    const double e = x[k] - lp.log_pred;
    nll += 0.5 * e * e * ivp + std::log(par.sigma_process);
  }

  const double obs0 =
      std::log(data.index[0]) - std::log(par.q) - std::log(B0(par));
  nll += 0.5 * obs0 * obs0 * ivi + std::log(par.sigma_index);

  for (int k = 0; k < n_state; ++k) {
    const int t = k + 1;
    const double obs = std::log(data.index[t]) - std::log(par.q) - x[k];
    nll += 0.5 * obs * obs * ivi + std::log(par.sigma_index);
  }

  return nll;
}

Eigen::VectorXd grad_x(const Data &data, const FixedParameters &par,
                       const Eigen::VectorXd &x) {
  const int n_state = static_cast<int>(x.size());
  const double ivp = inv_var_process(par);
  const double ivi = inv_var_index(par);

  Eigen::VectorXd g = Eigen::VectorXd::Zero(n_state);

  for (int k = 0; k < n_state; ++k) {
    const double log_B_prev = (k == 0) ? std::log(B0(par)) : x[k - 1];
    const auto lp_k = log_pred_derivatives(log_B_prev, data.catch_mt[k], par);
    const double e_k = x[k] - lp_k.log_pred;
    g[k] += e_k * ivp;

    const int t = k + 1;
    const double obs = std::log(data.index[t]) - std::log(par.q) - x[k];
    g[k] += -obs * ivi;

    if (k + 1 < n_state) {
      const auto lp_next =
          log_pred_derivatives(x[k], data.catch_mt[k + 1], par);
      const double e_next = x[k + 1] - lp_next.log_pred;
      g[k] += -e_next * lp_next.d1 * ivp;
    }
  }

  return g;
}

void hessian_tridiagonal(const Data &data, const FixedParameters &par,
                         const Eigen::VectorXd &x, Eigen::VectorXd &diag,
                         Eigen::VectorXd &offdiag) {
  const int n_state = static_cast<int>(x.size());
  const double ivp = inv_var_process(par);
  const double ivi = inv_var_index(par);

  diag = Eigen::VectorXd::Zero(n_state);
  offdiag = Eigen::VectorXd::Zero(std::max(0, n_state - 1));

  for (int k = 0; k < n_state; ++k) {
    double d = ivp + ivi;

    if (k + 1 < n_state) {
      const auto lp_next =
          log_pred_derivatives(x[k], data.catch_mt[k + 1], par);
      const double e_next = x[k + 1] - lp_next.log_pred;
      d += (lp_next.d1 * lp_next.d1 - e_next * lp_next.d2) * ivp;
      offdiag[k] = -lp_next.d1 * ivp;
    }

    diag[k] = std::max(d, 1e-8);
  }
}

Eigen::VectorXd solve_tridiagonal_ldlt(const Eigen::VectorXd &diag,
                                       const Eigen::VectorXd &offdiag,
                                       const Eigen::VectorXd &rhs) {
  const int n = static_cast<int>(diag.size());

  Eigen::VectorXd D(n);
  Eigen::VectorXd Lsub(std::max(0, n - 1));

  D[0] = diag[0];
  for (int i = 1; i < n; ++i) {
    Lsub[i - 1] = offdiag[i - 1] / D[i - 1];
    D[i] = diag[i] - Lsub[i - 1] * offdiag[i - 1];
    D[i] = std::max(D[i], 1e-10);
  }

  Eigen::VectorXd y(n);
  y[0] = rhs[0];
  for (int i = 1; i < n; ++i) {
    y[i] = rhs[i] - Lsub[i - 1] * y[i - 1];
  }

  Eigen::VectorXd z(n);
  for (int i = 0; i < n; ++i) {
    z[i] = y[i] / D[i];
  }

  Eigen::VectorXd x(n);
  x[n - 1] = z[n - 1];
  for (int i = n - 2; i >= 0; --i) {
    x[i] = z[i] - Lsub[i] * x[i + 1];
  }

  return x;
}

double logdet_tridiagonal(const Eigen::VectorXd &diag,
                          const Eigen::VectorXd &offdiag) {
  const int n = static_cast<int>(diag.size());
  if (n == 0)
    return 0.0;

  double d_prev = diag[0];
  if (!(d_prev > 0.0))
    return std::numeric_limits<double>::quiet_NaN();

  double logdet = std::log(d_prev);
  for (int i = 1; i < n; ++i) {
    const double d = diag[i] - offdiag[i - 1] * offdiag[i - 1] / d_prev;
    if (!(d > 0.0))
      return std::numeric_limits<double>::quiet_NaN();
    logdet += std::log(d);
    d_prev = d;
  }

  return logdet;
}

NewtonResult optimize_x_newton(const Data &data, const FixedParameters &par,
                               const Eigen::VectorXd &initial_x) {
  NewtonResult out;
  out.xhat = initial_x;

  double f = joint_x(data, par, out.xhat);

  for (int iter = 0; iter < 30; ++iter) {
    const Eigen::VectorXd g = grad_x(data, par, out.xhat);
    out.grad_norm = g.norm();

    if (out.grad_norm < 1e-7) {
      out.converged = true;
      out.iterations = iter;
      out.joint = f;
      return out;
    }

    Eigen::VectorXd diag;
    Eigen::VectorXd offdiag;
    hessian_tridiagonal(data, par, out.xhat, diag, offdiag);

    const Eigen::VectorXd step = solve_tridiagonal_ldlt(diag, offdiag, -g);

    double alpha = 1.0;
    bool accepted = false;

    for (int ls = 0; ls < 30; ++ls) {
      const Eigen::VectorXd trial = out.xhat + alpha * step;
      const double f_trial = joint_x(data, par, trial);

      if (std::isfinite(f_trial) && f_trial < f) {
        out.xhat = trial;
        f = f_trial;
        accepted = true;
        break;
      }

      alpha *= 0.5;
    }

    if (!accepted) {
      out.iterations = iter + 1;
      out.joint = f;
      return out;
    }
  }

  out.joint = f;
  out.iterations = 30;
  return out;
}

double laplace_objective(const Data &data, const FixedParameters &par,
                         const Eigen::VectorXd &xhat, double *out_joint,
                         double *out_logdet) {
  Eigen::VectorXd diag;
  Eigen::VectorXd offdiag;
  hessian_tridiagonal(data, par, xhat, diag, offdiag);

  const double joint = joint_x(data, par, xhat);
  const double logdet = logdet_tridiagonal(diag, offdiag);
  const double n = static_cast<double>(xhat.size());

  if (out_joint)
    *out_joint = joint;
  if (out_logdet)
    *out_logdet = logdet;

  return joint + 0.5 * logdet - 0.5 * n * std::log(2.0 * M_PI);
}

double estimate_q_from_x(const Data &data, const FixedParameters &par_without_q,
                         const Eigen::VectorXd &xhat) {
  std::vector<double> log_biomass;
  log_biomass.reserve(data.index.size());

  log_biomass.push_back(std::log(B0(par_without_q)));
  for (int k = 0; k < xhat.size(); ++k) {
    log_biomass.push_back(xhat[k]);
  }

  double sum = 0.0;
  int n = 0;

  for (std::size_t t = 0; t < data.index.size(); ++t) {
    if (data.index[t] > 0.0 && std::isfinite(data.index[t])) {
      sum += std::log(data.index[t]) - log_biomass[t];
      ++n;
    }
  }

  if (n == 0) {
    return par_without_q.q;
  }

  return std::exp(sum / static_cast<double>(n));
}

struct OuterObjective {
  const Data &data;
  Eigen::VectorXd cached_x;
  int inner_iterations_total = 0;
  int evals = 0;

  double value(const Eigen::VectorXd &theta) {
    ++evals;

    const FixedParameters par = theta_to_par(theta);
    if (cached_x.size() == 0) {
      cached_x = deterministic_initial_x(data, par);
    }

    NewtonResult nr = optimize_x_newton(data, par, cached_x);
    cached_x = nr.xhat;
    inner_iterations_total += nr.iterations;

    return laplace_objective(data, par, nr.xhat, nullptr, nullptr);
  }

  double operator()(const Eigen::VectorXd &theta, Eigen::VectorXd &grad) {
    const double f0 = value(theta);
    grad = Eigen::VectorXd::Zero(theta.size());

    for (int j = 0; j < theta.size(); ++j) {
      const double h = 1e-4 * (1.0 + std::abs(theta[j]));
      Eigen::VectorXd tp = theta;
      Eigen::VectorXd tm = theta;
      tp[j] += h;
      tm[j] -= h;
      grad[j] = (value(tp) - value(tm)) / (2.0 * h);
    }

    return f0;
  }
};

FitResult fit_model(const Data &data) {
  // Stable public example:
  // Fit q and latent biomass states while holding biological/process
  // quantities fixed near the synthetic generating values.
  FixedParameters par;
  par.r = 0.34;
  par.K = 950.0;
  par.q = 0.0010;
  par.sigma_process = 0.10;
  par.sigma_index = 0.08;
  par.B0_frac = 0.82;

  Eigen::VectorXd x = deterministic_initial_x(data, par);

  int total_inner_iterations = 0;
  int evals = 0;

  for (int outer = 0; outer < 20; ++outer) {
    NewtonResult nr = optimize_x_newton(data, par, x);
    x = nr.xhat;
    total_inner_iterations += nr.iterations;
    ++evals;

    const double old_q = par.q;
    par.q = estimate_q_from_x(data, par, x);
    par.q = std::min(0.0050, std::max(0.0002, par.q));

    const double rel_change = std::abs(par.q - old_q) / std::max(old_q, 1e-12);

    if (rel_change < 1e-10 && nr.grad_norm < 1e-6) {
      break;
    }
  }

  NewtonResult final_nr = optimize_x_newton(data, par, x);
  total_inner_iterations += final_nr.iterations;

  double joint = 0.0;
  double logdet = 0.0;
  const double objective =
      laplace_objective(data, par, final_nr.xhat, &joint, &logdet);

  FitResult out;
  out.par = par;
  out.theta = par_to_theta(par);
  out.xhat = final_nr.xhat;
  out.objective = objective;
  out.joint = joint;
  out.logdet = logdet;
  out.grad_norm = final_nr.grad_norm;
  out.outer_iterations = evals;
  out.inner_iterations = total_inner_iterations;
  out.objective_evaluations = evals;
  return out;
}

double project_next_biomass(const double B, const double catch_mt,
                            const FixedParameters &par) {
  const double production = par.r * B * (1.0 - B / par.K);
  return std::max(B + production - catch_mt, 1e-9);
}

void write_fit_summary(const std::string &path, const FitResult &fit) {
  std::ofstream out(path);
  out << "quantity,value\n";
  out << "objective," << fit.objective << "\n";
  out << "joint," << fit.joint << "\n";
  out << "logdet," << fit.logdet << "\n";
  out << "grad_norm," << fit.grad_norm << "\n";
  out << "outer_iterations," << fit.outer_iterations << "\n";
  out << "inner_iterations," << fit.inner_iterations << "\n";
  out << "objective_evaluations," << fit.objective_evaluations << "\n";
  out << "r," << fit.par.r << "\n";
  out << "K," << fit.par.K << "\n";
  out << "q," << fit.par.q << "\n";
  out << "sigma_process," << fit.par.sigma_process << "\n";
  out << "sigma_index," << fit.par.sigma_index << "\n";
  out << "B0," << B0(fit.par) << "\n";
}

void write_fit_trajectory(const std::string &path, const Data &data,
                          const FitResult &fit) {
  std::ofstream out(path);
  out << "year,phase,catch_mt,index,biomass_hat,index_hat\n";
  out << 1 << ",history," << data.catch_mt[0] << "," << data.index[0] << ","
      << B0(fit.par) << "," << fit.par.q * B0(fit.par) << "\n";

  for (int k = 0; k < fit.xhat.size(); ++k) {
    const int year = k + 2;
    const double B = std::exp(fit.xhat[k]);
    out << year << ",history," << data.catch_mt[static_cast<std::size_t>(k + 1)]
        << "," << data.index[static_cast<std::size_t>(k + 1)] << "," << B << ","
        << fit.par.q * B << "\n";
  }
}

void write_projection_scenarios(const std::string &path, const FitResult &fit,
                                const int years) {
  const double last_B = std::exp(fit.xhat[fit.xhat.size() - 1]);

  struct Scenario {
    std::string name;
    double catch_mt;
  };

  const std::vector<Scenario> scenarios = {
      {"status_quo", 80.0},
      {"reduced_catch", 60.0},
      {"high_catch", 100.0},
  };

  std::ofstream out(path);
  out << "scenario,projection_year,catch_mt,biomass,biomass_frac_K\n";

  for (const auto &scenario : scenarios) {
    double B = last_B;
    for (int y = 1; y <= years; ++y) {
      B = project_next_biomass(B, scenario.catch_mt, fit.par);
      out << scenario.name << "," << y << "," << scenario.catch_mt << "," << B
          << "," << B / fit.par.K << "\n";
    }
  }
}

} // namespace

int main() {
  std::cout << "Synthetic opakapaka-style fit + projection example\n";
  std::cout << "==================================================\n\n";
  std::cout
      << "Synthetic and public-data-safe. Not an official assessment.\n\n";

  const Data data = make_synthetic_data();
  const FitResult fit = fit_model(data);

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "Fit summary\n";
  std::cout << "-----------\n";
  std::cout << "objective             " << fit.objective << "\n";
  std::cout << "joint                 " << fit.joint << "\n";
  std::cout << "logdet                " << fit.logdet << "\n";
  std::cout << "grad_norm             " << fit.grad_norm << "\n";
  std::cout << "outer_iterations      " << fit.outer_iterations << "\n";
  std::cout << "inner_iterations      " << fit.inner_iterations << "\n";
  std::cout << "objective_evaluations " << fit.objective_evaluations << "\n\n";

  std::cout << "Estimated parameters\n";
  std::cout << "--------------------\n";
  std::cout << "r              " << fit.par.r << " (fixed)\n";
  std::cout << "K              " << fit.par.K << " (fixed)\n";
  std::cout << "q              " << fit.par.q << " (estimated)\n";
  std::cout << "sigma_process  " << fit.par.sigma_process << " (fixed)\n";
  std::cout << "sigma_index    " << fit.par.sigma_index << " (fixed)\n";
  std::cout << "B0             " << B0(fit.par) << " (fixed fraction)\n\n";

  const std::string outdir = "examples/opakapaka_projection/outputs";
  write_fit_summary(outdir + "/synthetic_fit_summary.csv", fit);
  write_fit_trajectory(outdir + "/synthetic_fit_trajectory.csv", data, fit);
  write_projection_scenarios(outdir + "/synthetic_projection_scenarios.csv",
                             fit, 20);

  std::cout << "Wrote outputs:\n";
  std::cout << "  " << outdir << "/synthetic_fit_summary.csv\n";
  std::cout << "  " << outdir << "/synthetic_fit_trajectory.csv\n";
  std::cout << "  " << outdir << "/synthetic_projection_scenarios.csv\n";

  return 0;
}
