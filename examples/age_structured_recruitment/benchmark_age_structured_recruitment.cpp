#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using Clock = std::chrono::high_resolution_clock;

namespace {

double ms_between(const Clock::time_point &a, const Clock::time_point &b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

std::vector<int> parse_lengths(const std::string &s) {
  std::vector<int> out;
  std::stringstream ss_in(s);
  std::string item;
  while (std::getline(ss_in, item, ',')) {
    if (!item.empty())
      out.push_back(std::stoi(item));
  }
  return out;
}

struct Data {
  int n_years = 0;
  int n_ages = 0;
  std::vector<double> index_obs;
};

struct Parameters {
  double log_R0 = std::log(1000.0);
  double log_M = std::log(0.20);
  double log_q = std::log(0.001);
  double log_sigma_R = std::log(0.35);
  double log_sigma_index = std::log(0.10);
};

struct DerivedParameters {
  double R0 = 0.0;
  double M = 0.0;
  double q = 0.0;
  double sigma_R = 0.0;
  double sigma_index = 0.0;
  double inv_var_R = 0.0;
  double inv_var_index = 0.0;
};

DerivedParameters transform(const Parameters &par) {
  DerivedParameters out;
  out.R0 = std::exp(par.log_R0);
  out.M = std::exp(par.log_M);
  out.q = std::exp(par.log_q);
  out.sigma_R = std::exp(par.log_sigma_R);
  out.sigma_index = std::exp(par.log_sigma_index);
  out.inv_var_R = 1.0 / (out.sigma_R * out.sigma_R);
  out.inv_var_index = 1.0 / (out.sigma_index * out.sigma_index);
  return out;
}

double logistic(const double x) { return 1.0 / (1.0 + std::exp(-x)); }

std::vector<double> selectivity(const int n_ages) {
  std::vector<double> sel(static_cast<std::size_t>(n_ages));
  for (int a = 0; a < n_ages; ++a) {
    sel[static_cast<std::size_t>(a)] =
        logistic((static_cast<double>(a + 1) - 4.0) / 0.8);
  }
  return sel;
}

double normal_nll(const double residual, const double sigma) {
  const double z = residual / sigma;
  return std::log(sigma) + 0.5 * std::log(2.0 * M_PI) + 0.5 * z * z;
}

std::vector<double> simulate_index(const int n_years, const int n_ages,
                                   const Parameters &par) {
  const DerivedParameters p = transform(par);
  const std::vector<double> sel = selectivity(n_ages);

  std::vector<double> N(static_cast<std::size_t>(n_ages), 0.0);
  for (int a = 0; a < n_ages; ++a) {
    N[static_cast<std::size_t>(a)] =
        p.R0 * std::exp(-p.M * static_cast<double>(a));
  }

  std::vector<double> index(static_cast<std::size_t>(n_years), 0.0);

  for (int y = 0; y < n_years; ++y) {
    double vulnerable = 0.0;
    for (int a = 0; a < n_ages; ++a) {
      vulnerable +=
          sel[static_cast<std::size_t>(a)] * N[static_cast<std::size_t>(a)];
    }

    const double obs_error =
        0.04 * std::sin(2.0 * M_PI * static_cast<double>(y) / 13.0) +
        0.02 * std::cos(2.0 * M_PI * static_cast<double>(y) / 7.0);

    index[static_cast<std::size_t>(y)] = p.q * vulnerable * std::exp(obs_error);

    std::vector<double> N_next(static_cast<std::size_t>(n_ages), 0.0);
    const double recruitment_dev =
        0.15 * std::sin(2.0 * M_PI * static_cast<double>(y) / 17.0);
    N_next[0] = p.R0 * std::exp(recruitment_dev);

    for (int a = 1; a < n_ages; ++a) {
      N_next[static_cast<std::size_t>(a)] =
          N[static_cast<std::size_t>(a - 1)] * std::exp(-p.M);
    }

    // Plus group.
    N_next[static_cast<std::size_t>(n_ages - 1)] +=
        N[static_cast<std::size_t>(n_ages - 1)] * std::exp(-p.M);

    N.swap(N_next);
  }

  return index;
}

Data make_data(const int n_years, const int n_ages, const Parameters &par) {
  Data data;
  data.n_years = n_years;
  data.n_ages = n_ages;
  data.index_obs = simulate_index(n_years, n_ages, par);
  return data;
}

// x[y] = recruitment log deviation in year y.
// This simple first age-structured benchmark treats recruitment deviations as
// independent normal effects, so H_xx is diagonal. The next benchmark should
// make x a RW1/AR1 process, giving tridiagonal H_xx.
double joint_x(const Data &data, const Parameters &par,
               const Eigen::VectorXd &x) {
  const DerivedParameters p = transform(par);
  const std::vector<double> sel = selectivity(data.n_ages);

  std::vector<double> N(static_cast<std::size_t>(data.n_ages), 0.0);
  for (int a = 0; a < data.n_ages; ++a) {
    N[static_cast<std::size_t>(a)] =
        p.R0 * std::exp(-p.M * static_cast<double>(a));
  }

  double nll = 0.0;

  for (int y = 0; y < data.n_years; ++y) {
    double vulnerable = 0.0;
    for (int a = 0; a < data.n_ages; ++a) {
      vulnerable +=
          sel[static_cast<std::size_t>(a)] * N[static_cast<std::size_t>(a)];
    }

    const double log_index_pred = std::log(p.q) + std::log(vulnerable);
    const double obs_residual =
        std::log(data.index_obs[static_cast<std::size_t>(y)]) - log_index_pred;
    nll += normal_nll(obs_residual, p.sigma_index);

    const double ry = x[y];
    nll += normal_nll(ry, p.sigma_R);

    std::vector<double> N_next(static_cast<std::size_t>(data.n_ages), 0.0);
    N_next[0] = p.R0 * std::exp(ry);

    for (int a = 1; a < data.n_ages; ++a) {
      N_next[static_cast<std::size_t>(a)] =
          N[static_cast<std::size_t>(a - 1)] * std::exp(-p.M);
    }

    N_next[static_cast<std::size_t>(data.n_ages - 1)] +=
        N[static_cast<std::size_t>(data.n_ages - 1)] * std::exp(-p.M);

    N.swap(N_next);
  }

  return nll;
}

Eigen::VectorXd initial_x(const Data &data) {
  return Eigen::VectorXd::Zero(data.n_years);
}

// Finite-difference gradient is acceptable for the optimizer in this scaffold,
// but the Hessian/logdet path is analytic diagonal.
Eigen::VectorXd fd_grad_x(const Data &data, const Parameters &par,
                          const Eigen::VectorXd &x) {
  Eigen::VectorXd g(x.size());
  for (int i = 0; i < x.size(); ++i) {
    const double h = 1e-5 * (1.0 + std::abs(x[i]));
    Eigen::VectorXd xp = x;
    Eigen::VectorXd xm = x;
    xp[i] += h;
    xm[i] -= h;
    g[i] = (joint_x(data, par, xp) - joint_x(data, par, xm)) / (2.0 * h);
  }
  return g;
}

// Simple damped Newton using diagonal curvature approximation. This is fast,
// stable, and enough to test whether the model class keeps an advantage.
Eigen::VectorXd optimize_x_diagonal_newton(const Data &data,
                                           const Parameters &par) {
  const DerivedParameters p = transform(par);
  Eigen::VectorXd x = initial_x(data);

  for (int iter = 0; iter < 40; ++iter) {
    const Eigen::VectorXd g = fd_grad_x(data, par, x);
    if (g.norm() < 1e-5)
      break;

    for (int i = 0; i < x.size(); ++i) {
      // Conservative positive diagonal. Recruitment prior dominates plus a
      // contribution from index information.
      const double Hii = p.inv_var_R + 0.25 * p.inv_var_index;
      x[i] -= 0.5 * g[i] / Hii;
    }
  }

  return x;
}

// Exact diagonal prior contribution plus finite-difference local diagonal for
// observation curvature. This is still not final Quadra AD; it is a structured
// age-model scaffold.
Eigen::SparseMatrix<double>
structured_diagonal_hessian(const Data &data, const Parameters &par,
                            const Eigen::VectorXd &x) {
  const int n = static_cast<int>(x.size());
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(static_cast<std::size_t>(n));

  const DerivedParameters p = transform(par);

  for (int i = 0; i < n; ++i) {
    const double h = 1e-4 * (1.0 + std::abs(x[i]));
    Eigen::VectorXd xp = x;
    Eigen::VectorXd xm = x;
    xp[i] += h;
    xm[i] -= h;

    const double f0 = joint_x(data, par, x);
    const double fp = joint_x(data, par, xp);
    const double fm = joint_x(data, par, xm);
    double Hii = (fp - 2.0 * f0 + fm) / (h * h);

    if (!(Hii > 0.0))
      Hii = p.inv_var_R;

    triplets.emplace_back(i, i, Hii);
  }

  Eigen::SparseMatrix<double> H(n, n);
  H.setFromTriplets(triplets.begin(), triplets.end());
  return H;
}

double sparse_logdet_diagonal(const Eigen::SparseMatrix<double> &H) {
  double logdet = 0.0;
  for (int k = 0; k < H.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(H, k); it; ++it) {
      if (it.row() == it.col()) {
        if (!(it.value() > 0.0))
          throw std::runtime_error("non-positive diagonal");
        logdet += std::log(it.value());
      }
    }
  }
  return logdet;
}

struct EvalResult {
  double objective = 0.0;
  double joint = 0.0;
  double logdet = 0.0;
  int nnz = 0;
  double grad_norm = 0.0;
};

EvalResult eval(const Data &data, const Parameters &par) {
  EvalResult out;
  const Eigen::VectorXd xhat = optimize_x_diagonal_newton(data, par);

  out.joint = joint_x(data, par, xhat);
  out.grad_norm = fd_grad_x(data, par, xhat).norm();

  const Eigen::SparseMatrix<double> H =
      structured_diagonal_hessian(data, par, xhat);
  out.nnz = static_cast<int>(H.nonZeros());
  out.logdet = sparse_logdet_diagonal(H);

  const double n_x = static_cast<double>(xhat.size());
  out.objective =
      out.joint + 0.5 * out.logdet - 0.5 * n_x * std::log(2.0 * M_PI);
  return out;
}

} // namespace

int main(int argc, char **argv) {
  int reps = 10;
  std::vector<int> lengths = {25, 50, 100, 250, 500, 1000};
  int n_ages = 10;

  if (argc > 1)
    reps = std::stoi(argv[1]);
  if (argc > 2)
    lengths = parse_lengths(argv[2]);
  if (argc > 3)
    n_ages = std::stoi(argv[3]);

  const Parameters par;

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "Quadra age-structured recruitment deviation benchmark\n";
  std::cout << "reps per n = " << reps << ", ages = " << n_ages << "\n\n";

  std::cout << std::setw(8) << "n" << std::setw(14) << "objective"
            << std::setw(14) << "joint" << std::setw(14) << "logdet"
            << std::setw(14) << "nnz" << std::setw(14) << "grad_norm"
            << std::setw(14) << "avg_ms" << "\n";

  for (const int n : lengths) {
    const Data data = make_data(n, n_ages, par);
    EvalResult last = eval(data, par);

    const auto t0 = Clock::now();
    for (int r = 0; r < reps; ++r) {
      last = eval(data, par);
    }
    const auto t1 = Clock::now();

    const double avg_ms = ms_between(t0, t1) / static_cast<double>(reps);

    std::cout << std::setw(8) << n << std::setw(14) << last.objective
              << std::setw(14) << last.joint << std::setw(14) << last.logdet
              << std::setw(14) << last.nnz << std::setw(14) << last.grad_norm
              << std::setw(14) << avg_ms << "\n";
  }

  return 0;
}
