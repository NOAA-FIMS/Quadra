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

double logistic(const double x) { return 1.0 / (1.0 + std::exp(-x)); }

double normal_nll(const double residual, const double sigma) {
  const double z = residual / sigma;
  return std::log(sigma) + 0.5 * std::log(2.0 * M_PI) + 0.5 * z * z;
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

struct Derived {
  double R0 = 0.0;
  double M = 0.0;
  double q = 0.0;
  double sigma_R = 0.0;
  double sigma_index = 0.0;
};

Derived transform(const Parameters &par) {
  Derived d;
  d.R0 = std::exp(par.log_R0);
  d.M = std::exp(par.log_M);
  d.q = std::exp(par.log_q);
  d.sigma_R = std::exp(par.log_sigma_R);
  d.sigma_index = std::exp(par.log_sigma_index);
  return d;
}

std::vector<double> selectivity(const int n_ages) {
  std::vector<double> sel(static_cast<std::size_t>(n_ages));
  for (int a = 0; a < n_ages; ++a) {
    sel[static_cast<std::size_t>(a)] =
        logistic((static_cast<double>(a + 1) - 4.0) / 0.8);
  }
  return sel;
}

std::vector<double> simulate_index(const int n_years, const int n_ages,
                                   const Parameters &par) {
  const Derived p = transform(par);
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

double joint_x(const Data &data, const Parameters &par,
               const Eigen::VectorXd &x) {
  const Derived p = transform(par);
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

    const double obs_residual =
        std::log(data.index_obs[static_cast<std::size_t>(y)]) -
        (std::log(p.q) + std::log(vulnerable));

    nll += normal_nll(obs_residual, p.sigma_index);
    nll += normal_nll(x[y], p.sigma_R);

    std::vector<double> N_next(static_cast<std::size_t>(data.n_ages), 0.0);
    N_next[0] = p.R0 * std::exp(x[y]);

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

class ObjX {
public:
  ObjX(const Data &data, const Parameters &par) : data_(data), par_(par) {}

  double operator()(const Eigen::VectorXd &x, Eigen::VectorXd &grad) {
    const double f = joint_x(data_, par_, x);
    grad = fd_grad_x(data_, par_, x);
    return f;
  }

private:
  const Data &data_;
  const Parameters &par_;
};

Eigen::VectorXd optimize_x(const Data &data, const Parameters &par) {
  Eigen::VectorXd x = Eigen::VectorXd::Zero(data.n_years);

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
    if (!(std::isfinite(joint_x(data, par, x)) && gnorm < 1e-3)) {
      throw;
    }
  }

  return x;
}

Eigen::MatrixXd fd_dense_hessian_xx(const Data &data, const Parameters &par,
                                    const Eigen::VectorXd &x) {
  const int n = static_cast<int>(x.size());
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(n, n);
  const double f0 = joint_x(data, par, x);

  for (int i = 0; i < n; ++i) {
    const double hi = 1e-4 * (1.0 + std::abs(x[i]));

    Eigen::VectorXd xp = x;
    Eigen::VectorXd xm = x;
    xp[i] += hi;
    xm[i] -= hi;

    H(i, i) = (joint_x(data, par, xp) - 2.0 * f0 + joint_x(data, par, xm)) /
              (hi * hi);

    for (int j = i + 1; j < n; ++j) {
      const double hj = 1e-4 * (1.0 + std::abs(x[j]));

      Eigen::VectorXd xpp = x;
      Eigen::VectorXd xpm = x;
      Eigen::VectorXd xmp = x;
      Eigen::VectorXd xmm = x;

      xpp[i] += hi;
      xpp[j] += hj;
      xpm[i] += hi;
      xpm[j] -= hj;
      xmp[i] -= hi;
      xmp[j] += hj;
      xmm[i] -= hi;
      xmm[j] -= hj;

      const double hij = (joint_x(data, par, xpp) - joint_x(data, par, xpm) -
                          joint_x(data, par, xmp) + joint_x(data, par, xmm)) /
                         (4.0 * hi * hj);

      H(i, j) = hij;
      H(j, i) = hij;
    }
  }

  return H;
}

Eigen::SparseMatrix<double> fd_banded_hessian_xx(const Data &data,
                                                 const Parameters &par,
                                                 const Eigen::VectorXd &x,
                                                 const int bandwidth) {
  const int n = static_cast<int>(x.size());
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(static_cast<std::size_t>(n * (2 * bandwidth + 1)));

  const double f0 = joint_x(data, par, x);

  for (int i = 0; i < n; ++i) {
    const double hi = 1e-4 * (1.0 + std::abs(x[i]));

    Eigen::VectorXd xp = x;
    Eigen::VectorXd xm = x;
    xp[i] += hi;
    xm[i] -= hi;

    const double hii =
        (joint_x(data, par, xp) - 2.0 * f0 + joint_x(data, par, xm)) /
        (hi * hi);

    if (std::abs(hii) > 1e-12) {
      triplets.emplace_back(i, i, hii);
    }

    const int jmax = std::min(n - 1, i + bandwidth);
    for (int j = i + 1; j <= jmax; ++j) {
      const double hj = 1e-4 * (1.0 + std::abs(x[j]));

      Eigen::VectorXd xpp = x;
      Eigen::VectorXd xpm = x;
      Eigen::VectorXd xmp = x;
      Eigen::VectorXd xmm = x;

      xpp[i] += hi;
      xpp[j] += hj;
      xpm[i] += hi;
      xpm[j] -= hj;
      xmp[i] -= hi;
      xmp[j] += hj;
      xmm[i] -= hi;
      xmm[j] -= hj;

      const double hij = (joint_x(data, par, xpp) - joint_x(data, par, xpm) -
                          joint_x(data, par, xmp) + joint_x(data, par, xmm)) /
                         (4.0 * hi * hj);

      if (std::abs(hij) > 1e-12) {
        triplets.emplace_back(i, j, hij);
        triplets.emplace_back(j, i, hij);
      }
    }
  }

  Eigen::SparseMatrix<double> H(n, n);
  H.setFromTriplets(triplets.begin(), triplets.end());
  return H;
}

double dense_logdet_cholesky(const Eigen::MatrixXd &H) {
  Eigen::LLT<Eigen::MatrixXd> llt(H);
  if (llt.info() != Eigen::Success) {
    throw std::runtime_error("dense Hxx is not SPD");
  }

  const auto &L = llt.matrixL();
  double logdet = 0.0;
  for (int i = 0; i < H.rows(); ++i) {
    logdet += 2.0 * std::log(L(i, i));
  }

  return logdet;
}

double sparse_logdet_ldlt(const Eigen::SparseMatrix<double> &H) {
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
  ldlt.compute(H);

  if (ldlt.info() != Eigen::Success) {
    throw std::runtime_error("sparse LDLT failed");
  }

  const auto &D = ldlt.vectorD();

  double logdet = 0.0;
  for (int i = 0; i < D.size(); ++i) {
    if (!(D[i] > 0.0)) {
      throw std::runtime_error("sparse Hxx is not positive definite");
    }
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

double laplace_from_parts(const double joint, const double logdet,
                          const int n_random) {
  return joint + 0.5 * logdet -
         0.5 * static_cast<double>(n_random) * std::log(2.0 * M_PI);
}

EvalResult eval_dense(const Data &data, const Parameters &par,
                      const Eigen::VectorXd &xhat) {
  EvalResult out;
  out.joint = joint_x(data, par, xhat);
  out.grad_norm = fd_grad_x(data, par, xhat).norm();

  const Eigen::MatrixXd H = fd_dense_hessian_xx(data, par, xhat);
  out.nnz = static_cast<int>(H.rows() * H.cols());
  out.logdet = dense_logdet_cholesky(H);
  out.objective =
      laplace_from_parts(out.joint, out.logdet, static_cast<int>(xhat.size()));
  return out;
}

EvalResult eval_banded(const Data &data, const Parameters &par,
                       const Eigen::VectorXd &xhat, const int bandwidth) {
  EvalResult out;
  out.joint = joint_x(data, par, xhat);
  out.grad_norm = fd_grad_x(data, par, xhat).norm();

  const Eigen::SparseMatrix<double> H =
      fd_banded_hessian_xx(data, par, xhat, bandwidth);
  out.nnz = static_cast<int>(H.nonZeros());
  out.logdet = sparse_logdet_ldlt(H);
  out.objective =
      laplace_from_parts(out.joint, out.logdet, static_cast<int>(xhat.size()));
  return out;
}

} // namespace

int main(int argc, char **argv) {
  int reps = 3;
  std::vector<int> lengths = {25, 50};
  int n_ages = 10;
  int bandwidth = 9;

  if (argc > 1)
    reps = std::stoi(argv[1]);
  if (argc > 2)
    lengths = parse_lengths(argv[2]);
  if (argc > 3)
    n_ages = std::stoi(argv[3]);
  if (argc > 4)
    bandwidth = std::stoi(argv[4]);

  const Parameters par;

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "Quadra age-structured banded Laplace validation\n";
  std::cout << "reps per n = " << reps << ", ages = " << n_ages
            << ", bandwidth = " << bandwidth << "\n\n";

  std::cout << std::setw(8) << "n" << std::setw(14) << "dense_obj"
            << std::setw(14) << "band_obj" << std::setw(14) << "diff"
            << std::setw(14) << "band_nnz" << std::setw(14) << "grad_norm"
            << std::setw(14) << "band_ms" << "\n";

  for (const int n : lengths) {
    const Data data = make_data(n, n_ages, par);
    const Eigen::VectorXd xhat = optimize_x(data, par);

    const EvalResult dense = eval_dense(data, par, xhat);
    EvalResult band = eval_banded(data, par, xhat, bandwidth);

    const auto t0 = Clock::now();
    for (int r = 0; r < reps; ++r) {
      band = eval_banded(data, par, xhat, bandwidth);
    }
    const auto t1 = Clock::now();

    const double band_ms = ms_between(t0, t1) / static_cast<double>(reps);

    std::cout << std::setw(8) << n << std::setw(14) << dense.objective
              << std::setw(14) << band.objective << std::setw(14)
              << (band.objective - dense.objective) << std::setw(14) << band.nnz
              << std::setw(14) << band.grad_norm << std::setw(14) << band_ms
              << "\n";
  }

  return 0;
}
