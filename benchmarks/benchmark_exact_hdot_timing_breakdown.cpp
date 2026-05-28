#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../core/laplace/full_exact_laplace_gradient_hdot.hpp"
#include "../core/laplace/had_quadra_sparse_exact_hdot_provider.hpp"

DECLARE_ADGRAPH()

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
using Clock = std::chrono::steady_clock;

struct Result {
  double hdot_all_dirs_ms = 0.0;
  double trace_all_dirs_ms = 0.0;
  double full_gradient_ms = 0.0;
  double hdot_fraction = 0.0;
  double trace_fraction = 0.0;
  int pattern_nnz = 0;
};

struct BandedObjective {
  int m;

  template <class T> T operator()(const std::vector<T> &x) const {
    const T mu = x[0];
    const T log_sigma = x[1];
    const T log_lambda0 = x[2];
    const T log_lambda_rw = x[3];

    const T inv_sigma2 = exp(T(-2.0) * log_sigma);
    const T lambda0 = exp(log_lambda0);
    const T lambda_rw = exp(log_lambda_rw);

    T nll = T(0.0);

    for (int i = 0; i < m; ++i) {
      const double xd = static_cast<double>(i + 1);
      const T y =
          T(1.35 + 0.15 * std::sin(0.7 * xd) + 0.03 * std::cos(1.3 * xd));
      const T u = x[4 + i];
      const T resid = y - mu - u;

      nll = nll + T(0.5) * resid * resid * inv_sigma2 + log_sigma +
            T(0.5 * std::log(2.0 * kPi));

      nll = nll + T(0.5) * lambda0 * u * u;
    }

    for (int i = 1; i < m; ++i) {
      const T ui = x[4 + i];
      const T uim1 = x[4 + i - 1];
      const T diff = ui - uim1;
      nll = nll + T(0.5) * lambda_rw * diff * diff;
    }

    return nll;
  }
};

quadra::laplace::RandomHessianPattern tridiagonal_pattern(int m) {
  quadra::laplace::RandomHessianPattern pattern;
  pattern.reserve(static_cast<size_t>(2 * m - 1));

  for (int i = 0; i < m; ++i) {
    pattern.emplace_back(i, i);
    if (i > 0) {
      pattern.emplace_back(i, i - 1);
    }
  }

  return pattern;
}

Eigen::VectorXd make_y(int m) {
  Eigen::VectorXd y(m);
  for (int i = 0; i < m; ++i) {
    const double x = static_cast<double>(i + 1);
    y[i] = 1.35 + 0.15 * std::sin(0.7 * x) + 0.03 * std::cos(1.3 * x);
  }
  return y;
}

Eigen::MatrixXd Huu_dense(const Eigen::VectorXd &theta, int m) {
  const double inv_sigma2 = std::exp(-2.0 * theta[1]);
  const double lambda0 = std::exp(theta[2]);
  const double lambda_rw = std::exp(theta[3]);

  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(m, m);

  for (int i = 0; i < m; ++i) {
    H(i, i) += inv_sigma2 + lambda0;
  }

  for (int i = 1; i < m; ++i) {
    H(i, i) += lambda_rw;
    H(i - 1, i - 1) += lambda_rw;
    H(i, i - 1) -= lambda_rw;
    H(i - 1, i) -= lambda_rw;
  }

  return H;
}

Eigen::VectorXd uhat_for(const Eigen::VectorXd &theta, int m) {
  const double mu = theta[0];
  const double inv_sigma2 = std::exp(-2.0 * theta[1]);
  const Eigen::VectorXd y = make_y(m);
  const Eigen::MatrixXd H = Huu_dense(theta, m);

  Eigen::VectorXd rhs(m);
  for (int i = 0; i < m; ++i) {
    rhs[i] = inv_sigma2 * (y[i] - mu);
  }

  Eigen::LDLT<Eigen::MatrixXd> ldlt(H);
  return ldlt.solve(rhs);
}

Eigen::VectorXd joint_grad(const Eigen::VectorXd &theta,
                           const Eigen::VectorXd &uhat) {
  const int m = static_cast<int>(uhat.size());
  const double mu = theta[0];
  const double inv_sigma2 = std::exp(-2.0 * theta[1]);
  const double lambda0 = std::exp(theta[2]);
  const double lambda_rw = std::exp(theta[3]);
  const Eigen::VectorXd y = make_y(m);

  Eigen::VectorXd g(4);
  g.setZero();

  for (int i = 0; i < m; ++i) {
    const double resid = y[i] - mu - uhat[i];
    g[0] += -resid * inv_sigma2;
    g[1] += 1.0 - resid * resid * inv_sigma2;
    g[2] += 0.5 * lambda0 * uhat[i] * uhat[i];
  }

  for (int i = 1; i < m; ++i) {
    const double diff = uhat[i] - uhat[i - 1];
    g[3] += 0.5 * lambda_rw * diff * diff;
  }

  return g;
}

template <class Fn> double mean_ms(Fn &&fn, int reps) {
  const auto start = Clock::now();
  for (int r = 0; r < reps; ++r) {
    volatile double sink = fn();
    (void)sink;
  }
  const auto end = Clock::now();
  std::chrono::duration<double, std::milli> elapsed = end - start;
  return elapsed.count() / static_cast<double>(reps);
}

double trace_dense_solve(const Eigen::MatrixXd &H,
                         const Eigen::SparseMatrix<double> &Hdot) {
  Eigen::LDLT<Eigen::MatrixXd> ldlt(H);
  const Eigen::MatrixXd rhs = Eigen::MatrixXd(Hdot);
  return ldlt.solve(rhs).trace();
}

Result run_case(int m, int reps) {
  Eigen::VectorXd theta(4);
  theta << 1.1, std::log(0.45), std::log(0.75), std::log(1.35);

  const Eigen::VectorXd uhat = uhat_for(theta, m);
  const Eigen::VectorXd gj = joint_grad(theta, uhat);
  const Eigen::MatrixXd H = Huu_dense(theta, m);

  BandedObjective objective{m};

  auto sparse_provider =
      quadra::laplace::make_had_quadra_sparse_exact_hdot_provider(
          objective, 4, m, tridiagonal_pattern(m), 0.0);

  std::vector<Eigen::SparseMatrix<double>> Hdots;
  Hdots.reserve(4);

  for (int j = 0; j < 4; ++j) {
    Hdots.push_back(sparse_provider.sparse(theta, uhat, j));
  }

  auto Hfn = [](const Eigen::VectorXd &th, const Eigen::VectorXd &uh) {
    return Huu_dense(th, static_cast<int>(uh.size()));
  };

  Result out;
  out.pattern_nnz = 2 * m - 1;

  out.hdot_all_dirs_ms = mean_ms(
      [&]() {
        double s = 0.0;
        for (int j = 0; j < 4; ++j) {
          Eigen::SparseMatrix<double> Hdot =
              sparse_provider.sparse(theta, uhat, j);
          s += static_cast<double>(Hdot.nonZeros());
        }
        return s;
      },
      reps);

  out.trace_all_dirs_ms = mean_ms(
      [&]() {
        double s = 0.0;
        for (int j = 0; j < 4; ++j) {
          s += trace_dense_solve(H, Hdots[static_cast<size_t>(j)]);
        }
        return s;
      },
      reps);

  out.full_gradient_ms = mean_ms(
      [&]() {
        return quadra::laplace::full_exact_laplace_gradient_hdot(
                   gj, Hfn, sparse_provider, theta, uhat)
            .sum();
      },
      reps);

  if (out.full_gradient_ms > 0.0) {
    out.hdot_fraction = out.hdot_all_dirs_ms / out.full_gradient_ms;
    out.trace_fraction = out.trace_all_dirs_ms / out.full_gradient_ms;
  }

  return out;
}

} // namespace

int main(int argc, char **argv) {
  int reps = 10;
  if (argc > 1) {
    reps = std::stoi(argv[1]);
  }

  const std::vector<int> dims = {10, 25, 50, 100, 250};

  std::cout << "Exact sparse Hdot timing breakdown benchmark\n";
  std::cout << "reps per case = " << reps << "\n\n";

  std::cout << std::setw(8) << "m" << std::setw(12) << "pattern"
            << std::setw(20) << "Hdot all dirs ms" << std::setw(20)
            << "trace all dirs ms" << std::setw(20) << "full grad ms"
            << std::setw(16) << "Hdot/full" << std::setw(16) << "trace/full"
            << "\n";

  std::cout << std::scientific << std::setprecision(6);

  for (int m : dims) {
    const Result r = run_case(m, reps);

    std::cout << std::setw(8) << m << std::setw(12) << r.pattern_nnz
              << std::setw(20) << r.hdot_all_dirs_ms << std::setw(20)
              << r.trace_all_dirs_ms << std::setw(20) << r.full_gradient_ms
              << std::setw(16) << r.hdot_fraction << std::setw(16)
              << r.trace_fraction << "\n";
  }

  std::cout << "\nBenchmark complete.\n";
  return 0;
}
