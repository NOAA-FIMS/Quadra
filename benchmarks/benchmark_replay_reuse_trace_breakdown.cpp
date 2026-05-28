#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../core/laplace/had_quadra_replay_reuse_sparse_hdot_provider.hpp"

DECLARE_ADGRAPH()

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
using Clock = std::chrono::steady_clock;

struct Result {
  double hdot_ms = 0.0;
  double factor_once_trace_ms = 0.0;
  double factor_each_trace_ms = 0.0;
  double full_gradient_ms = 0.0;
  double hdot_full = 0.0;
  double trace_full = 0.0;
  double factor_reuse_speedup = 0.0;
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
            T(0.5 * std::log(2.0 * kPi)) + T(0.5) * lambda0 * u * u;
    }

    for (int i = 1; i < m; ++i) {
      const T diff = x[4 + i] - x[4 + i - 1];
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

Eigen::MatrixXd Huu(const Eigen::VectorXd &theta, const Eigen::VectorXd &uhat) {
  const int m = static_cast<int>(uhat.size());
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
  const Eigen::MatrixXd H = Huu(theta, Eigen::VectorXd::Zero(m));

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

double
trace_factor_once(const Eigen::MatrixXd &H,
                  const std::vector<Eigen::SparseMatrix<double>> &Hdots) {
  Eigen::LDLT<Eigen::MatrixXd> ldlt(H);
  if (ldlt.info() != Eigen::Success) {
    throw std::runtime_error("LDLT failed.");
  }

  double s = 0.0;
  for (const auto &Hdot : Hdots) {
    if (Hdot.nonZeros() == 0) {
      continue;
    }
    const Eigen::MatrixXd rhs = Eigen::MatrixXd(Hdot);
    s += ldlt.solve(rhs).trace();
  }
  return s;
}

double
trace_factor_each(const Eigen::MatrixXd &H,
                  const std::vector<Eigen::SparseMatrix<double>> &Hdots) {
  double s = 0.0;
  for (const auto &Hdot : Hdots) {
    if (Hdot.nonZeros() == 0) {
      continue;
    }
    Eigen::LDLT<Eigen::MatrixXd> ldlt(H);
    if (ldlt.info() != Eigen::Success) {
      throw std::runtime_error("LDLT failed.");
    }
    const Eigen::MatrixXd rhs = Eigen::MatrixXd(Hdot);
    s += ldlt.solve(rhs).trace();
  }
  return s;
}

Result run_case(int m, int reps) {
  Eigen::VectorXd theta(4);
  theta << 1.1, std::log(0.45), std::log(0.75), std::log(1.35);

  const Eigen::VectorXd uhat = uhat_for(theta, m);
  const Eigen::VectorXd gj = joint_grad(theta, uhat);
  const Eigen::MatrixXd H = Huu(theta, uhat);

  BandedObjective objective{m};

  auto provider =
      quadra::laplace::make_had_quadra_replay_reuse_sparse_hdot_provider(
          objective, 4, m, tridiagonal_pattern(m), std::vector<int>{1, 2, 3},
          0.0);

  std::vector<Eigen::SparseMatrix<double>> Hdots =
      provider.compute_all_sparse(theta, uhat);

  auto Hfn = [](const Eigen::VectorXd &th, const Eigen::VectorXd &uh) {
    return Huu(th, uh);
  };

  Result out;

  out.hdot_ms = mean_ms(
      [&]() {
        const auto tmp = provider.compute_all_sparse(theta, uhat);
        double s = 0.0;
        for (const auto &Hdot : tmp) {
          s += static_cast<double>(Hdot.nonZeros());
        }
        return s;
      },
      reps);

  out.factor_once_trace_ms =
      mean_ms([&]() { return trace_factor_once(H, Hdots); }, reps);

  out.factor_each_trace_ms =
      mean_ms([&]() { return trace_factor_each(H, Hdots); }, reps);

  out.full_gradient_ms = mean_ms(
      [&]() {
        return quadra::laplace::full_exact_laplace_gradient_replay_reuse_hdot(
                   gj, Hfn, provider, theta, uhat)
            .sum();
      },
      reps);

  out.hdot_full = out.hdot_ms / out.full_gradient_ms;
  out.trace_full = out.factor_once_trace_ms / out.full_gradient_ms;
  out.factor_reuse_speedup =
      out.factor_each_trace_ms / out.factor_once_trace_ms;

  return out;
}

} // namespace

int main(int argc, char **argv) {
  int reps = 20;
  if (argc > 1) {
    reps = std::stoi(argv[1]);
  }

  const std::vector<int> dims = {10, 25, 50, 100, 250};

  std::cout << "Replay-reuse trace breakdown benchmark\n";
  std::cout << "reps per case = " << reps << "\n\n";

  std::cout << std::setw(8) << "m" << std::setw(18) << "Hdot ms"
            << std::setw(20) << "trace once ms" << std::setw(20)
            << "trace each ms" << std::setw(18) << "full grad ms"
            << std::setw(14) << "Hdot/full" << std::setw(14) << "trace/full"
            << std::setw(18) << "factor speedup" << "\n";

  std::cout << std::scientific << std::setprecision(6);

  for (int m : dims) {
    const Result r = run_case(m, reps);

    std::cout << std::setw(8) << m << std::setw(18) << r.hdot_ms
              << std::setw(20) << r.factor_once_trace_ms << std::setw(20)
              << r.factor_each_trace_ms << std::setw(18) << r.full_gradient_ms
              << std::setw(14) << r.hdot_full << std::setw(14) << r.trace_full
              << std::setw(18) << r.factor_reuse_speedup << "\n";
  }

  std::cout << "\nBenchmark complete.\n";
  return 0;
}
