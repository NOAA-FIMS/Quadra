#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "../core/laplace/sparse_laplace_evaluation_result.hpp"

DECLARE_ADGRAPH()

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
using Clock = std::chrono::steady_clock;

struct SparseRw1Objective {
  int m;

  template <class T> T operator()(const std::vector<T> &x) const {
    const T mu = x[0];
    const T log_sigma = x[1];
    const T log_lambda0 = x[2];
    const T log_lambda_rw = x[3];
    const T log_beta = x[4];

    const T inv_sigma2 = exp(T(-2.0) * log_sigma);
    const T lambda0 = exp(log_lambda0);
    const T lambda_rw = exp(log_lambda_rw);
    const T beta = exp(log_beta);

    T nll = T(0.0);

    for (int i = 0; i < m; ++i) {
      const double xd = static_cast<double>(i + 1);
      const T y =
          T(0.6 + 0.10 * std::sin(0.21 * xd) + 0.07 * std::cos(0.47 * xd));
      const T u = x[5 + i];
      const T resid = y - mu - u;

      nll = nll + T(0.5) * resid * resid * inv_sigma2 + log_sigma +
            T(0.5 * std::log(2.0 * kPi)) + T(0.5) * lambda0 * u * u +
            beta * exp(u);
    }

    for (int i = 1; i < m; ++i) {
      const T diff = x[5 + i] - x[5 + i - 1];
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
    y[i] = 0.6 + 0.10 * std::sin(0.21 * x) + 0.07 * std::cos(0.47 * x);
  }
  return y;
}

Eigen::SparseMatrix<double> Huu_sparse_direct(const Eigen::VectorXd &theta,
                                              const Eigen::VectorXd &u) {
  const int m = static_cast<int>(u.size());
  const double inv_sigma2 = std::exp(-2.0 * theta[1]);
  const double lambda0 = std::exp(theta[2]);
  const double lambda_rw = std::exp(theta[3]);
  const double beta = std::exp(theta[4]);

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(static_cast<size_t>(3 * m - 2));

  for (int i = 0; i < m; ++i) {
    double diag = inv_sigma2 + lambda0 + beta * std::exp(u[i]);
    if (i > 0) {
      diag += lambda_rw;
    }
    if (i + 1 < m) {
      diag += lambda_rw;
    }
    triplets.emplace_back(i, i, diag);
  }

  for (int i = 1; i < m; ++i) {
    triplets.emplace_back(i, i - 1, -lambda_rw);
    triplets.emplace_back(i - 1, i, -lambda_rw);
  }

  Eigen::SparseMatrix<double> H(m, m);
  H.setFromTriplets(triplets.begin(), triplets.end());
  H.makeCompressed();
  return H;
}

Eigen::MatrixXd Huu_dense_for_mode_solve(const Eigen::VectorXd &theta,
                                         const Eigen::VectorXd &u) {
  return Eigen::MatrixXd(Huu_sparse_direct(theta, u));
}

Eigen::VectorXd random_gradient(const Eigen::VectorXd &theta,
                                const Eigen::VectorXd &u) {
  const int m = static_cast<int>(u.size());
  const double mu = theta[0];
  const double inv_sigma2 = std::exp(-2.0 * theta[1]);
  const double lambda0 = std::exp(theta[2]);
  const double lambda_rw = std::exp(theta[3]);
  const double beta = std::exp(theta[4]);
  const Eigen::VectorXd y = make_y(m);

  Eigen::VectorXd g = Eigen::VectorXd::Zero(m);

  for (int i = 0; i < m; ++i) {
    const double resid = y[i] - mu - u[i];
    g[i] += -resid * inv_sigma2 + lambda0 * u[i] + beta * std::exp(u[i]);
  }

  for (int i = 1; i < m; ++i) {
    const double diff = u[i] - u[i - 1];
    g[i] += lambda_rw * diff;
    g[i - 1] -= lambda_rw * diff;
  }

  return g;
}

Eigen::VectorXd solve_uhat(const Eigen::VectorXd &theta, int m) {
  Eigen::VectorXd u = Eigen::VectorXd::Zero(m);

  for (int iter = 0; iter < 80; ++iter) {
    const Eigen::VectorXd g = random_gradient(theta, u);
    Eigen::LDLT<Eigen::MatrixXd> ldlt(Huu_dense_for_mode_solve(theta, u));
    const Eigen::VectorXd step = ldlt.solve(g);
    u -= step;
    if (step.lpNorm<Eigen::Infinity>() < 1.0e-12) {
      break;
    }
  }

  return u;
}

double joint_objective_double(const Eigen::VectorXd &theta,
                              const Eigen::VectorXd &uhat) {
  const int m = static_cast<int>(uhat.size());
  SparseRw1Objective objective{m};

  std::vector<double> x(static_cast<size_t>(5 + m));
  for (int j = 0; j < 5; ++j) {
    x[static_cast<size_t>(j)] = theta[j];
  }
  for (int i = 0; i < m; ++i) {
    x[static_cast<size_t>(5 + i)] = uhat[i];
  }

  return objective(x);
}

Eigen::VectorXd joint_envelope_gradient(const Eigen::VectorXd &theta,
                                        const Eigen::VectorXd &uhat) {
  const int m = static_cast<int>(uhat.size());
  const double mu = theta[0];
  const double inv_sigma2 = std::exp(-2.0 * theta[1]);
  const double lambda0 = std::exp(theta[2]);
  const double lambda_rw = std::exp(theta[3]);
  const double beta = std::exp(theta[4]);
  const Eigen::VectorXd y = make_y(m);

  Eigen::VectorXd g = Eigen::VectorXd::Zero(5);

  for (int i = 0; i < m; ++i) {
    const double resid = y[i] - mu - uhat[i];
    const double expu = std::exp(uhat[i]);

    g[0] += -resid * inv_sigma2;
    g[1] += 1.0 - resid * resid * inv_sigma2;
    g[2] += 0.5 * lambda0 * uhat[i] * uhat[i];
    g[4] += beta * expu;
  }

  for (int i = 1; i < m; ++i) {
    const double diff = uhat[i] - uhat[i - 1];
    g[3] += 0.5 * lambda_rw * diff * diff;
  }

  return g;
}

Eigen::VectorXd f_u_theta_column(const Eigen::VectorXd &theta,
                                 const Eigen::VectorXd &uhat, int theta_index) {
  const int m = static_cast<int>(uhat.size());
  const double inv_sigma2 = std::exp(-2.0 * theta[1]);
  const double lambda0 = std::exp(theta[2]);
  const double lambda_rw = std::exp(theta[3]);
  const double beta = std::exp(theta[4]);
  const Eigen::VectorXd y = make_y(m);

  Eigen::VectorXd col = Eigen::VectorXd::Zero(m);

  if (theta_index == 0) {
    col.array() = inv_sigma2;
    return col;
  }
  if (theta_index == 1) {
    for (int i = 0; i < m; ++i) {
      const double resid = y[i] - theta[0] - uhat[i];
      col[i] = 2.0 * resid * inv_sigma2;
    }
    return col;
  }
  if (theta_index == 2) {
    for (int i = 0; i < m; ++i) {
      col[i] = lambda0 * uhat[i];
    }
    return col;
  }
  if (theta_index == 3) {
    for (int i = 1; i < m; ++i) {
      const double diff = uhat[i] - uhat[i - 1];
      col[i] += lambda_rw * diff;
      col[i - 1] -= lambda_rw * diff;
    }
    return col;
  }
  if (theta_index == 4) {
    for (int i = 0; i < m; ++i) {
      col[i] = beta * std::exp(uhat[i]);
    }
    return col;
  }

  return col;
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

struct EvalResult {
  double ms = 0.0;
  Eigen::VectorXd gradient;
};

EvalResult evaluate_active_set(int m, int reps,
                               const std::vector<int> &active_directions) {
  Eigen::VectorXd theta(5);
  theta << 0.55, std::log(0.65), std::log(0.55), std::log(0.90), std::log(0.25);

  const Eigen::VectorXd uhat = solve_uhat(theta, m);
  const double joint = joint_objective_double(theta, uhat);
  const Eigen::VectorXd gj = joint_envelope_gradient(theta, uhat);
  const Eigen::SparseMatrix<double> H = Huu_sparse_direct(theta, uhat);
  const auto pattern = tridiagonal_pattern(m);

  SparseRw1Objective objective{m};

  auto cross = [theta, uhat](int theta_index) -> Eigen::VectorXd {
    return f_u_theta_column(theta, uhat, theta_index);
  };

  quadra::laplace::SparseLaplaceExactGradientEvaluationInputs inputs;
  inputs.theta = theta;
  inputs.uhat = uhat;
  inputs.Huu = H;
  inputs.joint_objective = joint;
  inputs.joint_envelope_gradient = gj;

  EvalResult out;

  const auto once =
      quadra::laplace::evaluate_sparse_laplace_with_exact_gradient(
          objective, cross, 5, m, pattern, active_directions, inputs);
  out.gradient = once.gradient;

  out.ms = mean_ms(
      [&]() {
        const auto result =
            quadra::laplace::evaluate_sparse_laplace_with_exact_gradient(
                objective, cross, 5, m, pattern, active_directions, inputs);
        return result.objective + result.gradient.sum();
      },
      reps);

  return out;
}

struct Row {
  double all_ms = 0.0;
  double no_mu_ms = 0.0;
  double structural_ms = 0.0;
  double no_mu_speedup = 0.0;
  double structural_speedup = 0.0;
  double no_mu_diff = 0.0;
  double structural_diff = 0.0;
};

Row run_case(int m, int reps) {
  const auto all =
      evaluate_active_set(m, reps, std::vector<int>{0, 1, 2, 3, 4});

  const auto no_mu = evaluate_active_set(m, reps, std::vector<int>{1, 2, 3, 4});

  const auto structural =
      evaluate_active_set(m, reps, std::vector<int>{2, 3, 4});

  Row out;
  out.all_ms = all.ms;
  out.no_mu_ms = no_mu.ms;
  out.structural_ms = structural.ms;
  out.no_mu_speedup = out.all_ms / out.no_mu_ms;
  out.structural_speedup = out.all_ms / out.structural_ms;
  out.no_mu_diff = (all.gradient - no_mu.gradient).cwiseAbs().maxCoeff();
  out.structural_diff =
      (all.gradient - structural.gradient).cwiseAbs().maxCoeff();

  return out;
}

} // namespace

int main(int argc, char **argv) {
  int reps = 20;
  if (argc > 1) {
    reps = std::stoi(argv[1]);
  }

  const std::vector<int> dims = {10, 25, 50, 100, 250, 500};

  std::cout << "Sparse RW1 active-direction pruning benchmark\n";
  std::cout << "reps per case = " << reps << "\n";
  std::cout
      << "theta = [mu, log_sigma, log_lambda0, log_lambda_rw, log_beta]\n\n";

  std::cout << std::setw(8) << "m" << std::setw(14) << "all ms" << std::setw(14)
            << "no_mu ms" << std::setw(16) << "struct ms" << std::setw(14)
            << "no_mu spd" << std::setw(14) << "struct spd" << std::setw(16)
            << "no_mu diff" << std::setw(16) << "struct diff" << "\n";

  std::cout << std::scientific << std::setprecision(6);

  for (int m : dims) {
    const Row r = run_case(m, reps);

    std::cout << std::setw(8) << m << std::setw(14) << r.all_ms << std::setw(14)
              << r.no_mu_ms << std::setw(16) << r.structural_ms << std::setw(14)
              << r.no_mu_speedup << std::setw(14) << r.structural_speedup
              << std::setw(16) << r.no_mu_diff << std::setw(16)
              << r.structural_diff << "\n";
  }

  std::cout << "\nBenchmark complete.\n";
  return 0;
}
