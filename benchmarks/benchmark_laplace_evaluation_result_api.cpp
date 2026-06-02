#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../core/laplace/laplace_evaluation_result.hpp"

DECLARE_ADGRAPH()

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
using Clock = std::chrono::steady_clock;

struct NonlinearDiagObjective {
  int m;

  template <class T> T operator()(const std::vector<T> &x) const {
    const T mu = x[0];
    const T log_sigma = x[1];
    const T log_lambda0 = x[2];
    const T log_beta = x[3];

    const T inv_sigma2 = exp(T(-2.0) * log_sigma);
    const T lambda0 = exp(log_lambda0);
    const T beta = exp(log_beta);

    T nll = T(0.0);

    for (int i = 0; i < m; ++i) {
      const double xd = static_cast<double>(i + 1);
      const T y =
          T(0.8 + 0.08 * std::sin(0.4 * xd) + 0.05 * std::cos(0.9 * xd));
      const T u = x[4 + i];
      const T resid = y - mu - u;

      nll = nll + T(0.5) * resid * resid * inv_sigma2 + log_sigma +
            T(0.5 * std::log(2.0 * kPi)) + T(0.5) * lambda0 * u * u +
            beta * exp(u);
    }

    return nll;
  }
};

quadra::laplace::RandomHessianPattern diagonal_pattern(int m) {
  quadra::laplace::RandomHessianPattern pattern;
  for (int i = 0; i < m; ++i) {
    pattern.emplace_back(i, i);
  }
  return pattern;
}

Eigen::VectorXd make_y(int m) {
  Eigen::VectorXd y(m);
  for (int i = 0; i < m; ++i) {
    const double x = static_cast<double>(i + 1);
    y[i] = 0.8 + 0.08 * std::sin(0.4 * x) + 0.05 * std::cos(0.9 * x);
  }
  return y;
}

Eigen::VectorXd solve_uhat(const Eigen::VectorXd &theta, int m) {
  const double mu = theta[0];
  const double inv_sigma2 = std::exp(-2.0 * theta[1]);
  const double lambda0 = std::exp(theta[2]);
  const double beta = std::exp(theta[3]);
  const Eigen::VectorXd y = make_y(m);

  Eigen::VectorXd u = Eigen::VectorXd::Zero(m);

  for (int iter = 0; iter < 80; ++iter) {
    double max_step = 0.0;
    for (int i = 0; i < m; ++i) {
      const double resid = y[i] - mu - u[i];
      const double expu = std::exp(u[i]);
      const double g = -resid * inv_sigma2 + lambda0 * u[i] + beta * expu;
      const double H = inv_sigma2 + lambda0 + beta * expu;
      const double step = g / H;
      u[i] -= step;
      max_step = std::max(max_step, std::abs(step));
    }
    if (max_step < 1.0e-13) {
      break;
    }
  }

  return u;
}

Eigen::MatrixXd Huu(const Eigen::VectorXd &theta, const Eigen::VectorXd &uhat) {
  const int m = static_cast<int>(uhat.size());
  const double inv_sigma2 = std::exp(-2.0 * theta[1]);
  const double lambda0 = std::exp(theta[2]);
  const double beta = std::exp(theta[3]);

  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(m, m);
  for (int i = 0; i < m; ++i) {
    H(i, i) = inv_sigma2 + lambda0 + beta * std::exp(uhat[i]);
  }
  return H;
}

Eigen::VectorXd f_u_theta_column(const Eigen::VectorXd &theta,
                                 const Eigen::VectorXd &uhat, int theta_index) {
  const int m = static_cast<int>(uhat.size());
  const double inv_sigma2 = std::exp(-2.0 * theta[1]);
  const double lambda0 = std::exp(theta[2]);
  const double beta = std::exp(theta[3]);
  const Eigen::VectorXd y = make_y(m);

  Eigen::VectorXd col(m);
  col.setZero();

  for (int i = 0; i < m; ++i) {
    const double resid = y[i] - theta[0] - uhat[i];
    const double expu = std::exp(uhat[i]);

    if (theta_index == 0) {
      col[i] = inv_sigma2;
    } else if (theta_index == 1) {
      col[i] = 2.0 * resid * inv_sigma2;
    } else if (theta_index == 2) {
      col[i] = lambda0 * uhat[i];
    } else if (theta_index == 3) {
      col[i] = beta * expu;
    }
  }

  return col;
}

Eigen::VectorXd joint_envelope_gradient(const Eigen::VectorXd &theta,
                                        const Eigen::VectorXd &uhat) {
  const int m = static_cast<int>(uhat.size());
  const double mu = theta[0];
  const double inv_sigma2 = std::exp(-2.0 * theta[1]);
  const double lambda0 = std::exp(theta[2]);
  const double beta = std::exp(theta[3]);
  const Eigen::VectorXd y = make_y(m);

  Eigen::VectorXd g(4);
  g.setZero();

  for (int i = 0; i < m; ++i) {
    const double resid = y[i] - mu - uhat[i];
    const double expu = std::exp(uhat[i]);

    g[0] += -resid * inv_sigma2;
    g[1] += 1.0 - resid * resid * inv_sigma2;
    g[2] += 0.5 * lambda0 * uhat[i] * uhat[i];
    g[3] += beta * expu;
  }

  return g;
}

double joint_objective_double(const Eigen::VectorXd &theta,
                              const Eigen::VectorXd &uhat) {
  const int m = static_cast<int>(uhat.size());
  NonlinearDiagObjective objective{m};

  std::vector<double> x(static_cast<size_t>(4 + m));
  for (int j = 0; j < 4; ++j) {
    x[static_cast<size_t>(j)] = theta[j];
  }
  for (int i = 0; i < m; ++i) {
    x[static_cast<size_t>(4 + i)] = uhat[i];
  }

  return objective(x);
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

struct Result {
  double eval_ms = 0.0;
  double objective = 0.0;
  double grad_norm = 0.0;
};

Result run_case(int m, int reps) {
  Eigen::VectorXd theta(4);
  theta << 0.65, std::log(0.55), std::log(0.85), std::log(0.30);

  const Eigen::VectorXd uhat = solve_uhat(theta, m);
  const Eigen::MatrixXd H = Huu(theta, uhat);
  const Eigen::VectorXd gj = joint_envelope_gradient(theta, uhat);
  const double joint = joint_objective_double(theta, uhat);

  NonlinearDiagObjective objective{m};

  auto cross = [theta, uhat](int theta_index) -> Eigen::VectorXd {
    return f_u_theta_column(theta, uhat, theta_index);
  };

  Result out;

  out.eval_ms = mean_ms(
      [&]() {
        const auto result =
            quadra::laplace::evaluate_laplace_objective_and_gradient(
                objective, cross, 4, m, diagonal_pattern(m),
                std::vector<int>{0, 1, 2, 3}, theta, uhat, H, joint, gj);
        return result.objective + result.gradient.sum();
      },
      reps);

  const auto result = quadra::laplace::evaluate_laplace_objective_and_gradient(
      objective, cross, 4, m, diagonal_pattern(m), std::vector<int>{0, 1, 2, 3},
      theta, uhat, H, joint, gj);

  out.objective = result.objective;
  out.grad_norm = result.gradient.norm();

  return out;
}

} // namespace

int main(int argc, char **argv) {
  int reps = 20;
  if (argc > 1) {
    reps = std::stoi(argv[1]);
  }

  const std::vector<int> dims = {10, 25, 50, 100, 250};

  std::cout << "LaplaceEvaluationResult API benchmark\n";
  std::cout << "reps per case = " << reps << "\n\n";

  std::cout << std::setw(8) << "m" << std::setw(18) << "eval ms"
            << std::setw(18) << "objective" << std::setw(18) << "grad norm"
            << "\n";

  std::cout << std::scientific << std::setprecision(6);

  for (int m : dims) {
    const Result r = run_case(m, reps);

    std::cout << std::setw(8) << m << std::setw(18) << r.eval_ms
              << std::setw(18) << r.objective << std::setw(18) << r.grad_norm
              << "\n";
  }

  std::cout << "\nBenchmark complete.\n";
  return 0;
}
