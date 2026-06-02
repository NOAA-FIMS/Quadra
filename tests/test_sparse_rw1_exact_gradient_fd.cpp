#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../core/laplace/laplace_evaluator_exact_gradient_integration.hpp"

DECLARE_ADGRAPH()

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

// theta = [mu, log_sigma, log_lambda0, log_lambda_rw, log_beta]
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

Eigen::MatrixXd Huu(const Eigen::VectorXd &theta, const Eigen::VectorXd &u) {
  const int m = static_cast<int>(u.size());
  const double inv_sigma2 = std::exp(-2.0 * theta[1]);
  const double lambda0 = std::exp(theta[2]);
  const double lambda_rw = std::exp(theta[3]);
  const double beta = std::exp(theta[4]);

  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(m, m);

  for (int i = 0; i < m; ++i) {
    H(i, i) += inv_sigma2 + lambda0 + beta * std::exp(u[i]);
  }

  for (int i = 1; i < m; ++i) {
    H(i, i) += lambda_rw;
    H(i - 1, i - 1) += lambda_rw;
    H(i, i - 1) -= lambda_rw;
    H(i - 1, i) -= lambda_rw;
  }

  return H;
}

Eigen::VectorXd solve_uhat(const Eigen::VectorXd &theta, int m) {
  Eigen::VectorXd u = Eigen::VectorXd::Zero(m);

  for (int iter = 0; iter < 80; ++iter) {
    const Eigen::VectorXd g = random_gradient(theta, u);
    const Eigen::MatrixXd H = Huu(theta, u);

    Eigen::LDLT<Eigen::MatrixXd> ldlt(H);
    if (ldlt.info() != Eigen::Success) {
      throw std::runtime_error("LDLT failed in solve_uhat.");
    }

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

  throw std::out_of_range("theta_index out of range.");
}

double laplace_objective(const Eigen::VectorXd &theta, int m) {
  const Eigen::VectorXd uhat = solve_uhat(theta, m);
  const Eigen::MatrixXd H = Huu(theta, uhat);

  Eigen::LDLT<Eigen::MatrixXd> ldlt(H);
  if (ldlt.info() != Eigen::Success) {
    throw std::runtime_error("LDLT failed in laplace_objective.");
  }

  return joint_objective_double(theta, uhat) +
         0.5 * ldlt.vectorD().array().log().sum();
}

Eigen::VectorXd central_fd_gradient(const Eigen::VectorXd &theta, int m,
                                    double step = 1.0e-5) {
  Eigen::VectorXd g(theta.size());

  for (int j = 0; j < theta.size(); ++j) {
    const double h = step * std::max(1.0, std::abs(theta[j]));
    Eigen::VectorXd plus = theta;
    Eigen::VectorXd minus = theta;

    plus[j] += h;
    minus[j] -= h;

    g[j] =
        (laplace_objective(plus, m) - laplace_objective(minus, m)) / (2.0 * h);
  }

  return g;
}

void run_case(int m, const Eigen::VectorXd &theta) {
  const Eigen::VectorXd uhat = solve_uhat(theta, m);
  const Eigen::MatrixXd H = Huu(theta, uhat);
  const double joint = joint_objective_double(theta, uhat);
  const Eigen::VectorXd gj = joint_envelope_gradient(theta, uhat);

  SparseRw1Objective objective{m};

  auto cross = [theta, uhat](int theta_index) -> Eigen::VectorXd {
    return f_u_theta_column(theta, uhat, theta_index);
  };

  quadra::laplace::LaplaceExactGradientEvaluationInputs inputs;
  inputs.theta = theta;
  inputs.uhat = uhat;
  inputs.Huu = H;
  inputs.joint_objective = joint;
  inputs.joint_envelope_gradient = gj;

  const auto result = quadra::laplace::evaluate_laplace_with_exact_gradient(
      objective, cross, 5, m, tridiagonal_pattern(m),
      std::vector<int>{0, 1, 2, 3, 4}, inputs);

  const Eigen::VectorXd fd = central_fd_gradient(theta, m);
  const double max_abs = (result.gradient - fd).cwiseAbs().maxCoeff();

  std::cout << "\n=== sparse RW1 exact-gradient FD validation: m = " << m
            << " ===\n";
  std::cout << std::scientific << std::setprecision(10);
  std::cout << "exact = " << result.gradient.transpose() << "\n";
  std::cout << "fd    = " << fd.transpose() << "\n";
  std::cout << "diff  = " << (result.gradient - fd).cwiseAbs().transpose()
            << "\n";
  std::cout << "max abs diff = " << max_abs << "\n";

  if (max_abs > 2.0e-6) {
    throw std::runtime_error("sparse RW1 exact gradient failed FD validation.");
  }
}

} // namespace

int main() {
  Eigen::VectorXd theta(5);

  theta << 0.55, std::log(0.65), std::log(0.55), std::log(0.90), std::log(0.25);
  run_case(8, theta);
  run_case(20, theta);

  theta << 0.45, std::log(0.80), std::log(0.70), std::log(1.20), std::log(0.35);
  run_case(15, theta);

  std::cout << "\nsparse RW1 exact-gradient FD validation passed\n";
  return 0;
}
