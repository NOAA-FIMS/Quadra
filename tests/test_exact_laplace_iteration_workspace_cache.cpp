#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../core/laplace/exact_laplace_iteration_workspace.hpp"

DECLARE_ADGRAPH()

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

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

void run_case(int m) {
  Eigen::VectorXd theta(4);
  theta << 0.65, std::log(0.55), std::log(0.85), std::log(0.30);

  const Eigen::VectorXd uhat = solve_uhat(theta, m);
  const Eigen::MatrixXd H = Huu(theta, uhat);
  const Eigen::VectorXd gj = joint_envelope_gradient(theta, uhat);

  NonlinearDiagObjective objective{m};

  auto cross = [theta, uhat](int theta_index) -> Eigen::VectorXd {
    return f_u_theta_column(theta, uhat, theta_index);
  };

  auto workspace = quadra::laplace::make_exact_laplace_iteration_workspace(
      objective, cross, 4, m, diagonal_pattern(m), std::vector<int>{0, 1, 2, 3},
      theta, uhat, H);

  if (workspace.trace_terms_prepared()) {
    throw std::runtime_error("cache should not be prepared initially.");
  }

  const Eigen::VectorXd g1 = workspace.gradient(gj);

  if (!workspace.trace_terms_prepared()) {
    throw std::runtime_error("cache should be prepared after gradient.");
  }

  const Eigen::VectorXd logdet = workspace.logdet_gradient_contribution();
  const Eigen::VectorXd g2 = workspace.gradient(gj);
  const Eigen::VectorXd g3 = gj + logdet;

  if ((g1 - g2).cwiseAbs().maxCoeff() > 0.0) {
    throw std::runtime_error("cached repeated gradient changed.");
  }

  if ((g1 - g3).cwiseAbs().maxCoeff() > 1.0e-14) {
    throw std::runtime_error(
        "gradient does not equal envelope + cached logdet.");
  }

  const auto &Hdots = workspace.cached_Hdots();
  if (static_cast<int>(Hdots.size()) != 4) {
    throw std::runtime_error("wrong number of cached Hdot matrices.");
  }
}

} // namespace

int main() {
  run_case(10);
  run_case(25);
  run_case(50);

  std::cout << "exact Laplace iteration workspace cache tests passed\n";
  return 0;
}
