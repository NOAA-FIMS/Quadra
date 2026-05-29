#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../core/laplace/laplace_evaluator_exact_gradient_integration.hpp"
#include "../core/laplace/sparse_laplace_evaluation_result.hpp"

DECLARE_ADGRAPH()

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

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

Eigen::MatrixXd Huu_dense(const Eigen::VectorXd &theta,
                          const Eigen::VectorXd &u) {
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
    Eigen::LDLT<Eigen::MatrixXd> ldlt(Huu_dense(theta, u));
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

void run_case(int m) {
  Eigen::VectorXd theta(5);
  theta << 0.55, std::log(0.65), std::log(0.55), std::log(0.90), std::log(0.25);

  const Eigen::VectorXd uhat = solve_uhat(theta, m);

  const Eigen::MatrixXd dense = Huu_dense(theta, uhat);
  const Eigen::SparseMatrix<double> direct = Huu_sparse_direct(theta, uhat);
  const Eigen::MatrixXd direct_dense = Eigen::MatrixXd(direct);

  const double H_diff = (dense - direct_dense).cwiseAbs().maxCoeff();

  if (H_diff > 1.0e-14) {
    std::cerr << "H_diff = " << H_diff << "\n";
    throw std::runtime_error("direct sparse Huu does not match dense Huu.");
  }

  const double joint = joint_objective_double(theta, uhat);
  const Eigen::VectorXd gj = joint_envelope_gradient(theta, uhat);

  SparseRw1Objective objective{m};
  auto cross = [theta, uhat](int theta_index) -> Eigen::VectorXd {
    return f_u_theta_column(theta, uhat, theta_index);
  };

  quadra::laplace::LaplaceExactGradientEvaluationInputs dense_inputs;
  dense_inputs.theta = theta;
  dense_inputs.uhat = uhat;
  dense_inputs.Huu = dense;
  dense_inputs.joint_objective = joint;
  dense_inputs.joint_envelope_gradient = gj;

  const auto dense_result =
      quadra::laplace::evaluate_laplace_with_exact_gradient(
          objective, cross, 5, m, tridiagonal_pattern(m),
          std::vector<int>{0, 1, 2, 3, 4}, dense_inputs);

  quadra::laplace::SparseLaplaceExactGradientEvaluationInputs sparse_inputs;
  sparse_inputs.theta = theta;
  sparse_inputs.uhat = uhat;
  sparse_inputs.Huu = direct;
  sparse_inputs.joint_objective = joint;
  sparse_inputs.joint_envelope_gradient = gj;

  const auto sparse_result =
      quadra::laplace::evaluate_sparse_laplace_with_exact_gradient(
          objective, cross, 5, m, tridiagonal_pattern(m),
          std::vector<int>{0, 1, 2, 3, 4}, sparse_inputs);

  const double obj_diff =
      std::abs(dense_result.objective - sparse_result.objective);
  const double grad_diff =
      (dense_result.gradient - sparse_result.gradient).cwiseAbs().maxCoeff();

  std::cout << "\n=== direct sparse RW1 Huu validation: m = " << m << " ===\n";
  std::cout << std::scientific << std::setprecision(10);
  std::cout << "H diff         = " << H_diff << "\n";
  std::cout << "objective diff = " << obj_diff << "\n";
  std::cout << "gradient diff  = " << grad_diff << "\n";

  if (obj_diff > 1.0e-9) {
    throw std::runtime_error(
        "direct sparse objective differs from dense objective.");
  }

  if (grad_diff > 1.0e-8) {
    throw std::runtime_error(
        "direct sparse gradient differs from dense gradient.");
  }
}

} // namespace

int main() {
  run_case(10);
  run_case(25);
  run_case(50);

  std::cout << "\ndirect sparse RW1 Huu builder tests passed\n";
  return 0;
}
