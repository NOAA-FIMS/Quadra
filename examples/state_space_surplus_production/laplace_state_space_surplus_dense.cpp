#include "state_space_surplus_production.hpp"

#include <Eigen/Dense>
#include <LBFGS.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace ss = quadra_examples::state_space_surplus_production;

namespace {

std::vector<double> to_std_vector(const Eigen::VectorXd& x) {
  std::vector<double> out(static_cast<std::size_t>(x.size()));
  for (int i = 0; i < x.size(); ++i) {
    out[static_cast<std::size_t>(i)] = x[i];
  }
  return out;
}

Eigen::VectorXd to_eigen_vector(const std::vector<double>& x) {
  Eigen::VectorXd out(static_cast<int>(x.size()));
  for (std::size_t i = 0; i < x.size(); ++i) {
    out[static_cast<int>(i)] = x[i];
  }
  return out;
}

double joint_u(const ss::Data& data,
               const ss::Parameters& par,
               const Eigen::VectorXd& u) {
  return ss::joint_objective(data, par, to_std_vector(u));
}

Eigen::VectorXd finite_difference_gradient_u(const ss::Data& data,
                                             const ss::Parameters& par,
                                             const Eigen::VectorXd& u) {
  Eigen::VectorXd grad(u.size());

  for (int i = 0; i < u.size(); ++i) {
    const double step = 1e-5 * (1.0 + std::abs(u[i]));

    Eigen::VectorXd plus = u;
    Eigen::VectorXd minus = u;
    plus[i] += step;
    minus[i] -= step;

    const double f_plus = joint_u(data, par, plus);
    const double f_minus = joint_u(data, par, minus);

    grad[i] = (f_plus - f_minus) / (2.0 * step);
  }

  return grad;
}

class RandomEffectsObjective {
 public:
  RandomEffectsObjective(const ss::Data& data, const ss::Parameters& par)
      : data_(data), par_(par) {}

  double operator()(const Eigen::VectorXd& u, Eigen::VectorXd& grad) {
    const double f = joint_u(data_, par_, u);
    grad = finite_difference_gradient_u(data_, par_, u);
    return f;
  }

 private:
  const ss::Data& data_;
  const ss::Parameters& par_;
};

struct UHatResult {
  Eigen::VectorXd u_hat;
  double joint = 0.0;
  double grad_norm = 0.0;
  int iterations = 0;
  bool converged = false;
  bool accepted_line_search_failure = false;
};

UHatResult optimize_u_hat(const ss::Data& data, const ss::Parameters& par) {
  Eigen::VectorXd u = to_eigen_vector(ss::zero_random_effects(data));

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
  RandomEffectsObjective objective(data, par);

  UHatResult out;
  out.u_hat = u;

  try {
    out.iterations = solver.minimize(objective, u, out.joint);
    out.converged = true;
  } catch (const std::exception& e) {
    out.joint = joint_u(data, par, u);
    out.grad_norm = finite_difference_gradient_u(data, par, u).norm();

    if (std::isfinite(out.joint) && out.grad_norm < 1e-3) {
      out.converged = true;
      out.accepted_line_search_failure = true;
      std::cout << "u_hat solve accepted line-search termination: "
                << e.what() << "\n";
    } else {
      throw;
    }
  }

  out.u_hat = u;
  out.joint = joint_u(data, par, u);
  out.grad_norm = finite_difference_gradient_u(data, par, u).norm();

  return out;
}

Eigen::MatrixXd finite_difference_hessian_uu(const ss::Data& data,
                                             const ss::Parameters& par,
                                             const Eigen::VectorXd& u) {
  const int n = static_cast<int>(u.size());
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(n, n);

  const double f0 = joint_u(data, par, u);

  for (int i = 0; i < n; ++i) {
    const double hi = 1e-4 * (1.0 + std::abs(u[i]));

    Eigen::VectorXd up = u;
    Eigen::VectorXd um = u;
    up[i] += hi;
    um[i] -= hi;

    const double fp = joint_u(data, par, up);
    const double fm = joint_u(data, par, um);

    H(i, i) = (fp - 2.0 * f0 + fm) / (hi * hi);

    for (int j = i + 1; j < n; ++j) {
      const double hj = 1e-4 * (1.0 + std::abs(u[j]));

      Eigen::VectorXd upp = u;
      Eigen::VectorXd upm = u;
      Eigen::VectorXd ump = u;
      Eigen::VectorXd umm = u;

      upp[i] += hi;
      upp[j] += hj;

      upm[i] += hi;
      upm[j] -= hj;

      ump[i] -= hi;
      ump[j] += hj;

      umm[i] -= hi;
      umm[j] -= hj;

      const double fpp = joint_u(data, par, upp);
      const double fpm = joint_u(data, par, upm);
      const double fmp = joint_u(data, par, ump);
      const double fmm = joint_u(data, par, umm);

      const double hij = (fpp - fpm - fmp + fmm) / (4.0 * hi * hj);

      H(i, j) = hij;
      H(j, i) = hij;
    }
  }

  return H;
}

struct LaplaceResult {
  UHatResult uhat;
  Eigen::MatrixXd Huu;
  double logdet_Huu = 0.0;
  double laplace_correction = 0.0;
  double objective = 0.0;
  double min_eigenvalue = 0.0;
  bool hessian_spd = false;
};

LaplaceResult evaluate_dense_laplace(const ss::Data& data,
                                     const ss::Parameters& par) {
  LaplaceResult out;

  out.uhat = optimize_u_hat(data, par);
  out.Huu = finite_difference_hessian_uu(data, par, out.uhat.u_hat);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(out.Huu);
  out.min_eigenvalue = es.eigenvalues().minCoeff();

  Eigen::LLT<Eigen::MatrixXd> llt(out.Huu);
  out.hessian_spd = (llt.info() == Eigen::Success);

  if (!out.hessian_spd) {
    throw std::runtime_error("Huu is not SPD at u_hat");
  }

  const auto& L = llt.matrixL();

  out.logdet_Huu = 0.0;
  for (int i = 0; i < out.Huu.rows(); ++i) {
    out.logdet_Huu += 2.0 * std::log(L(i, i));
  }

  const double n_u = static_cast<double>(out.uhat.u_hat.size());
  out.laplace_correction =
      0.5 * out.logdet_Huu - 0.5 * n_u * std::log(2.0 * M_PI);

  out.objective = out.uhat.joint + out.laplace_correction;

  return out;
}

void print_u_hat(const Eigen::VectorXd& u) {
  std::cout << std::setw(8) << "t"
            << std::setw(16) << "u_hat"
            << "\n";

  for (int i = 0; i < u.size(); ++i) {
    std::cout << std::setw(8) << i
              << std::setw(16) << u[i]
              << "\n";
  }

  std::cout << "\n";
}

}  // namespace

int main() {
  const ss::Data data = ss::make_demo_data();
  const ss::Parameters par = ss::make_demo_parameters();

  std::cout << std::fixed << std::setprecision(6);

  std::cout << "Dense Laplace state-space surplus production example\n";
  std::cout << "====================================================\n\n";

  const double joint_at_zero =
      ss::joint_objective(data, par, ss::zero_random_effects(data));

  std::cout << "joint(theta, u=0) = " << joint_at_zero << "\n\n";

  const LaplaceResult result = evaluate_dense_laplace(data, par);

  std::cout << "u_hat solve\n";
  std::cout << "  converged  = " << (result.uhat.converged ? "yes" : "no") << "\n";
  std::cout << "  accepted line-search termination = "
            << (result.uhat.accepted_line_search_failure ? "yes" : "no") << "\n";
  std::cout << "  iterations = " << result.uhat.iterations << "\n";
  std::cout << "  joint(theta, u_hat) = " << result.uhat.joint << "\n";
  std::cout << "  grad_norm = " << result.uhat.grad_norm << "\n\n";

  std::cout << "Huu diagnostics\n";
  std::cout << "  n_u            = " << result.Huu.rows() << "\n";
  std::cout << "  min eigenvalue = " << result.min_eigenvalue << "\n";
  std::cout << "  SPD            = " << (result.hessian_spd ? "yes" : "no") << "\n";
  std::cout << "  logdet(Huu)    = " << result.logdet_Huu << "\n\n";

  std::cout << "Laplace objective\n";
  std::cout << "  joint          = " << result.uhat.joint << "\n";
  std::cout << "  correction     = " << result.laplace_correction << "\n";
  std::cout << "  marginal nll   = " << result.objective << "\n\n";

  print_u_hat(result.uhat.u_hat);

  std::cout << "State-space report at u_hat\n";
  std::cout << "===========================\n";
  ss::print_report(data, par, to_std_vector(result.uhat.u_hat));

  return 0;
}
