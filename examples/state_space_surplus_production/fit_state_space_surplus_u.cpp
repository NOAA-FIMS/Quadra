#include "state_space_surplus_production.hpp"

#include <Eigen/Dense>
#include <LBFGS.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

namespace ss = quadra_examples::state_space_surplus_production;

namespace {

std::vector<double> to_std_vector(const Eigen::VectorXd &x) {
  std::vector<double> out(static_cast<std::size_t>(x.size()));
  for (int i = 0; i < x.size(); ++i) {
    out[static_cast<std::size_t>(i)] = x[i];
  }
  return out;
}

Eigen::VectorXd to_eigen_vector(const std::vector<double> &x) {
  Eigen::VectorXd out(static_cast<int>(x.size()));
  for (std::size_t i = 0; i < x.size(); ++i) {
    out[static_cast<int>(i)] = x[i];
  }
  return out;
}

double objective_u(const ss::Data &data, const ss::Parameters &par,
                   const Eigen::VectorXd &u) {
  try {
    const std::vector<double> u_std = to_std_vector(u);
    return ss::joint_objective(data, par, u_std);
  } catch (...) {
    return std::numeric_limits<double>::infinity();
  }
}

Eigen::VectorXd finite_difference_gradient_u(const ss::Data &data,
                                             const ss::Parameters &par,
                                             const Eigen::VectorXd &u) {
  Eigen::VectorXd grad(u.size());

  for (int i = 0; i < u.size(); ++i) {
    const double step = 1e-5 * (1.0 + std::abs(u[i]));

    Eigen::VectorXd plus = u;
    Eigen::VectorXd minus = u;
    plus[i] += step;
    minus[i] -= step;

    const double f_plus = objective_u(data, par, plus);
    const double f_minus = objective_u(data, par, minus);

    if (!std::isfinite(f_plus) || !std::isfinite(f_minus)) {
      grad[i] = 0.0;
    } else {
      grad[i] = (f_plus - f_minus) / (2.0 * step);
    }
  }

  return grad;
}

class RandomEffectsObjective {
public:
  RandomEffectsObjective(const ss::Data &data, const ss::Parameters &par)
      : data_(data), par_(par) {}

  double operator()(const Eigen::VectorXd &u, Eigen::VectorXd &grad) {
    const double f = objective_u(data_, par_, u);
    grad = finite_difference_gradient_u(data_, par_, u);
    return f;
  }

private:
  const ss::Data &data_;
  const ss::Parameters &par_;
};

struct FitResult {
  Eigen::VectorXd u_hat;
  double objective = 0.0;
  double grad_norm = 0.0;
  int iterations = 0;
  bool converged = false;
  bool accepted_line_search_failure = false;
};

FitResult fit_random_effects(const ss::Data &data, const ss::Parameters &par) {
  const std::vector<double> u0_std = ss::zero_random_effects(data);
  Eigen::VectorXd u = to_eigen_vector(u0_std);

  LBFGSpp::LBFGSParam<double> param;
  param.epsilon = 1e-6;
  param.max_iterations = 500;
  param.max_linesearch = 100;
  param.m = 8;
  param.ftol = 1e-4;
  param.wolfe = 0.9;
  param.min_step = 1e-20;
  param.max_step = 1.0;

  LBFGSpp::LBFGSSolver<double> solver(param);
  RandomEffectsObjective objective(data, par);

  FitResult result;
  result.u_hat = u;

  try {
    result.iterations = solver.minimize(objective, u, result.objective);
    result.converged = true;
  } catch (const std::exception &e) {
    result.objective = objective_u(data, par, u);
    result.grad_norm = finite_difference_gradient_u(data, par, u).norm();

    if (std::isfinite(result.objective) && result.grad_norm < 1e-3) {
      result.converged = true;
      result.accepted_line_search_failure = true;
      std::cout << "LBFGS++ terminated during line search after reaching "
                << "an acceptable random-effects optimum.\n";
      std::cout << "  message    = " << e.what() << "\n";
      std::cout << "  objective  = " << result.objective << "\n";
      std::cout << "  grad_norm  = " << result.grad_norm << "\n\n";
    } else {
      std::cerr << "LBFGS++ random-effects optimization failed: " << e.what()
                << "\n";
      std::cerr << "Current objective = " << result.objective << "\n";
      std::cerr << "Current grad_norm = " << result.grad_norm << "\n";
      throw;
    }
  }

  result.u_hat = u;
  result.objective = objective_u(data, par, u);
  result.grad_norm = finite_difference_gradient_u(data, par, u).norm();

  return result;
}

void print_u_summary(const Eigen::VectorXd &u) {
  double mean = 0.0;
  double min_value = std::numeric_limits<double>::infinity();
  double max_value = -std::numeric_limits<double>::infinity();
  double ssq = 0.0;

  for (int i = 0; i < u.size(); ++i) {
    mean += u[i];
    ssq += u[i] * u[i];
    min_value = std::min(min_value, u[i]);
    max_value = std::max(max_value, u[i]);
  }

  mean /= static_cast<double>(u.size());

  double sd = 0.0;
  for (int i = 0; i < u.size(); ++i) {
    sd += (u[i] - mean) * (u[i] - mean);
  }
  const int sd_denominator = std::max(1, static_cast<int>(u.size()) - 1);
  sd = std::sqrt(sd / static_cast<double>(sd_denominator));

  std::cout << "Random effects summary\n";
  std::cout << "  n          = " << u.size() << "\n";
  std::cout << "  mean       = " << mean << "\n";
  std::cout << "  sd         = " << sd << "\n";
  std::cout << "  min        = " << min_value << "\n";
  std::cout << "  max        = " << max_value << "\n";
  std::cout << "  norm       = " << std::sqrt(ssq) << "\n\n";

  std::cout << std::setw(8) << "t" << std::setw(16) << "u_hat" << "\n";

  for (int i = 0; i < u.size(); ++i) {
    std::cout << std::setw(8) << i << std::setw(16) << u[i] << "\n";
  }

  std::cout << "\n";
}

} // namespace

int main() {
  const ss::Data data = ss::make_demo_data();
  const ss::Parameters par = ss::make_demo_parameters();

  const std::vector<double> u0 = ss::zero_random_effects(data);
  const double initial_objective = ss::joint_objective(data, par, u0);

  std::cout << std::fixed << std::setprecision(6);

  std::cout << "Fit state-space surplus production random effects\n";
  std::cout << "=================================================\n\n";
  std::cout << "Fixed effects are held constant.\n";
  std::cout << "initial joint objective at u = 0: " << initial_objective
            << "\n\n";

  const FitResult fit = fit_random_effects(data, par);
  const std::vector<double> u_hat = to_std_vector(fit.u_hat);

  std::cout << "Random-effects fit summary\n";
  std::cout << "  converged  = " << (fit.converged ? "yes" : "no") << "\n";
  std::cout << "  accepted line-search termination = "
            << (fit.accepted_line_search_failure ? "yes" : "no") << "\n";
  std::cout << "  iterations = " << fit.iterations << "\n";
  std::cout << "  objective  = " << fit.objective << "\n";
  std::cout << "  grad_norm  = " << fit.grad_norm << "\n";
  std::cout << "  improvement = " << (initial_objective - fit.objective)
            << "\n\n";

  print_u_summary(fit.u_hat);

  std::cout << "State-space report at u_hat\n";
  std::cout << "===========================\n";
  ss::print_report(data, par, u_hat);

  return 0;
}
