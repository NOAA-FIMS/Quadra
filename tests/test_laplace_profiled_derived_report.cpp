#include "../core/laplace/laplace_profiled_derived_report.hpp"
#include "../core/model/parameter_transform.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <vector>

DECLARE_ADGRAPH()

int main() {
  std::vector<double> theta{2.0, -0.5};
  std::vector<double> u{1.5, -2.0};

  Eigen::MatrixXd theta_covariance(2, 2);
  theta_covariance << 0.04, 0.01, 0.01, 0.09;

  quadra::LaplaceImplicitWorkspace workspace;
  workspace.success_m = true;
  workspace.message_m = "test workspace";

  workspace.du_dtheta_m = Eigen::MatrixXd(2, 2);
  workspace.du_dtheta_m << 0.5, -0.1, 0.2, 0.3;

  std::vector<quadra::ProfiledDerivedQuantity> quantities;

  quantities.push_back({"q1", [](const std::vector<quadra::AD> &theta,
                                 const std::vector<quadra::AD> &u) {
                          return theta[0] * theta[0] + 2.0 * u[0] - 3.0 * u[1];
                        }});

  quantities.push_back({"q2", [](const std::vector<quadra::AD> &theta,
                                 const std::vector<quadra::AD> &u) {
                          return exp(theta[1]) + theta[0] * u[1];
                        }});

  const auto report = quadra::compute_laplace_profiled_derived_report(
      quantities, theta, u, theta_covariance, workspace);

  if (!report.success_m) {
    std::cerr << "FAIL: profiled report failed: " << report.message_m << "\n";
    return 1;
  }

  if (report.names_m.size() != 2 || report.delta_m.estimate_m.size() != 2 ||
      report.delta_m.jacobian_m.rows() != 2 ||
      report.delta_m.jacobian_m.cols() != 2) {
    std::cerr << "FAIL: unexpected report dimensions\n";
    return 1;
  }

  if (!report.delta_m.covariance_m.allFinite() ||
      !report.delta_m.correlation_m.allFinite()) {
    std::cerr << "FAIL: non-finite report covariance/correlation\n";
    return 1;
  }

  std::cout << "PASS: Laplace profiled derived report\n";
  std::cout << "  names: " << report.names_m[0] << ", " << report.names_m[1]
            << "\n";
  std::cout << "  estimates: " << report.delta_m.estimate_m.transpose() << "\n";
  std::cout << "  profiled Jacobian:\n" << report.delta_m.jacobian_m << "\n";
  std::cout << "  covariance:\n" << report.delta_m.covariance_m << "\n";
  std::cout << "  correlation:\n" << report.delta_m.correlation_m << "\n";

  return 0;
}
