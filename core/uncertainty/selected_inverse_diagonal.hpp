#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>

#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace quadra {
namespace uncertainty {

struct SelectedInverseDiagonalResult {
  std::vector<double> variance;
  std::vector<double> standard_error;
  bool success = false;
  std::string message;
};

// Level-1 conditional random-effect covariance utility.
//
// For a positive-definite random-effect Hessian H_uu, the conditional
// covariance is approximately inv(H_uu). This conservative implementation
// extracts diag(inv(H_uu)) by solving H_uu x_i = e_i and reading x_i[i].
//
// This is deliberately simple and validation-friendly. Later, this can be
// replaced with a true selected-inverse algorithm for very large models.
inline SelectedInverseDiagonalResult selected_inverse_diagonal_from_spd_hessian(
    const Eigen::SparseMatrix<double>& hessian,
    double min_variance = 0.0) {
  SelectedInverseDiagonalResult out;

  const int n = static_cast<int>(hessian.rows());
  if (hessian.rows() != hessian.cols()) {
    out.message = "Hessian is not square";
    return out;
  }

  out.variance.assign(static_cast<std::size_t>(n),
                      std::numeric_limits<double>::quiet_NaN());
  out.standard_error.assign(static_cast<std::size_t>(n),
                            std::numeric_limits<double>::quiet_NaN());

  if (n == 0) {
    out.success = true;
    out.message = "empty Hessian";
    return out;
  }

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
  ldlt.analyzePattern(hessian);
  ldlt.factorize(hessian);

  if (ldlt.info() != Eigen::Success) {
    out.message = "SimplicialLDLT factorization failed";
    return out;
  }

  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n);
  for (int i = 0; i < n; ++i) {
    rhs.setZero();
    rhs[i] = 1.0;

    const Eigen::VectorXd sol = ldlt.solve(rhs);
    if (ldlt.info() != Eigen::Success) {
      out.message = "SimplicialLDLT solve failed";
      return out;
    }

    double v = sol[i];

    if (std::isfinite(v) && v < 0.0 && std::abs(v) <= 1.0e-12) {
      v = 0.0;
    }

    if (std::isfinite(v) && v < min_variance) {
      v = min_variance;
    }

    out.variance[static_cast<std::size_t>(i)] = v;
    out.standard_error[static_cast<std::size_t>(i)] =
        (std::isfinite(v) && v >= 0.0) ? std::sqrt(v)
                                       : std::numeric_limits<double>::quiet_NaN();
  }

  out.success = true;
  out.message = "ok";
  return out;
}

}  // namespace uncertainty
}  // namespace quadra
