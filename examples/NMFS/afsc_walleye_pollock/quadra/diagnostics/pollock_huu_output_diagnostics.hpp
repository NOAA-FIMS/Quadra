#pragma once

#include "../model/pollock_laplace_helpers.hpp"
#include "../model/pollock_model.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

namespace pollock_example {

void pollock_write_huu_matrix(const std::string &path, PollockModel &model,
                              quadra::ParameterVector params,
                              const quadra::OptResult &fit) {
  std::ofstream out(path);
  out << std::setprecision(15);

  const auto random_idx = quadra::build_random_index(params);
  const std::size_t n = random_idx.size();

  out << "row";
  for (std::size_t j = 0; j < n; ++j) {
    out << ",u" << (j + 1);
  }
  out << "\n";

  Eigen::MatrixXd dense = pollock_fd_huu(model, params, fit);

  for (std::size_t i = 0; i < n; ++i) {
    out << "u" << (i + 1);
    for (std::size_t j = 0; j < n; ++j) {
      out << ","
          << dense(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
    }
    out << "\n";
  }
}

void pollock_write_huu_pattern_compare(const std::string &path,
                                       PollockModel &model,
                                       quadra::ParameterVector params,
                                       const quadra::OptResult &fit,
                                       double tol = 1.0e-8) {
  std::ofstream out(path);
  out << "field,value\n";

  const auto fixed_idx = quadra::build_fixed_index(params);
  const auto random_idx = quadra::build_random_index(params);
  const std::size_t n = random_idx.size();

  out << "random_effects," << n << "\n";
  out << "fd_tol," << tol << "\n";
  out << "quadra_pattern_available," << (fit.pattern.available ? "yes" : "no")
      << "\n";
  out << "quadra_pattern_detected_structure," << fit.pattern.detected_structure
      << "\n";
  out << "quadra_pattern_nonzeros_reported," << fit.pattern.nonzeros << "\n";

  if (n == 0 || fit.par.size() != fixed_idx.size() || fit.u_hat.size() != n) {
    out << "available,no\n";
    out << "reason,missing random effects or size mismatch\n";
    return;
  }

  Eigen::MatrixXd Hfd = pollock_fd_huu(model, params, fit);

  std::size_t fd_nonzeros_all = 0;
  std::size_t fd_nonzeros_upper = 0;
  std::size_t fd_nonzeros_diag = 0;
  double max_abs_fd = 0.0;
  double min_abs_fd_nonzero = std::numeric_limits<double>::infinity();

  for (Eigen::Index i = 0; i < Hfd.rows(); ++i) {
    for (Eigen::Index j = 0; j < Hfd.cols(); ++j) {
      const double av = std::abs(Hfd(i, j));
      max_abs_fd = std::max(max_abs_fd, av);
      if (av > tol) {
        ++fd_nonzeros_all;
        min_abs_fd_nonzero = std::min(min_abs_fd_nonzero, av);
        if (i <= j) {
          ++fd_nonzeros_upper;
        }
        if (i == j) {
          ++fd_nonzeros_diag;
        }
      }
    }
  }

  const std::size_t fd_nonzeros_offdiag_all =
      fd_nonzeros_all >= fd_nonzeros_diag ? fd_nonzeros_all - fd_nonzeros_diag
                                          : 0;
  const std::size_t fd_nonzeros_offdiag_upper =
      fd_nonzeros_upper >= fd_nonzeros_diag
          ? fd_nonzeros_upper - fd_nonzeros_diag
          : 0;

  out << "available,yes\n";
  out << "fd_nonzeros_all," << fd_nonzeros_all << "\n";
  out << "fd_nonzeros_upper_including_diag," << fd_nonzeros_upper << "\n";
  out << "fd_nonzeros_diag," << fd_nonzeros_diag << "\n";
  out << "fd_nonzeros_offdiag_all," << fd_nonzeros_offdiag_all << "\n";
  out << "fd_nonzeros_offdiag_upper," << fd_nonzeros_offdiag_upper << "\n";
  out << "fd_density_all,"
      << (n == 0 ? 0.0
                 : static_cast<double>(fd_nonzeros_all) /
                       static_cast<double>(n * n))
      << "\n";
  out << "fd_density_upper,"
      << (n == 0 ? 0.0
                 : static_cast<double>(fd_nonzeros_upper) /
                       static_cast<double>((n * (n + 1)) / 2))
      << "\n";
  out << "max_abs_fd," << max_abs_fd << "\n";
  out << "min_abs_fd_nonzero,"
      << (std::isfinite(min_abs_fd_nonzero) ? min_abs_fd_nonzero : 0.0) << "\n";
  out << "note,OptPatternInfo does not currently expose individual pattern "
         "entries; this compares reported Quadra count to finite-difference "
         "numerical sparsity.\n";

  std::ofstream detail("examples/NMFS/afsc_walleye_pollock/outputs/"
                       "walleye_pollock_huu_pattern_compare_detail.csv");
  detail << "i,j,fd_nonzero,fd_value,abs_fd_value,band_distance\n";
  detail << std::setprecision(15);

  for (Eigen::Index i = 0; i < Hfd.rows(); ++i) {
    for (Eigen::Index j = 0; j < Hfd.cols(); ++j) {
      const double v = Hfd(i, j);
      const double av = std::abs(v);
      if (av > tol) {
        detail << (i + 1) << "," << (j + 1) << ",yes," << v << "," << av << ","
               << std::abs(i - j) << "\n";
      }
    }
  }
}

} // namespace pollock_example

// Compatibility aliases for current walleye_pollock.cpp call sites.
using pollock_example::pollock_write_huu_matrix;
using pollock_example::pollock_write_huu_pattern_compare;
