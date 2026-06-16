#pragma once

#include "../model/pollock_model.hpp"

#include "../../../../../core/laplace/laplace_structure_report.hpp"
#include "../../../../../core/optimizer.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

#ifdef WALLEYE_POLLOCK_HUU_DIAGNOSTICS

namespace pollock_example {

void pollock_write_huu_diagnostics(const std::string &path, PollockModel &model,
                                   quadra::ParameterVector &params,
                                   const quadra::OptResult &fit) {
  std::ofstream out(path);
  out << std::setprecision(15);
  out << "field,value\n";
  out << "random_effects," << fit.u_hat.size() << "\n";

  if (fit.u_hat.empty()) {
    out << "available,no\n";
    out << "reason,no random effects\n";
    return;
  }

  try {
    const auto fixed_idx = quadra::build_fixed_index(params);
    const auto random_idx = quadra::build_random_index(params);

    for (std::size_t k = 0; k < fixed_idx.size() && k < fit.par.size(); ++k) {
      params.params[static_cast<std::size_t>(fixed_idx[k])].value = fit.par[k];
    }

    for (std::size_t k = 0; k < random_idx.size() && k < fit.u_hat.size();
         ++k) {
      params.params[static_cast<std::size_t>(random_idx[k])].value =
          fit.u_hat[k];
    }

    had::ADGraph graph;
    quadra::ADScope scope(graph);

    std::vector<quadra::AD> p_full;
    p_full.reserve(static_cast<std::size_t>(params.size()));
    for (int i = 0; i < params.size(); ++i) {
      p_full.emplace_back(
          quadra::AD(params.params[static_cast<std::size_t>(i)].value));
    }

    quadra::AD nll = model(p_full);
    scope.backward(nll);

    const auto &pattern = quadra::get_pattern(scope, p_full, random_idx);
    Eigen::SparseMatrix<double> H =
        quadra::extract_sparse_hessian(scope, p_full, random_idx, pattern);

    Eigen::MatrixXd dense = Eigen::MatrixXd(H);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(dense);

    out << "available,yes\n";
    out << "pattern_entries," << pattern.size() << "\n";
    out << "hessian_nonzeros," << H.nonZeros() << "\n";
    out << "min_diagonal," << dense.diagonal().minCoeff() << "\n";
    out << "max_diagonal," << dense.diagonal().maxCoeff() << "\n";

    if (es.info() == Eigen::Success) {
      const auto evals = es.eigenvalues();
      out << "eigen_success,yes\n";
      out << "min_eigenvalue," << evals.minCoeff() << "\n";
      out << "max_eigenvalue," << evals.maxCoeff() << "\n";
      out << "positive_definite," << (evals.minCoeff() > 0.0 ? "yes" : "no")
          << "\n";
      if (std::abs(evals.minCoeff()) > 0.0) {
        out << "condition_number_abs,"
            << std::abs(evals.maxCoeff()) / std::abs(evals.minCoeff()) << "\n";
      } else {
        out << "condition_number_abs,inf\n";
      }
    } else {
      out << "eigen_success,no\n";
    }
  } catch (const std::exception &e) {
    out << "available,no\n";
    out << "reason," << e.what() << "\n";
  }
}

} // namespace pollock_example

using pollock_example::pollock_write_huu_diagnostics;

#endif // WALLEYE_POLLOCK_HUU_DIAGNOSTICS
