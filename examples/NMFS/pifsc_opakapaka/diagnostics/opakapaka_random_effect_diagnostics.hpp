#pragma once

#include "../quadra/opakapaka_model.hpp"

#include "../../../../core/uncertainty/reporting.hpp"
#include "../../../../core/uncertainty/selected_inverse_diagonal.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

namespace opakapaka_example {

inline void
write_random_effect_uncertainty_csv(const std::string &path,
                                    const std::vector<double> &u_hat,
                                    const Eigen::SparseMatrix<double> &h_uu) {
  const auto diag =
      quadra::uncertainty::selected_inverse_diagonal_from_spd_hessian(h_uu);

  std::ofstream out(path);
  out << "effect,mode,conditional_se,conditional_variance,note\n";

  for (std::size_t i = 0; i < u_hat.size(); ++i) {
    double se = std::numeric_limits<double>::quiet_NaN();
    double var = std::numeric_limits<double>::quiet_NaN();
    std::string note = diag.message;

    if (diag.success && i < diag.standard_error.size() &&
        i < diag.variance.size()) {
      se = diag.standard_error[i];
      var = diag.variance[i];
      note = "selected_inverse_diagonal";
    }

    out << "log_B[" << i << "]," << u_hat[i] << "," << se << "," << var << ","
        << note << "\n";
  }
}

template <class Model>
Eigen::SparseMatrix<double> compute_final_random_effect_hessian(
    Model &model, quadra::ParameterVector &params,
    quadra::LaplaceOptions & /*opts*/, const quadra::OptResult &fit) {
  // QUADRA_OPAKAPAKA_HUU_ADSCOPE_REPAIR_V1
  //
  // LaplaceResult currently stores value/gradients only. For conditional
  // random-effect SEs, rebuild the fitted AD vector, evaluate the model,
  // propagate adjoints, discover the sparse Hessian pattern, and extract H_uu
  // using Quadra's sparse Hessian extraction API.

  const std::size_t n_fixed = fit.par.size();
  const std::size_t n_random = fit.u_hat.size();
  const std::size_t n_total = n_fixed + n_random;

  std::vector<int> random_idx;
  random_idx.reserve(n_random);
  for (std::size_t i = 0; i < n_random; ++i) {
    random_idx.push_back(static_cast<int>(n_fixed + i));
  }

  // QUADRA_OPAKAPAKA_HUU_CURRENT_API_REPAIR_V1
  had::ADGraph graph;
  quadra::ADScope scope(graph);

  std::vector<quadra::AD> p_full;
  p_full.reserve(n_total);

  for (std::size_t i = 0; i < n_fixed; ++i) {
    p_full.emplace_back(quadra::AD(fit.par.at(i)));
  }
  for (std::size_t i = 0; i < n_random; ++i) {
    p_full.emplace_back(quadra::AD(fit.u_hat.at(i)));
  }

  quadra::AD nll = model(p_full);
  scope.backward(nll);

  const auto &pattern = quadra::get_pattern(scope, p_full, random_idx);
  auto h_uu =
      quadra::extract_sparse_hessian(scope, p_full, random_idx, pattern);

  h_uu.makeCompressed();
  return h_uu;
}

} // namespace opakapaka_example

// Compatibility aliases for the current Opakapaka driver.
using opakapaka_example::compute_final_random_effect_hessian;
using opakapaka_example::write_random_effect_uncertainty_csv;
