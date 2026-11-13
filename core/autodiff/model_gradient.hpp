#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

#include "../autodiff.hpp"
#include "../model/model_context.hpp"
#include "../model/model_result.hpp"

namespace quadra {

/// First-order objective, coordinate gradient, norm and reports; no Hessian is
/// constructed.
struct GradientResult {
  double objective_value_m = 0.0;
  double gradient_norm_m = 0.0;
  std::vector<double> parameters_m;
  std::vector<double> gradient_m;
  std::vector<ReportValue> reports_m;
};

/// Return the Euclidean norm of an ordinary coordinate vector.
inline double euclidean_norm(const std::vector<double> &x) {
  double out = 0.0;
  for (double xi : x) {
    out += xi * xi;
  }
  return std::sqrt(out);
}

/// Evaluate a model on a fresh AD tape and reverse first-order adjoints only.
/// @param model Initialized through its Quadra interface on a fresh evaluation
/// context.
/// @param parameters Independent model coordinates in parameter_names order.
/// @return Objective value, coordinate gradient/norm, input coordinates and
/// report values.
/// @throws std::invalid_argument If parameters is empty; model evaluation
/// errors propagate.
/// @details Reversing only first-order adjoints is intentional. Building a
/// Hessian here would add memory/work that GradientResult cannot return. This
/// method does not change the model target or silently integrate any random
/// effects.
template <class Model>
GradientResult evaluate_gradient(Model &model,
                                 const std::vector<double> &parameters) {
  if (parameters.empty()) {
    throw std::invalid_argument(
        "evaluate_gradient: parameter vector cannot be empty");
  }

  TapeContext tape;
  ADScope scope(tape.graph);

  ModelReportContext ctx;
  model.initialize(ctx);

  std::vector<AD> x = to_ad(parameters);
  AD objective = model.template evaluate<AD>(x, ctx);

  // GradientResult contains no curvature; avoid building the full Hessian.
  scope.backward_first_order(objective);
  Eigen::VectorXd g_eigen = extract_gradient(x);

  GradientResult result;
  result.objective_value_m = value_of(objective);
  result.parameters_m = parameters;
  result.reports_m = ctx.reports().values();
  result.gradient_m.resize(static_cast<size_t>(g_eigen.size()));

  for (int i = 0; i < g_eigen.size(); ++i) {
    result.gradient_m[static_cast<size_t>(i)] = g_eigen[i];
  }

  result.gradient_norm_m = euclidean_norm(result.gradient_m);
  return result;
}

} // namespace quadra
