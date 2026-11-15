/** @file optimizer_scaling.hpp
 * Constant linear optimizer scaling; statistical target and original gradients
 * remain unchanged. A covariance used here is an optimization metric, not an
 * uncertainty estimate for the new fit.
 */
#pragma once
#include "core/model/quadra_model.hpp"
#include "csv.hpp"
#include "delta_method.hpp"
#include <fstream>

namespace swordfish {
/// Read a finite positive-definite optimizer metric with exact coordinate
/// labels.
inline Eigen::MatrixXd
read_optimizer_metric(const std::string &path,
                      const std::vector<std::string> &names) {
  std::ifstream in(path);
  std::string line;
  if (!csv_line(in, line))
    throw std::invalid_argument("Missing optimizer metric");
  const auto header = csv_fields(line);
  if (header.size() != names.size() + 1 || header[0] != "quantity")
    throw std::invalid_argument("Optimizer metric header");
  for (size_t i = 0; i < names.size(); ++i)
    if (header[i + 1] != names[i])
      throw std::invalid_argument("Optimizer metric coordinate order");
  Eigen::MatrixXd matrix(names.size(), names.size());
  for (size_t i = 0; i < names.size(); ++i) {
    if (!csv_line(in, line))
      throw std::invalid_argument("Incomplete optimizer metric");
    const auto row = csv_fields(line);
    if (row.size() != header.size() || row[0] != names[i])
      throw std::invalid_argument("Optimizer metric row order");
    for (size_t j = 0; j < names.size(); ++j)
      matrix(i, j) = csv_number(row[j + 1]);
  }
  if (csv_line(in, line))
    throw std::invalid_argument("Trailing optimizer metric rows");
  validate_parameter_covariance(matrix);
  return matrix;
}

/// Optimize z while evaluating the unchanged objective at x=center+factor*z.
/// This constant MLE coordinate change adds no statistical density Jacobian.
template <class Model>
struct OptimizerScaling : quadra::QuadraModel<OptimizerScaling<Model>> {
  Model &model;
  std::vector<double> center;
  Eigen::MatrixXd factor;
  /// Borrow the objective and copy its starting coordinates and linear factor.
  OptimizerScaling(Model &m, std::vector<double> x, Eigen::MatrixXd scaling)
      : model(m), center(std::move(x)), factor(std::move(scaling)) {
    if (factor.rows() != int(center.size()) ||
        factor.cols() != int(center.size()) || !factor.allFinite())
      throw std::invalid_argument("Optimizer scaling dimensions");
  }
  /// Preserve canonical coordinate order through the optimizer adapter.
  std::vector<std::string> parameter_names_impl() const {
    return model.parameter_names_impl();
  }
  /// Expand scaled optimizer coordinates into original model coordinates.
  template <class T> std::vector<T> expand(const std::vector<T> &z) const {
    if (z.size() != center.size())
      throw std::invalid_argument("Scaled coordinate dimension");
    std::vector<T> x;
    for (size_t i = 0; i < z.size(); ++i) {
      T value = T(center[i]);
      for (size_t j = 0; j < z.size(); ++j)
        value += T(factor(i, j)) * z[j];
      x.push_back(value);
    }
    return x;
  }
  /// Evaluate the same scalar likelihood after constant optimizer scaling.
  template <class T> T operator()(const std::vector<T> &z) const {
    return model(expand(z));
  }
  /// Expose the scaled objective through the regular Quadra model interface.
  template <class T>
  T evaluate_impl(const std::vector<T> &z, quadra::ModelReportContext &) const {
    return (*this)(z);
  }
};
} // namespace swordfish
