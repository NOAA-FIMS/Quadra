/**
 * @file summarize_curvature.cpp
 * @details
 * Invert validated local curvature without ridge regularization. Input rows
 * identify perturbed coordinates and columns identify gradient components, so
 * load places them in Hessian-column order. Symmetrization follows
 * asymmetry/step checks; singular or indefinite curvature is rejected.
 */
#include "swordfish/csv.hpp"
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
/// Split a generated Hessian/coordinate record into comma-delimited fields for
/// curvature summarization.
std::vector<std::string> split(const std::string &s) {
  return swordfish::csv_fields(s);
}
/// Load a square generated Hessian by labeled perturbation rows and gradient
/// columns, checking row/header correspondence.
Eigen::MatrixXd load(const std::string &path) {
  std::ifstream in(path);
  std::string line;
  std::getline(in, line);
  const auto names = split(line);
  if (names.size() < 2 || names[0] != "parameter")
    throw std::runtime_error("Hessian header");
  int n = names.size() - 1;
  Eigen::MatrixXd a(n, n);
  for (int j = 0; j < n; ++j) {
    if (!std::getline(in, line))
      throw std::runtime_error("Incomplete Hessian");
    const auto cells = split(line);
    if (cells.size() != names.size() || cells[0] != names[j + 1])
      throw std::runtime_error("Hessian labels mismatch");
    for (int i = 0; i < n; ++i)
      a(i, j) = swordfish::csv_number(cells[i + 1]);
  }
  if (std::getline(in, line) || !a.allFinite())
    throw std::runtime_error("Trailing/nonfinite Hessian values");
  return a;
}
/// Run summarize curvature; return zero on success and nonzero on a failed
/// check or invalid input.
int main() {
  try {
    auto a = load("build/hessian_coarse.csv"),
         b = load("build/hessian_fine.csv");
    if (a.rows() != b.rows())
      throw std::runtime_error("Hessian dimensions differ");
    const double asym = (b - b.transpose()).norm() / b.norm(),
                 difference = (a - b).norm() / b.norm();
    const Eigen::MatrixXd h = (b + b.transpose()) * .5;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(h);
    if (eig.info() != Eigen::Success || eig.eigenvalues().minCoeff() <= 0 ||
        asym > 1e-5 || difference > 1e-5)
      throw std::runtime_error(
          "Curvature not sufficiently stable and positive definite");
    const Eigen::MatrixXd covariance =
        eig.eigenvectors() * eig.eigenvalues().cwiseInverse().asDiagonal() *
        eig.eigenvectors().transpose();
    const Eigen::MatrixXd hc = (a + a.transpose()) * .5;
    const Eigen::MatrixXd cc =
        hc.ldlt().solve(Eigen::MatrixXd::Identity(h.rows(), h.cols()));
    double se_change = 0;
    for (int i = 0; i < h.rows(); ++i)
      se_change = std::max(
          se_change, std::abs(std::sqrt(cc(i, i) / covariance(i, i)) - 1));
    if (se_change > 1e-3)
      throw std::runtime_error("Parameter SEs sensitive to Hessian step");
    std::ofstream summary("build/curvature_summary.json");
    summary << std::setprecision(16)
            << "{\"method\":\"central differences of exact AD gradients; logit "
               "coordinates; joint "
               "likelihood\",\"steps\":[0.0001,0.00005],\"relative_asymmetry\":"
            << asym << ",\"relative_step_difference\":" << difference
            << ",\"minimum_eigenvalue\":" << eig.eigenvalues().minCoeff()
            << ",\"maximum_eigenvalue\":" << eig.eigenvalues().maxCoeff()
            << ",\"condition_number\":"
            << eig.eigenvalues().maxCoeff() / eig.eigenvalues().minCoeff()
            << ",\"maximum_relative_SE_step_change\":" << se_change << "}\n";
    std::ifstream in("build/uncertainty_coordinates.csv");
    std::string line;
    std::getline(in, line);
    std::ofstream out("build/parameter_uncertainty.csv");
    out << "parameter,estimate,standard_error,lower95,upper95,coordinate_"
           "standard_error\n"
        << std::setprecision(16);
    std::vector<std::string> names;
    int i = 0;
    while (std::getline(in, line)) {
      auto v = split(line);
      if (v.size() != 5 || i >= h.rows())
        throw std::runtime_error("Coordinate dimensions");
      names.push_back(v[0]);
      double x = swordfish::csv_number(v[1]),
             value = swordfish::csv_number(v[2]),
             lo = swordfish::csv_number(v[3]), hi = swordfish::csv_number(v[4]);
      const double s = std::sqrt(covariance(i, i)),
                   q = (value - lo) / (hi - lo);
      auto transform = [&](double z) {
        return lo + (hi - lo) / (1 + std::exp(-z));
      };
      out << v[0] << ',' << value << ',' << (hi - lo) * q * (1 - q) * s << ','
          << transform(x - 1.95996398454005 * s) << ','
          << transform(x + 1.95996398454005 * s) << ',' << s << '\n';
      ++i;
    }
    if (i != h.rows())
      throw std::runtime_error("Incomplete coordinates");
    std::ofstream cov("build/coordinate_covariance.csv");
    cov << "parameter";
    for (const auto &name : names)
      cov << ',' << name;
    cov << '\n' << std::setprecision(17);
    for (int j = 0; j < i; ++j) {
      cov << names[j];
      for (int k = 0; k < i; ++k)
        cov << ',' << covariance(j, k);
      cov << '\n';
    }
    std::ifstream jac("build/derived_jacobian.csv");
    if (std::getline(jac, line)) {
      const auto headers = split(line);
      const int m = headers.size() - 1;
      Eigen::MatrixXd j(m, i);
      for (int k = 0; k < i; ++k) {
        if (!std::getline(jac, line))
          throw std::runtime_error("Incomplete derived Jacobian");
        const auto row = split(line);
        if (row.size() != headers.size() || row[0] != names[k])
          throw std::runtime_error("Derived Jacobian order mismatch");
        for (int q = 0; q < m; ++q)
          j(q, k) = swordfish::csv_number(row[q + 1]);
      }
      const Eigen::MatrixXd v = j * covariance * j.transpose();
      std::ifstream estimates("build/derived_estimates.csv");
      std::getline(estimates, line);
      std::ofstream derived("build/derived_uncertainty.csv");
      derived << "quantity,estimate,standard_error,lower95,upper95\n"
              << std::setprecision(16);
      for (int q = 0; q < m; ++q) {
        std::getline(estimates, line);
        auto row = split(line);
        if (row.size() != 2 || row[0] != headers[q + 1])
          throw std::runtime_error("Derived estimate order mismatch");
        double value = swordfish::csv_number(row[1]), se = std::sqrt(v(q, q));
        derived << row[0] << ',' << value << ',' << se << ','
                << value - 1.95996398454005 * se << ','
                << value + 1.95996398454005 * se << '\n';
      }
    }
    if (!summary || !out || !cov)
      throw std::runtime_error("Curvature output failure");
    std::cout << "Positive definite Hessian; max relative SE step change "
              << se_change << '\n';
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
