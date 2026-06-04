#include "../core/laplace/model_analysis_report.hpp"

#include <Eigen/Sparse>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using quadra::laplace::LaplaceBackendKind;
using quadra::laplace::ModelAnalysisReport;
using quadra::laplace::SolverRecommendation;
using quadra::laplace::StructureDetectorOptions;
using quadra::laplace::analyze_hessian_structure;

StructureDetectorOptions detector_options() {
  StructureDetectorOptions opts;
  opts.prefer_dense_for_small_matrices = false;
  opts.dense_size_cutoff = 0;
  opts.banded_width_cutoff = 64;
  opts.dense_fill_ratio = 0.75;
  return opts;
}

Eigen::SparseMatrix<double> make_diagonal(const int n) {
  std::vector<Eigen::Triplet<double>> t;
  for (int i = 0; i < n; ++i) {
    t.emplace_back(i, i, 2.0 + 0.01 * i);
  }

  Eigen::SparseMatrix<double> H(n, n);
  H.setFromTriplets(t.begin(), t.end());
  H.makeCompressed();
  return H;
}

Eigen::SparseMatrix<double> make_tridiagonal(const int n) {
  std::vector<Eigen::Triplet<double>> t;
  for (int i = 0; i < n; ++i) {
    t.emplace_back(i, i, 4.0 + 0.01 * i);

    if (i > 0) {
      const double e = -0.2 + 0.001 * i;
      t.emplace_back(i, i - 1, e);
      t.emplace_back(i - 1, i, e);
    }
  }

  Eigen::SparseMatrix<double> H(n, n);
  H.setFromTriplets(t.begin(), t.end());
  H.makeCompressed();
  return H;
}

Eigen::SparseMatrix<double> make_banded(const int n, const int bandwidth) {
  std::vector<Eigen::Triplet<double>> t;
  for (int i = 0; i < n; ++i) {
    double diag = 10.0 + 0.01 * i;

    for (int d = 1; d <= bandwidth; ++d) {
      if (i - d < 0) continue;

      const double e =
          ((d % 2 == 0) ? 0.015 : -0.025) / static_cast<double>(d);

      t.emplace_back(i, i - d, e);
      t.emplace_back(i - d, i, e);
      diag += 2.0 * std::abs(e);
    }

    t.emplace_back(i, i, diag);
  }

  Eigen::SparseMatrix<double> H(n, n);
  H.setFromTriplets(t.begin(), t.end());
  H.makeCompressed();
  return H;
}

void test_diagonal_report() {
  const ModelAnalysisReport report =
      analyze_hessian_structure(make_diagonal(20), detector_options());

  if (!report.is_diagonal()) {
    throw std::runtime_error("diagonal report did not detect diagonal");
  }

  if (report.backend != LaplaceBackendKind::Diagonal) {
    throw std::runtime_error("diagonal report did not recommend diagonal backend");
  }

  if (report.solver != SolverRecommendation::Newton) {
    throw std::runtime_error("diagonal report did not recommend Newton");
  }

  const std::string text = report.summary();
  if (text.find("Model Analysis Report") == std::string::npos) {
    throw std::runtime_error("summary missing title");
  }
}

void test_tridiagonal_report() {
  const ModelAnalysisReport report =
      analyze_hessian_structure(make_tridiagonal(50), detector_options());

  if (!report.is_tridiagonal()) {
    throw std::runtime_error("tridiagonal report did not detect tridiagonal");
  }

  if (report.backend != LaplaceBackendKind::Tridiagonal) {
    throw std::runtime_error("tridiagonal report did not recommend tridiagonal backend");
  }

  if (report.bandwidth != 1) {
    throw std::runtime_error("tridiagonal report wrong bandwidth");
  }
}

void test_banded_report() {
  const ModelAnalysisReport report =
      analyze_hessian_structure(make_banded(80, 5), detector_options());

  if (!report.is_banded()) {
    throw std::runtime_error("banded report did not detect banded");
  }

  if (report.backend != LaplaceBackendKind::Banded) {
    throw std::runtime_error("banded report did not recommend banded backend");
  }

  if (report.bandwidth != 5) {
    throw std::runtime_error("banded report wrong bandwidth");
  }
}

int main() {
  test_diagonal_report();
  test_tridiagonal_report();
  test_banded_report();

  std::cout << "model analysis report tests passed\n";
  return 0;
}
