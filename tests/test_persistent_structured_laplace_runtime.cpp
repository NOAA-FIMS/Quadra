#include "../core/laplace/hessian_structure.hpp"
#include "../core/laplace/persistent_structured_runtime.hpp"

#include <Eigen/Sparse>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using quadra::laplace::LaplaceBackendKind;
using quadra::laplace::LogDetSparseLDLT;
using quadra::laplace::PersistentStructuredLaplaceRuntime;
using quadra::laplace::StructureDetectorOptions;

void expect_close(const double a, const double b, const char *msg) {
  const double diff = std::abs(a - b);
  const double scale = 1.0 + std::max(std::abs(a), std::abs(b));

  if (diff > 1e-9 * scale) {
    std::cerr << msg << ": a=" << a << " b=" << b << " diff=" << diff << "\n";
    throw std::runtime_error(msg);
  }
}

StructureDetectorOptions detector_options() {
  StructureDetectorOptions opts;
  opts.prefer_dense_for_small_matrices = false;
  opts.dense_size_cutoff = 0;
  opts.banded_width_cutoff = 64;
  opts.dense_fill_ratio = 0.75;
  return opts;
}

Eigen::SparseMatrix<double> make_tridiagonal_scaled(const int n,
                                                    const double scale) {
  std::vector<Eigen::Triplet<double>> t;

  for (int i = 0; i < n; ++i) {
    t.emplace_back(i, i, scale * (4.0 + 0.01 * i));

    if (i > 0) {
      const double e = scale * (-0.2 + 0.001 * i);
      t.emplace_back(i, i - 1, e);
      t.emplace_back(i - 1, i, e);
    }
  }

  Eigen::SparseMatrix<double> H(n, n);
  H.setFromTriplets(t.begin(), t.end());
  H.makeCompressed();
  return H;
}

Eigen::SparseMatrix<double> make_banded_scaled(const int n, const int bandwidth,
                                               const double scale) {
  std::vector<Eigen::Triplet<double>> t;

  for (int i = 0; i < n; ++i) {
    double diag = scale * (10.0 + 0.01 * i);

    for (int d = 1; d <= bandwidth; ++d) {
      if (i - d < 0)
        continue;

      const double e =
          scale * (((d % 2 == 0) ? 0.015 : -0.025) / static_cast<double>(d));

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

void test_tridiagonal_runtime_reuses_structure() {
  PersistentStructuredLaplaceRuntime runtime(detector_options());

  const auto H1 = make_tridiagonal_scaled(50, 1.0);
  const auto H2 = make_tridiagonal_scaled(50, 1.03);

  const auto r1 = runtime.evaluate(H1);
  if (!r1.detected_structure || r1.initialized_before_call) {
    throw std::runtime_error("first tridiagonal call detection flags wrong");
  }
  if (r1.recommendation.backend != LaplaceBackendKind::Tridiagonal) {
    throw std::runtime_error("expected tridiagonal backend");
  }
  expect_close(r1.logdet, LogDetSparseLDLT(H1), "first tridiagonal runtime");

  const auto r2 = runtime.evaluate(H2);
  if (r2.detected_structure || !r2.initialized_before_call) {
    throw std::runtime_error("second tridiagonal call reuse flags wrong");
  }
  if (r2.recommendation.backend != r1.recommendation.backend) {
    throw std::runtime_error("tridiagonal backend changed");
  }
  expect_close(r2.logdet, LogDetSparseLDLT(H2), "second tridiagonal runtime");
}

void test_banded_runtime_reuses_structure() {
  PersistentStructuredLaplaceRuntime runtime(detector_options());

  const int bw = 5;
  const auto H1 = make_banded_scaled(80, bw, 1.0);
  const auto H2 = make_banded_scaled(80, bw, 0.97);

  const auto r1 = runtime.evaluate(H1);
  if (!r1.detected_structure) {
    throw std::runtime_error("first banded call did not detect");
  }
  if (r1.recommendation.backend != LaplaceBackendKind::Banded) {
    throw std::runtime_error("expected banded backend");
  }
  if (r1.recommendation.bandwidth != bw) {
    throw std::runtime_error("unexpected banded width");
  }
  expect_close(r1.logdet, LogDetSparseLDLT(H1), "first banded runtime");

  const auto r2 = runtime.evaluate(H2);
  if (r2.detected_structure) {
    throw std::runtime_error("second banded call redetected");
  }
  if (r2.recommendation.backend != r1.recommendation.backend ||
      r2.recommendation.bandwidth != r1.recommendation.bandwidth) {
    throw std::runtime_error("banded recommendation changed");
  }
  expect_close(r2.logdet, LogDetSparseLDLT(H2), "second banded runtime");
}

void test_reset_forces_redetection() {
  PersistentStructuredLaplaceRuntime runtime(detector_options());
  const auto H = make_tridiagonal_scaled(30, 1.0);

  const auto r1 = runtime.evaluate(H);
  if (!r1.detected_structure) {
    throw std::runtime_error("initial call did not detect");
  }

  runtime.reset();
  if (runtime.initialized()) {
    throw std::runtime_error("runtime still initialized after reset");
  }

  const auto r2 = runtime.evaluate(H);
  if (!r2.detected_structure) {
    throw std::runtime_error("post-reset call did not redetect");
  }
}

void test_structure_expansion_forces_redetection() {
  PersistentStructuredLaplaceRuntime runtime(detector_options());
  const int n = 30;
  std::vector<Eigen::Triplet<double>> diagonal_entries;
  for (int i = 0; i < n; ++i) {
    diagonal_entries.emplace_back(i, i, 4.0 + 0.01 * i);
  }
  Eigen::SparseMatrix<double> diagonal(n, n);
  diagonal.setFromTriplets(diagonal_entries.begin(), diagonal_entries.end());

  const auto first = runtime.evaluate(diagonal);
  if (first.recommendation.backend != LaplaceBackendKind::Diagonal) {
    throw std::runtime_error("expected initial diagonal backend");
  }

  const auto tridiagonal = make_tridiagonal_scaled(n, 1.0);
  const auto second = runtime.evaluate(tridiagonal);
  if (!second.detected_structure ||
      second.recommendation.backend != LaplaceBackendKind::Tridiagonal) {
    throw std::runtime_error(
        "diagonal-to-tridiagonal expansion was not redetected");
  }
  expect_close(second.logdet, LogDetSparseLDLT(tridiagonal),
               "expanded tridiagonal runtime");

  const auto third = runtime.evaluate(make_tridiagonal_scaled(n, 1.03));
  if (third.detected_structure) {
    throw std::runtime_error("stable expanded structure was redetected");
  }
}

int main() {
  test_tridiagonal_runtime_reuses_structure();
  test_banded_runtime_reuses_structure();
  test_reset_forces_redetection();
  test_structure_expansion_forces_redetection();

  std::cout << "persistent structured Laplace runtime tests passed\n";
  return 0;
}
