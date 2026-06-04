#include "../core/laplace/structured_value_factory.hpp"
#include "../core/laplace/hessian_structure.hpp"

#include <Eigen/Sparse>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using quadra::laplace::BackendRecommendation;
using quadra::laplace::LaplaceBackendKind;
using quadra::laplace::LogDetSparseLDLT;
using quadra::laplace::StructureDetector;
using quadra::laplace::StructureDetectorOptions;
using quadra::laplace::StructuredValues;
using quadra::laplace::extract_structured_values;
using quadra::laplace::logdet_structured_values;
using quadra::laplace::update_structured_values_from_hessian;

void expect_close(const double a, const double b, const char* msg) {
  const double diff = std::abs(a - b);
  const double scale = 1.0 + std::max(std::abs(a), std::abs(b));

  if (diff > 1e-9 * scale) {
    std::cerr << msg << ": a=" << a << " b=" << b
              << " diff=" << diff << "\\n";
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

Eigen::SparseMatrix<double> make_diag(const int n, const double scale) {
  std::vector<Eigen::Triplet<double>> t;
  for (int i = 0; i < n; ++i) {
    t.emplace_back(i, i, scale * (2.0 + 0.01 * i));
  }
  Eigen::SparseMatrix<double> H(n, n);
  H.setFromTriplets(t.begin(), t.end());
  H.makeCompressed();
  return H;
}

Eigen::SparseMatrix<double> make_tri(const int n, const double scale) {
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

Eigen::SparseMatrix<double> make_banded(const int n,
                                        const int bandwidth,
                                        const double scale) {
  std::vector<Eigen::Triplet<double>> t;
  for (int i = 0; i < n; ++i) {
    double diag = scale * (10.0 + 0.01 * i);
    for (int d = 1; d <= bandwidth; ++d) {
      if (i - d < 0) continue;
      const double e =
          scale * (((d % 2 == 0) ? 0.015 : -0.025) /
                   static_cast<double>(d));
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

void test_case(const Eigen::SparseMatrix<double>& H1,
               const Eigen::SparseMatrix<double>& H2,
               const LaplaceBackendKind expected,
               const char* label) {
  StructureDetector detector(detector_options());
  const BackendRecommendation rec = detector.Analyze(H1);

  if (rec.backend != expected) {
    throw std::runtime_error(std::string(label) + " unexpected backend");
  }

  StructuredValues values = extract_structured_values(H1, rec);
  expect_close(logdet_structured_values(values), LogDetSparseLDLT(H1),
               (std::string(label) + " initial").c_str());

  update_structured_values_from_hessian(values, H2, rec);
  expect_close(logdet_structured_values(values), LogDetSparseLDLT(H2),
               (std::string(label) + " updated").c_str());
}

void test_reject_unsupported() {
  BackendRecommendation rec;
  rec.backend = LaplaceBackendKind::SparseLDLT;
  StructuredValues values;
  const auto H = make_tri(10, 1.0);

  bool threw = false;
  try {
    update_structured_values_from_hessian(values, H, rec);
  } catch (...) {
    threw = true;
  }

  if (!threw) {
    throw std::runtime_error("unsupported update was not rejected");
  }
}

int main() {
  test_case(make_diag(20, 1.0), make_diag(20, 1.05),
            LaplaceBackendKind::Diagonal, "diagonal native update");

  test_case(make_tri(50, 1.0), make_tri(50, 0.97),
            LaplaceBackendKind::Tridiagonal, "tridiagonal native update");

  test_case(make_banded(80, 5, 1.0), make_banded(80, 5, 1.02),
            LaplaceBackendKind::Banded, "banded native update");

  test_reject_unsupported();

  std::cout << "native structured value update tests passed\\n";
  return 0;
}
