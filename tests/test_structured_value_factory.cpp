#include "../core/laplace/hessian_structure.hpp"
#include "../core/laplace/structured_value_factory.hpp"

#include <Eigen/Sparse>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using quadra::laplace::BackendRecommendation;
using quadra::laplace::extract_structured_values;
using quadra::laplace::HessianStructure;
using quadra::laplace::LaplaceBackendKind;
using quadra::laplace::logdet_structured_values;
using quadra::laplace::LogDetSparseLDLT;
using quadra::laplace::StructureDetector;
using quadra::laplace::StructureDetectorOptions;
using quadra::laplace::StructureInfo;
using quadra::laplace::ToString;

void expect_close(const double a, const double b, const char *msg) {
  const double diff = std::abs(a - b);
  const double scale = 1.0 + std::max(std::abs(a), std::abs(b));
  if (diff > 1e-9 * scale) {
    std::cerr << msg << ": a=" << a << " b=" << b << " diff=" << diff << "\n";
    throw std::runtime_error(msg);
  }
}

const char *structure_name(const HessianStructure s) {
  switch (s) {
  case HessianStructure::Diagonal:
    return "diagonal";
  case HessianStructure::Tridiagonal:
    return "tridiagonal";
  case HessianStructure::Banded:
    return "banded";
  case HessianStructure::SparsePattern:
    return "sparse_pattern";
  case HessianStructure::Dense:
    return "dense";
  }
  return "unknown";
}

Eigen::SparseMatrix<double> make_diagonal() {
  std::vector<Eigen::Triplet<double>> t;
  t.emplace_back(0, 0, 2.0);
  t.emplace_back(1, 1, 3.0);
  t.emplace_back(2, 2, 4.0);
  t.emplace_back(3, 3, 5.0);
  Eigen::SparseMatrix<double> H(4, 4);
  H.setFromTriplets(t.begin(), t.end());
  H.makeCompressed();
  return H;
}

Eigen::SparseMatrix<double> make_tridiagonal() {
  const int n = 6;
  std::vector<Eigen::Triplet<double>> t;
  for (int i = 0; i < n; ++i) {
    t.emplace_back(i, i, 4.0 + 0.1 * i);
    if (i > 0) {
      const double e = -0.2 + 0.01 * i;
      t.emplace_back(i, i - 1, e);
      t.emplace_back(i - 1, i, e);
    }
  }
  Eigen::SparseMatrix<double> H(n, n);
  H.setFromTriplets(t.begin(), t.end());
  H.makeCompressed();
  return H;
}

Eigen::SparseMatrix<double> make_banded() {
  const int n = 8;
  std::vector<Eigen::Triplet<double>> t;
  for (int i = 0; i < n; ++i) {
    t.emplace_back(i, i, 5.0 + 0.1 * i);
    if (i > 0) {
      const double e1 = -0.25 + 0.01 * i;
      t.emplace_back(i, i - 1, e1);
      t.emplace_back(i - 1, i, e1);
    }
    if (i > 1) {
      const double e2 = 0.04 + 0.005 * i;
      t.emplace_back(i, i - 2, e2);
      t.emplace_back(i - 2, i, e2);
    }
  }
  Eigen::SparseMatrix<double> H(n, n);
  H.setFromTriplets(t.begin(), t.end());
  H.makeCompressed();
  return H;
}

StructureDetector make_detector() {
  StructureDetectorOptions opts;
  opts.prefer_dense_for_small_matrices = true;
  opts.dense_size_cutoff = 16;
  opts.banded_width_cutoff = 64;
  opts.dense_fill_ratio = 0.75;
  return StructureDetector(opts);
}

void analyze_case(const std::string &label,
                  const Eigen::SparseMatrix<double> &H,
                  const LaplaceBackendKind expected_backend) {
  StructureDetector detector = make_detector();
  const BackendRecommendation rec = detector.Analyze(H);
  const StructureInfo info = quadra::laplace::InspectHessianStructure(
      H, detector.options().structure_options);

  if (rec.backend != expected_backend) {
    std::cerr << label << " expected backend " << ToString(expected_backend)
              << " but got " << ToString(rec.backend)
              << " reason=" << rec.reason << "\n";
    throw std::runtime_error("unexpected backend recommendation");
  }

  const auto values = extract_structured_values(H, rec);
  const double structured_logdet = logdet_structured_values(values);
  const double sparse_logdet = LogDetSparseLDLT(H);
  const double diff = structured_logdet - sparse_logdet;
  expect_close(structured_logdet, sparse_logdet,
               (label + " structured factory logdet").c_str());

  std::cout << std::setw(14) << label << std::setw(8) << info.rows
            << std::setw(8) << info.nnz << std::setw(10) << info.max_bandwidth
            << std::setw(14) << info.fill_ratio << std::setw(16)
            << structure_name(info.detected) << std::setw(16)
            << ToString(rec.backend) << std::setw(24) << rec.reason
            << std::setw(18) << diff << "\n";
}

void test_unsupported_backend_rejected() {
  const Eigen::SparseMatrix<double> H = make_tridiagonal();
  BackendRecommendation rec;
  rec.backend = LaplaceBackendKind::SparseLDLT;
  bool threw = false;
  try {
    (void)extract_structured_values(H, rec);
  } catch (...) {
    threw = true;
  }
  if (!threw) {
    throw std::runtime_error("unsupported backend was not rejected");
  }
}

int main() {
  std::cout << std::fixed << std::setprecision(12);
  std::cout << "structured value factory structural analysis\n";
  std::cout << std::setw(14) << "case" << std::setw(8) << "rows" << std::setw(8)
            << "nnz" << std::setw(10) << "band" << std::setw(14) << "fill"
            << std::setw(16) << "structure" << std::setw(16) << "backend"
            << std::setw(24) << "reason" << std::setw(18) << "logdet_diff"
            << "\n";

  analyze_case("diagonal", make_diagonal(), LaplaceBackendKind::Diagonal);
  analyze_case("tridiagonal", make_tridiagonal(),
               LaplaceBackendKind::Tridiagonal);
  analyze_case("banded", make_banded(), LaplaceBackendKind::Banded);
  test_unsupported_backend_rejected();
  std::cout << "structured value factory tests passed\n";
  return 0;
}
