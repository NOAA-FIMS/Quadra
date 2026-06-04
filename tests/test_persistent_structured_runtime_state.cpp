#include "../core/laplace/persistent_structured_runtime.hpp"
#include "../core/laplace/hessian_structure.hpp"

#include <Eigen/Sparse>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using quadra::laplace::LaplaceBackendKind;
using quadra::laplace::LogDetSparseLDLT;
using quadra::laplace::PersistentStructuredRuntimeState;
using quadra::laplace::StructureDetectorOptions;

void expect_close(const double a, const double b, const char* msg) {
  const double diff = std::abs(a - b);
  const double scale = 1.0 + std::max(std::abs(a), std::abs(b));

  if (diff > 1e-9 * scale) {
    std::cerr << msg << ": a=" << a << " b=" << b
              << " diff=" << diff << "\n";
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

Eigen::SparseMatrix<double> make_banded_scaled(const int n,
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

void test_uninitialized_rejected() {
  PersistentStructuredRuntimeState state;

  bool threw = false;
  try {
    (void)state.logdet();
  } catch (...) {
    threw = true;
  }

  if (!threw) {
    throw std::runtime_error("uninitialized logdet was not rejected");
  }
}

void test_case(const Eigen::SparseMatrix<double>& H,
               const LaplaceBackendKind expected_backend,
               const char* label) {
  PersistentStructuredRuntimeState state;
  state.update_from_hessian(H, detector_options());

  if (!state.initialized) {
    throw std::runtime_error("state was not initialized");
  }

  if (state.backend_recommendation().backend != expected_backend) {
    throw std::runtime_error(std::string(label) + " backend mismatch");
  }

  expect_close(state.logdet(), LogDetSparseLDLT(H), label);
}

void test_update_values_only_tridiagonal() {
  PersistentStructuredRuntimeState state;

  const Eigen::SparseMatrix<double> H1 = make_tridiagonal_scaled(50, 1.0);
  const Eigen::SparseMatrix<double> H2 = make_tridiagonal_scaled(50, 1.05);

  state.update_from_hessian(H1, detector_options());
  const auto initial_backend = state.backend_recommendation().backend;
  const int initial_bandwidth = state.backend_recommendation().bandwidth;

  state.update_values_only(H2);

  if (state.backend_recommendation().backend != initial_backend) {
    throw std::runtime_error("update_values_only changed backend recommendation");
  }
  if (state.backend_recommendation().bandwidth != initial_bandwidth) {
    throw std::runtime_error("update_values_only changed bandwidth");
  }
  expect_close(state.logdet(), LogDetSparseLDLT(H2),
               "update_values_only tridiagonal");
}

void test_update_values_only_banded() {
  PersistentStructuredRuntimeState state;

  const Eigen::SparseMatrix<double> H1 = make_banded_scaled(80, 5, 1.0);
  const Eigen::SparseMatrix<double> H2 = make_banded_scaled(80, 5, 0.95);

  state.update_from_hessian(H1, detector_options());
  const auto initial_backend = state.backend_recommendation().backend;
  const int initial_bandwidth = state.backend_recommendation().bandwidth;

  state.update_values_only(H2);

  if (state.backend_recommendation().backend != initial_backend) {
    throw std::runtime_error("update_values_only changed backend recommendation");
  }
  if (state.backend_recommendation().bandwidth != initial_bandwidth) {
    throw std::runtime_error("update_values_only changed bandwidth");
  }
  expect_close(state.logdet(), LogDetSparseLDLT(H2),
               "update_values_only banded");
}

void test_update_values_only_rejects_uninitialized() {
  PersistentStructuredRuntimeState state;
  const Eigen::SparseMatrix<double> H = make_tridiagonal(10);
  bool threw = false;
  try {
    state.update_values_only(H);
  } catch (...) {
    threw = true;
  }
  if (!threw) {
    throw std::runtime_error("update_values_only did not reject uninitialized state");
  }
}

void test_clear() {
  PersistentStructuredRuntimeState state;
  state.update_from_hessian(make_tridiagonal(10), detector_options());

  if (!state.initialized) {
    throw std::runtime_error("state not initialized before clear");
  }

  state.clear();

  if (state.initialized) {
    throw std::runtime_error("state still initialized after clear");
  }
}

int main() {
  test_uninitialized_rejected();
  test_update_values_only_rejects_uninitialized();
  test_update_values_only_tridiagonal();
  test_update_values_only_banded();

  test_case(make_diagonal(20), LaplaceBackendKind::Diagonal,
            "persistent diagonal");
  test_case(make_tridiagonal(50), LaplaceBackendKind::Tridiagonal,
            "persistent tridiagonal");
  test_case(make_banded(80, 5), LaplaceBackendKind::Banded,
            "persistent banded");

  test_clear();

  std::cout << "persistent structured runtime state tests passed\n";
  return 0;
}
