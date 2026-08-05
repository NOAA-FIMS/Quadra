#include "core/laplace/persistent_structured_runtime.hpp"

#include <Eigen/Sparse>

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ql = quadra::laplace;
using Clock = std::chrono::steady_clock;

namespace {

Eigen::SparseMatrix<double> make_precision(const std::string &kind, const int n,
                                           const double scale) {
  std::vector<Eigen::Triplet<double>> offdiag;
  std::vector<double> row_sum(static_cast<std::size_t>(n), 0.0);

  auto add_edge = [&](const int i, const int j, const double value) {
    if (i == j)
      return;
    const double scaled = scale * value;
    offdiag.emplace_back(i, j, scaled);
    offdiag.emplace_back(j, i, scaled);
    row_sum[static_cast<std::size_t>(i)] += std::abs(scaled);
    row_sum[static_cast<std::size_t>(j)] += std::abs(scaled);
  };

  int bandwidth = 0;
  if (kind == "tridiagonal")
    bandwidth = 1;
  else if (kind.rfind("banded", 0) == 0)
    bandwidth = std::stoi(kind.substr(6));
  else if (kind != "irregular")
    throw std::invalid_argument("unknown structure: " + kind);

  if (bandwidth > 0) {
    for (int i = 0; i < n; ++i) {
      for (int distance = 1; distance <= bandwidth && i + distance < n;
           ++distance) {
        add_edge(i, i + distance, -0.08 / static_cast<double>(distance));
      }
    }
  } else {
    // Markov backbone plus deterministic long-range chords. The number of
    // edges remains O(n), but the large, irregular bandwidth forces the sparse
    // backend and provides a controlled departure from banded structure.
    for (int i = 0; i + 1 < n; ++i)
      add_edge(i, i + 1, -0.08);
    for (int i = 0; i < n; i += 7) {
      const int j = (i * 37 + 17) % n;
      add_edge(i, j, -0.025);
    }
  }

  std::vector<Eigen::Triplet<double>> triplets = offdiag;
  triplets.reserve(offdiag.size() + static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    triplets.emplace_back(i, i,
                          scale * 2.0 + row_sum[static_cast<std::size_t>(i)]);
  }

  Eigen::SparseMatrix<double> matrix(n, n);
  matrix.setFromTriplets(triplets.begin(), triplets.end());
  matrix.makeCompressed();
  return matrix;
}

double reference_logdet(const Eigen::SparseMatrix<double> &matrix) {
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> factor;
  factor.compute(matrix);
  if (factor.info() != Eigen::Success)
    throw std::runtime_error("reference sparse LDLT failed");
  double value = 0.0;
  for (int i = 0; i < factor.vectorD().size(); ++i)
    value += std::log(factor.vectorD()[i]);
  return value;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "usage: quadra_structure_transition CASE N REPS\n";
    return 1;
  }

  const std::string kind = argv[1];
  const int n = std::stoi(argv[2]);
  const int reps = std::stoi(argv[3]);

  ql::StructureDetectorOptions options;
  options.prefer_dense_for_small_matrices = false;
  options.dense_size_cutoff = 0;
  options.banded_width_cutoff = 64;
  options.dense_fill_ratio = 0.75;

  const Eigen::SparseMatrix<double> initial = make_precision(kind, n, 1.0);
  ql::PersistentStructuredLaplaceRuntime runtime(options);
  const auto first = runtime.evaluate(initial);
  Eigen::SparseMatrix<double> current = initial;

  volatile double accumulator = 0.0;
  const auto start = Clock::now();
  for (int rep = 0; rep < reps; ++rep) {
    const double scale = 1.0 + 1.0e-6 * static_cast<double>(rep + 1);
    for (int slot = 0; slot < current.nonZeros(); ++slot)
      current.valuePtr()[slot] = initial.valuePtr()[slot] * scale;
    accumulator += runtime.evaluate(current).logdet;
  }
  const auto stop = Clock::now();
  (void)accumulator;

  const double elapsed_ms =
      std::chrono::duration<double, std::milli>(stop - start).count() / reps;
  const double difference = first.logdet - reference_logdet(initial);
  const auto &rec = first.recommendation;
  const auto &state = runtime.state();

  std::cout << std::setprecision(17) << kind << "," << n << ","
            << initial.nonZeros() << "," << rec.bandwidth << ","
            << ql::ToString(rec.structure) << "," << ql::ToString(rec.backend)
            << "," << state.sparse_ldlt_workspace.symbolic_analysis_count()
            << "," << state.sparse_ldlt_workspace.numeric_factorization_count()
            << "," << elapsed_ms << "," << difference << "\n";
  return std::abs(difference) <= 1.0e-9 ? 0 : 2;
}
