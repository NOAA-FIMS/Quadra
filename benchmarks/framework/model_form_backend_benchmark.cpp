#include "core/laplace/laplace_backend_factory.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace ql = quadra::laplace;
using Clock = std::chrono::steady_clock;

namespace {

struct ModelCase {
  std::string form;
  Eigen::SparseMatrix<double> hessian;
  ql::LaplaceBackendKind expected;
};

Eigen::SparseMatrix<double>
precision_from_edges(int n, const std::vector<std::pair<int, int>> &edges,
                     bool nearly_dense = false) {
  std::vector<Eigen::Triplet<double>> entries;
  std::vector<double> row_sum(static_cast<std::size_t>(n), 0.0);
  for (const auto &edge : edges) {
    const int i = edge.first;
    const int j = edge.second;
    if (i == j)
      continue;
    const double value = -0.025 * (1.0 + ((i + j) % 3));
    entries.emplace_back(i, j, value);
    entries.emplace_back(j, i, value);
    row_sum[static_cast<std::size_t>(i)] += std::abs(value);
    row_sum[static_cast<std::size_t>(j)] += std::abs(value);
  }
  for (int i = 0; i < n; ++i)
    entries.emplace_back(i, i,
                         2.0 + row_sum[static_cast<std::size_t>(i)] +
                             (nearly_dense ? 0.5 : 0.0));
  Eigen::SparseMatrix<double> out(n, n);
  out.setFromTriplets(entries.begin(), entries.end());
  out.makeCompressed();
  return out;
}

std::vector<std::pair<int, int>> band_edges(int n, int width) {
  std::vector<std::pair<int, int>> edges;
  for (int i = 0; i < n; ++i)
    for (int d = 1; d <= width && i + d < n; ++d)
      edges.emplace_back(i, i + d);
  return edges;
}

std::vector<ModelCase> catalog(int n) {
  std::vector<ModelCase> out;
  out.push_back({"diagonal", precision_from_edges(n, {}),
                 ql::LaplaceBackendKind::Diagonal});
  out.push_back({"tridiagonal", precision_from_edges(n, band_edges(n, 1)),
                 ql::LaplaceBackendKind::Tridiagonal});
  out.push_back({"banded", precision_from_edges(n, band_edges(n, 4)),
                 ql::LaplaceBackendKind::Banded});

  std::vector<std::pair<int, int>> blocks;
  constexpr int block_size = 8;
  for (int start = 0; start < n; start += block_size)
    for (int i = start; i < std::min(n, start + block_size); ++i)
      for (int j = i + 1; j < std::min(n, start + block_size); ++j)
        blocks.emplace_back(i, j);
  out.push_back({"block_diagonal", precision_from_edges(n, blocks),
                 ql::LaplaceBackendKind::Banded});

  std::vector<std::pair<int, int>> arrowhead = band_edges(n, 1);
  for (int i = 1; i < n; ++i)
    arrowhead.emplace_back(0, i);
  out.push_back({"arrowhead", precision_from_edges(n, arrowhead),
                 ql::LaplaceBackendKind::SparseLDLT});

  std::vector<std::pair<int, int>> sparse = band_edges(n, 1);
  for (int i = 0; i < n; i += 7) {
    const int j = (i * 37 + 17) % n;
    if (i != j)
      sparse.emplace_back(std::min(i, j), std::max(i, j));
  }
  out.push_back({"general_sparse", precision_from_edges(n, sparse),
                 ql::LaplaceBackendKind::SparseLDLT});

  std::vector<std::pair<int, int>> nearly_dense;
  std::vector<std::pair<int, int>> dense;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      dense.emplace_back(i, j);
      if ((i + 3 * j) % 10 != 0)
        nearly_dense.emplace_back(i, j);
    }
  }
  out.push_back({"nearly_dense", precision_from_edges(n, nearly_dense, true),
                 ql::LaplaceBackendKind::DenseLDLT});
  out.push_back({"dense", precision_from_edges(n, dense, true),
                 ql::LaplaceBackendKind::DenseLDLT});
  return out;
}

double dense_logdet(const Eigen::SparseMatrix<double> &matrix) {
  const Eigen::MatrixXd dense(matrix);
  Eigen::LDLT<Eigen::MatrixXd> factor(dense);
  if (factor.info() != Eigen::Success || !factor.isPositive())
    throw std::runtime_error("dense reference factorization failed");
  return factor.vectorD().array().log().sum();
}

void emit(const ModelCase &model, const std::string &phase, int sample,
          const ql::BackendRecommendation &recommendation, double elapsed_ms,
          bool success, const std::string &message = "") {
  std::cout << std::setprecision(17)
            << "{\"schema_version\":1,\"benchmark\":"
               "\"model_form_backend\",\"model_form\":\""
            << model.form << "\",\"phase\":\"" << phase
            << "\",\"sample\":" << sample
            << ",\"dimension\":" << model.hessian.rows()
            << ",\"hessian_nnz\":" << model.hessian.nonZeros()
            << ",\"bandwidth\":" << recommendation.bandwidth
            << ",\"backend\":\"" << ql::ToString(recommendation.backend)
            << "\",\"elapsed_ms\":" << elapsed_ms
            << ",\"success\":" << (success ? "true" : "false")
            << ",\"message\":\"" << message << "\"}\n";
}

} // namespace

int main(int argc, char **argv) {
  const int n = argc > 1 ? std::stoi(argv[1]) : 128;
  const int repetitions = argc > 2 ? std::stoi(argv[2]) : 10;
  if (n < 16 || repetitions < 1)
    throw std::invalid_argument("N must be >= 16 and repetitions >= 1");

  ql::StructureDetectorOptions options;
  options.prefer_dense_for_small_matrices = false;
  options.dense_size_cutoff = 0;
  options.banded_width_cutoff = 8;
  options.dense_fill_ratio = 0.75;

  bool all_ok = true;
  for (const ModelCase &model : catalog(n)) {
    ql::BackendRecommendation recommendation;
    const auto cold_start = Clock::now();
    auto backend = ql::CreateLaplaceBackendForHessian(model.hessian,
                                                      &recommendation, options);
    backend->factorize(model.hessian);
    const double cold_ms =
        std::chrono::duration<double, std::milli>(Clock::now() - cold_start)
            .count();
    const double difference =
        std::abs(backend->logdet() - dense_logdet(model.hessian));
    const bool correct = recommendation.backend == model.expected &&
                         backend->is_spd() && difference <= 1.0e-9;
    all_ok = all_ok && correct;
    emit(model, "cold_total", 0, recommendation, cold_ms, correct,
         correct ? "" : "backend or logdet mismatch");

    for (int sample = 0; sample < repetitions; ++sample) {
      const auto start = Clock::now();
      backend->factorize(model.hessian);
      const double elapsed_ms =
          std::chrono::duration<double, std::milli>(Clock::now() - start)
              .count();
      emit(model, "factorization", sample, recommendation, elapsed_ms,
           backend->is_spd());
    }
  }
  return all_ok ? 0 : 2;
}
