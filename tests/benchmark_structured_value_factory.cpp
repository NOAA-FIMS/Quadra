#include "../core/laplace/hessian_structure.hpp"
#include "../core/laplace/structured_value_factory.hpp"

#include <Eigen/Sparse>

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using quadra::laplace::BackendRecommendation;
using quadra::laplace::extract_structured_values;
using quadra::laplace::logdet_structured_values;
using quadra::laplace::LogDetSparseLDLT;
using quadra::laplace::StructureDetector;
using quadra::laplace::StructureDetectorOptions;
using quadra::laplace::StructureInfo;

using Clock = std::chrono::high_resolution_clock;

double ms_between(const Clock::time_point &a, const Clock::time_point &b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

std::vector<int> parse_lengths(const std::string &s) {
  std::vector<int> out;
  std::stringstream ss(s);
  std::string item;

  while (std::getline(ss, item, ',')) {
    if (!item.empty())
      out.push_back(std::stoi(item));
  }

  return out;
}

Eigen::SparseMatrix<double> make_diagonal(const int n) {
  std::vector<Eigen::Triplet<double>> t;
  t.reserve(static_cast<std::size_t>(n));

  for (int i = 0; i < n; ++i) {
    t.emplace_back(i, i, 2.0 + 0.001 * i);
  }

  Eigen::SparseMatrix<double> H(n, n);
  H.setFromTriplets(t.begin(), t.end());
  H.makeCompressed();
  return H;
}

Eigen::SparseMatrix<double> make_tridiagonal(const int n) {
  std::vector<Eigen::Triplet<double>> t;
  t.reserve(static_cast<std::size_t>(3 * n));

  for (int i = 0; i < n; ++i) {
    t.emplace_back(i, i, 4.0 + 0.001 * i);

    if (i > 0) {
      const double e = -0.20 + 0.00001 * (i % 17);
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
  t.reserve(static_cast<std::size_t>(n * (2 * bandwidth + 1)));

  // Strictly diagonally dominant SPD-ish symmetric banded matrix.
  for (int i = 0; i < n; ++i) {
    double diag = 10.0 + 0.001 * i;

    for (int d = 1; d <= bandwidth; ++d) {
      if (i - d < 0)
        continue;

      const double e = ((d % 2 == 0) ? 0.015 : -0.025) / static_cast<double>(d);

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

StructureDetector make_detector() {
  StructureDetectorOptions opts;
  opts.prefer_dense_for_small_matrices = false;
  opts.dense_size_cutoff = 0;
  opts.banded_width_cutoff = 64;
  opts.dense_fill_ratio = 0.75;
  return StructureDetector(opts);
}

struct BenchResult {
  std::string name;
  int n = 0;
  int nnz = 0;
  int bandwidth = 0;
  double sparse_ms = 0.0;
  double structured_ms = 0.0;
  double speedup = 0.0;
  double logdet_diff = 0.0;
};

BenchResult bench_case(const std::string &name,
                       const Eigen::SparseMatrix<double> &H, const int reps) {
  StructureDetector detector = make_detector();
  const BackendRecommendation rec = detector.Analyze(H);
  const StructureInfo info = quadra::laplace::InspectHessianStructure(
      H, detector.options().structure_options);

  volatile double sparse_acc = 0.0;
  volatile double structured_acc = 0.0;

  const double sparse_ref = LogDetSparseLDLT(H);
  const auto values = extract_structured_values(H, rec);
  const double structured_ref = logdet_structured_values(values);

  const auto sparse0 = Clock::now();
  for (int r = 0; r < reps; ++r) {
    sparse_acc += LogDetSparseLDLT(H);
  }
  const auto sparse1 = Clock::now();

  const auto structured0 = Clock::now();
  for (int r = 0; r < reps; ++r) {
    const auto values_r = extract_structured_values(H, rec);
    structured_acc += logdet_structured_values(values_r);
  }
  const auto structured1 = Clock::now();

  (void)sparse_acc;
  (void)structured_acc;

  BenchResult out;
  out.name = name;
  out.n = static_cast<int>(H.rows());
  out.nnz = info.nnz;
  out.bandwidth = info.max_bandwidth;
  out.sparse_ms = ms_between(sparse0, sparse1) / static_cast<double>(reps);
  out.structured_ms =
      ms_between(structured0, structured1) / static_cast<double>(reps);
  out.speedup =
      out.structured_ms > 0.0 ? out.sparse_ms / out.structured_ms : 0.0;
  out.logdet_diff = structured_ref - sparse_ref;
  return out;
}

int main(int argc, char **argv) {
  int reps = 20;
  std::vector<int> lengths = {100, 500, 1000, 5000};

  if (argc > 1)
    reps = std::stoi(argv[1]);
  if (argc > 2)
    lengths = parse_lengths(argv[2]);

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "Structured value factory benchmark\n";
  std::cout << "reps per case = " << reps << "\n\n";

  std::cout << std::setw(14) << "case" << std::setw(8) << "n" << std::setw(10)
            << "nnz" << std::setw(10) << "band" << std::setw(14) << "sparse_ms"
            << std::setw(16) << "structured_ms" << std::setw(12) << "speedup"
            << std::setw(16) << "logdet_diff" << "\n";

  for (const int n : lengths) {
    const std::vector<std::pair<std::string, Eigen::SparseMatrix<double>>>
        cases = {
            {"diagonal", make_diagonal(n)},
            {"tridiagonal", make_tridiagonal(n)},
            {"banded2", make_banded(n, 2)},
            {"banded5", make_banded(n, 5)},
            {"banded10", make_banded(n, 10)},
        };

    for (const auto &item : cases) {
      const BenchResult b = bench_case(item.first, item.second, reps);

      std::cout << std::setw(14) << b.name << std::setw(8) << b.n
                << std::setw(10) << b.nnz << std::setw(10) << b.bandwidth
                << std::setw(14) << b.sparse_ms << std::setw(16)
                << b.structured_ms << std::setw(12) << b.speedup
                << std::setw(16) << b.logdet_diff << "\n";
    }
  }

  return 0;
}
