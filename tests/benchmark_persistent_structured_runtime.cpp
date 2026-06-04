#include "../core/laplace/persistent_structured_runtime.hpp"
#include "../core/laplace/structured_value_factory.hpp"
#include "../core/laplace/hessian_structure.hpp"

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
using quadra::laplace::PersistentStructuredLaplaceRuntime;
using quadra::laplace::PersistentStructuredRuntimeState;
using quadra::laplace::StructureDetector;
using quadra::laplace::StructureDetectorOptions;
using quadra::laplace::StructureInfo;
using quadra::laplace::extract_structured_values;
using quadra::laplace::logdet_structured_values;

using Clock = std::chrono::high_resolution_clock;

double ms_between(const Clock::time_point& a, const Clock::time_point& b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

std::vector<int> parse_lengths(const std::string& s) {
  std::vector<int> out;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, ',')) {
    if (!item.empty()) out.push_back(std::stoi(item));
  }
  return out;
}

StructureDetectorOptions detector_options() {
  StructureDetectorOptions opts;
  opts.prefer_dense_for_small_matrices = false;
  opts.dense_size_cutoff = 0;
  opts.banded_width_cutoff = 64;
  opts.dense_fill_ratio = 0.75;
  return opts;
}

Eigen::SparseMatrix<double> make_diagonal(const int n, const double scale) {
  std::vector<Eigen::Triplet<double>> t;
  t.reserve(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    t.emplace_back(i, i, scale * (2.0 + 0.001 * i));
  }
  Eigen::SparseMatrix<double> H(n, n);
  H.setFromTriplets(t.begin(), t.end());
  H.makeCompressed();
  return H;
}

Eigen::SparseMatrix<double> make_tridiagonal(const int n, const double scale) {
  std::vector<Eigen::Triplet<double>> t;
  t.reserve(static_cast<std::size_t>(3 * n));
  for (int i = 0; i < n; ++i) {
    t.emplace_back(i, i, scale * (4.0 + 0.001 * i));
    if (i > 0) {
      const double e = scale * (-0.20 + 0.00001 * (i % 17));
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
  t.reserve(static_cast<std::size_t>(n * (2 * bandwidth + 1)));

  for (int i = 0; i < n; ++i) {
    double diag = scale * (10.0 + 0.001 * i);

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

Eigen::SparseMatrix<double> make_case(const std::string& name,
                                      const int n,
                                      const double scale) {
  if (name == "diagonal") return make_diagonal(n, scale);
  if (name == "tridiagonal") return make_tridiagonal(n, scale);
  if (name == "banded2") return make_banded(n, 2, scale);
  if (name == "banded5") return make_banded(n, 5, scale);
  if (name == "banded10") return make_banded(n, 10, scale);
  throw std::invalid_argument("unknown benchmark case: " + name);
}

struct BenchResult {
  std::string name;
  int n = 0;
  int nnz = 0;
  int bandwidth = 0;

  double full_detect_ms = 0.0;
  double runtime_first_ms = 0.0;
  double runtime_update_ms = 0.0;
  double state_update_ms = 0.0;

  double runtime_speedup = 0.0;
  double state_speedup = 0.0;
  double logdet_diff = 0.0;
};

BenchResult bench_case(const std::string& name,
                       const int n,
                       const int reps) {
  const StructureDetectorOptions opts = detector_options();
  StructureDetector detector(opts);

  const Eigen::SparseMatrix<double> H0 = make_case(name, n, 1.0);
  const StructureInfo info =
      quadra::laplace::InspectHessianStructure(
          H0, opts.structure_options);

  // Reference for equality checks.
  const BackendRecommendation rec0 = detector.Analyze(H0);
  const auto values0 = extract_structured_values(H0, rec0);
  const double reference_logdet = logdet_structured_values(values0);

  volatile double full_acc = 0.0;
  volatile double runtime_acc = 0.0;
  volatile double state_acc = 0.0;

  // A. Full path each call.
  const auto full0 = Clock::now();
  for (int r = 0; r < reps; ++r) {
    const double scale = 1.0 + 1e-6 * static_cast<double>(r + 1);
    const Eigen::SparseMatrix<double> H = make_case(name, n, scale);
    const BackendRecommendation rec = detector.Analyze(H);
    const auto values = extract_structured_values(H, rec);
    full_acc += logdet_structured_values(values);
  }
  const auto full1 = Clock::now();

  // B/C. Persistent runtime.
  PersistentStructuredLaplaceRuntime runtime(opts);

  const auto first0 = Clock::now();
  const auto first_result = runtime.evaluate(H0);
  const auto first1 = Clock::now();

  const auto runtime0 = Clock::now();
  for (int r = 0; r < reps; ++r) {
    const double scale = 1.0 + 1e-6 * static_cast<double>(r + 1);
    const Eigen::SparseMatrix<double> H = make_case(name, n, scale);
    runtime_acc += runtime.evaluate(H).logdet;
  }
  const auto runtime1 = Clock::now();

  // D. State update_values_only directly.
  PersistentStructuredRuntimeState state;
  state.update_from_hessian(H0, opts);

  const auto state0 = Clock::now();
  for (int r = 0; r < reps; ++r) {
    const double scale = 1.0 + 1e-6 * static_cast<double>(r + 1);
    const Eigen::SparseMatrix<double> H = make_case(name, n, scale);
    state.update_values_only(H);
    state_acc += state.logdet();
  }
  const auto state1 = Clock::now();

  (void)full_acc;
  (void)runtime_acc;
  (void)state_acc;

  BenchResult out;
  out.name = name;
  out.n = n;
  out.nnz = info.nnz;
  out.bandwidth = info.max_bandwidth;
  out.full_detect_ms = ms_between(full0, full1) / static_cast<double>(reps);
  out.runtime_first_ms = ms_between(first0, first1);
  out.runtime_update_ms =
      ms_between(runtime0, runtime1) / static_cast<double>(reps);
  out.state_update_ms =
      ms_between(state0, state1) / static_cast<double>(reps);
  out.runtime_speedup = out.runtime_update_ms > 0.0
                            ? out.full_detect_ms / out.runtime_update_ms
                            : 0.0;
  out.state_speedup = out.state_update_ms > 0.0
                          ? out.full_detect_ms / out.state_update_ms
                          : 0.0;
  out.logdet_diff = first_result.logdet - reference_logdet;
  return out;
}

int main(int argc, char** argv) {
  int reps = 20;
  std::vector<int> lengths = {100, 500, 1000, 5000};

  if (argc > 1) reps = std::stoi(argv[1]);
  if (argc > 2) lengths = parse_lengths(argv[2]);

  const std::vector<std::string> cases = {
      "diagonal", "tridiagonal", "banded2", "banded5", "banded10"};

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "Persistent structured runtime benchmark\n";
  std::cout << "reps per case = " << reps << "\n\n";

  std::cout << std::setw(14) << "case"
            << std::setw(8) << "n"
            << std::setw(10) << "nnz"
            << std::setw(8) << "band"
            << std::setw(16) << "full_ms"
            << std::setw(16) << "first_ms"
            << std::setw(18) << "runtime_ms"
            << std::setw(16) << "state_ms"
            << std::setw(16) << "runtime_x"
            << std::setw(14) << "state_x"
            << std::setw(16) << "logdet_diff"
            << "\n";

  for (const int n : lengths) {
    for (const std::string& name : cases) {
      const BenchResult b = bench_case(name, n, reps);

      std::cout << std::setw(14) << b.name
                << std::setw(8) << b.n
                << std::setw(10) << b.nnz
                << std::setw(8) << b.bandwidth
                << std::setw(16) << b.full_detect_ms
                << std::setw(16) << b.runtime_first_ms
                << std::setw(18) << b.runtime_update_ms
                << std::setw(16) << b.state_update_ms
                << std::setw(16) << b.runtime_speedup
                << std::setw(14) << b.state_speedup
                << std::setw(16) << b.logdet_diff
                << "\n";
    }
  }

  return 0;
}
