#include "include/quadra/stats/laplace.hpp"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <sys/resource.h>

DECLARE_ADGRAPH();

namespace ql = quadra::laplace;
using Clock = std::chrono::steady_clock;

namespace {

struct QuadraticLaplaceModel {
  int random_size = 0;
  std::vector<std::pair<int, int>> edges;

  void initialize(quadra::ModelReportContext &context) { context.clear(); }

  template <class T>
  T evaluate(const std::vector<T> &parameters,
             quadra::ModelReportContext &) const {
    const T theta = parameters[0];
    T objective = T(0.0);
    for (int i = 0; i < random_size; ++i) {
      const T &u = parameters[static_cast<std::size_t>(i + 1)];
      objective += T(0.5) * u * u + T(0.01) * theta * u;
    }
    for (const auto &edge : edges) {
      const T difference =
          parameters[static_cast<std::size_t>(edge.first + 1)] -
          parameters[static_cast<std::size_t>(edge.second + 1)];
      objective += T(0.125) * difference * difference;
    }
    objective += T(0.5) * theta * theta;
    return objective;
  }
};

struct ModelCase {
  std::string form;
  QuadraticLaplaceModel model;
  ql::LaplaceBackendKind expected_backend;
};

std::vector<std::pair<int, int>> band_edges(int n, int width) {
  std::vector<std::pair<int, int>> edges;
  for (int i = 0; i < n; ++i)
    for (int d = 1; d <= width && i + d < n; ++d)
      edges.emplace_back(i, i + d);
  return edges;
}

std::vector<ModelCase> catalog(int n) {
  std::vector<ModelCase> out;
  out.push_back({"diagonal", {n, {}}, ql::LaplaceBackendKind::Diagonal});
  out.push_back({"tridiagonal", {n, band_edges(n, 1)},
                 ql::LaplaceBackendKind::Tridiagonal});
  out.push_back({"banded", {n, band_edges(n, 4)},
                 ql::LaplaceBackendKind::Banded});

  std::vector<std::pair<int, int>> blocks;
  constexpr int block_size = 8;
  for (int start = 0; start < n; start += block_size)
    for (int i = start; i < std::min(n, start + block_size); ++i)
      for (int j = i + 1; j < std::min(n, start + block_size); ++j)
        blocks.emplace_back(i, j);
  out.push_back({"block_diagonal", {n, blocks},
                 ql::LaplaceBackendKind::Banded});

  std::vector<std::pair<int, int>> sparse = band_edges(n, 1);
  for (int i = 0; i < n / 3; i += 7) {
    const int j = n - 1 - i;
    if (i != j)
      sparse.emplace_back(std::min(i, j), std::max(i, j));
  }
  out.push_back({"general_sparse", {n, sparse},
                 ql::LaplaceBackendKind::SparseLDLT});

  std::vector<std::pair<int, int>> dense;
  for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j)
      dense.emplace_back(i, j);
  out.push_back(
      {"dense", {n, dense}, ql::LaplaceBackendKind::DenseLDLT});
  return out;
}

std::size_t peak_rss_bytes() {
  rusage usage{};
  if (getrusage(RUSAGE_SELF, &usage) != 0)
    return 0;
#ifdef __APPLE__
  return static_cast<std::size_t>(usage.ru_maxrss);
#else
  return static_cast<std::size_t>(usage.ru_maxrss) * 1024u;
#endif
}

quadra::ParameterPartition partition(int n) {
  quadra::ParameterPartition out;
  out.fixed_indices_m.push_back(0);
  for (int i = 0; i < n; ++i)
    out.random_indices_m.push_back(static_cast<std::size_t>(i + 1));
  return out;
}

void emit(const ModelCase &model, const std::string &phase, int sample,
          const quadra::LaplaceObjectiveResult &result, double elapsed_ms,
          bool success, const std::string &message = "") {
  std::cout << std::setprecision(17)
            << "{\"schema_version\":1,\"benchmark\":"
               "\"laplace_model_catalog\",\"model_form\":\""
            << model.form << "\",\"phase\":\"" << phase
            << "\",\"sample\":" << sample << ",\"dimension\":"
            << model.model.random_size << ",\"hessian_nnz\":"
            << result.hessian_random_m.nonZeros() << ",\"bandwidth\":"
            << result.backend_m.bandwidth << ",\"backend\":\""
            << ql::ToString(result.backend_m.backend)
            << "\",\"elapsed_ms\":" << elapsed_ms
            << ",\"objective\":" << result.laplace_objective_m
            << ",\"gradient_norm\":" << result.gradient_norm_random_m
            << ",\"tape_rebuilds\":" << (result.tape_rebuilt_m ? 1 : 0)
            << ",\"peak_rss_bytes\":" << peak_rss_bytes()
            << ",\"success\":" << (success ? "true" : "false")
            << ",\"message\":\"" << message << "\"}\n";
}

} // namespace

int main(int argc, char **argv) {
  const int n = argc > 1 ? std::stoi(argv[1]) : 64;
  const int repetitions = argc > 2 ? std::stoi(argv[2]) : 10;
  const std::string requested_form = argc > 3 ? argv[3] : "all";
  if (n < 16 || repetitions < 1)
    throw std::invalid_argument("N must be >= 16 and repetitions >= 1");

  bool all_ok = true;
  bool ran_model = false;
  for (ModelCase model_case : catalog(n)) {
    if (requested_form != "all" && requested_form != model_case.form)
      continue;
    ran_model = true;
    std::vector<double> random(static_cast<std::size_t>(n), 0.25);
    const std::vector<double> fixed{0.1};
    quadra::LaplaceObjectiveOptions options;
    options.newton_m.gradient_tolerance_m = 1.0e-10;
    options.newton_m.max_iterations_m = 20;
    options.compute_mixed_derivatives_m = false;
    quadra::stats::LaplaceEvaluator<QuadraticLaplaceModel> evaluator(
        model_case.model, random, partition(n), options);

    std::ostringstream library_diagnostics;
    std::streambuf *normal_output = std::cout.rdbuf(library_diagnostics.rdbuf());
    const auto cold_start = Clock::now();
    quadra::LaplaceObjectiveResult result = evaluator.evaluate(fixed);
    const double cold_ms = std::chrono::duration<double, std::milli>(
                               Clock::now() - cold_start)
                               .count();
    std::cout.rdbuf(normal_output);
    const bool cold_ok = result.converged_m && result.logdet_ok_m &&
                         std::isfinite(result.laplace_objective_m) &&
                         result.gradient_norm_random_m <= 1.0e-8 &&
                         result.backend_m.backend ==
                             model_case.expected_backend;
    all_ok = all_ok && cold_ok;
    emit(model_case, "cold_total", 0, result, cold_ms, cold_ok,
         cold_ok ? "" : "Laplace convergence or backend mismatch");
    emit(model_case, "record", 0, result, result.tape_setup_ms_m, cold_ok);
    emit(model_case, "mode_solve", 0, result, result.mode_solve_ms_m,
         cold_ok);
    emit(model_case, "logdet", 0, result, result.logdet_ms_m, cold_ok);

    for (int sample = 0; sample < repetitions; ++sample) {
      normal_output = std::cout.rdbuf(library_diagnostics.rdbuf());
      const auto start = Clock::now();
      result = evaluator.evaluate(fixed);
      const double elapsed_ms = std::chrono::duration<double, std::milli>(
                                    Clock::now() - start)
                                    .count();
      std::cout.rdbuf(normal_output);
      const bool warm_ok = result.converged_m && result.logdet_ok_m &&
                           !result.tape_rebuilt_m &&
                           result.gradient_norm_random_m <= 1.0e-8 &&
                           result.backend_m.backend ==
                               model_case.expected_backend;
      all_ok = all_ok && warm_ok;
      emit(model_case, "warm_total", sample, result, elapsed_ms, warm_ok,
           warm_ok ? "" : "warm replay or backend mismatch");
      emit(model_case, "mode_solve", sample + 1, result,
           result.mode_solve_ms_m, warm_ok);
      emit(model_case, "logdet", sample + 1, result, result.logdet_ms_m,
           warm_ok);
    }
  }
  if (!ran_model)
    throw std::invalid_argument("unknown model form: " + requested_form);
  return all_ok ? 0 : 2;
}
