#include "../core/laplace/model_analysis_report.hpp"
#include "../core/model/quadra_model.hpp"
#include "../include/quadra/stats.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

DECLARE_ADGRAPH();

class TmbStyleGaussianStateSpace
    : public quadra::QuadraModel<TmbStyleGaussianStateSpace> {
public:
  explicit TmbStyleGaussianStateSpace(int n)
      : observations_(static_cast<std::size_t>(n)) {
    parameters_m.add("unconstrained_phi", 0.35,
                     quadra::ParameterTransform::Identity, false);
    parameters_m.add("log_process_sd", std::log(0.7),
                     quadra::ParameterTransform::Identity, false);
    parameters_m.add("log_observation_sd", std::log(0.45),
                     quadra::ParameterTransform::Identity, false);
    parameters_m.add("mean", 0.1, quadra::ParameterTransform::Identity, false);
    for (int i = 0; i < n; ++i) {
      parameters_m.add("state_" + std::to_string(i), 0.0,
                       quadra::ParameterTransform::Identity, true);
      const double time = static_cast<double>(i + 1);
      observations_[static_cast<std::size_t>(i)] =
          0.35 * std::sin(0.031 * time) + 0.2 * std::cos(0.077 * time);
    }
  }

  std::vector<std::string> parameter_names_impl() const {
    return parameters_m.names();
  }
  const quadra::ParameterSet &parameters() const { return parameters_m; }

  template <class T>
  T evaluate_impl(const std::vector<T> &parameters,
                  quadra::ModelReportContext &) const {
    const T phi = quadra::stats::correlation_from_unconstrained(parameters[0]);
    using std::exp;
    const T process_sd = exp(parameters[1]);
    const T observation_sd = exp(parameters[2]);
    const T mean = parameters[3];
    std::vector<T> states(parameters.begin() + 4, parameters.end());
    T nll = quadra::stats::ar1_stationary_nll(states, mean, phi, process_sd);
    for (std::size_t i = 0; i < observations_.size(); ++i) {
      nll += quadra::stats::normal_nll(T(observations_[i]), states[i],
                                       observation_sd);
    }
    return nll;
  }

private:
  std::vector<double> observations_;
  quadra::ParameterSet parameters_m;
};

template <class Function> double milliseconds(Function &&function) {
  const auto start = std::chrono::steady_clock::now();
  function();
  return std::chrono::duration<double, std::milli>(
             std::chrono::steady_clock::now() - start)
      .count();
}

double max_abs_difference(const std::vector<double> &left,
                          const std::vector<double> &right) {
  if (left.size() != right.size()) {
    return std::numeric_limits<double>::infinity();
  }
  double difference = 0.0;
  for (std::size_t i = 0; i < left.size(); ++i) {
    difference = std::max(difference, std::abs(left[i] - right[i]));
  }
  return difference;
}

int main() {
  std::cout << std::setprecision(17) << std::unitbuf;
  std::cout
      << "model,n,requested_workers,actual_workers,topology_owners,"
         "operation_owners,operation_count,avoided_operation_bytes,"
         "avoided_edge_destination_bytes,"
         "active_directions,setup_ms,"
         "warm_exact_ms,objective_phase_ms,factorization_phase_ms,"
         "sensitivity_phase_ms,hdot_phase_ms,trace_phase_ms,internal_total_ms,"
         "warm_speedup,hdot_speedup,setup_ratio,"
         "hessian_nnz,hessian_density,bandwidth,structure,backend,solver,"
         "complexity,symbolic_reuse,objective_rebuilds,hdot_rebuilds,"
         "objective_difference,gradient_max_abs_difference,success\n";

  const std::vector<double> fixed = {0.35, std::log(0.7), std::log(0.45), 0.1};
  for (int n : {100, 300, 1000}) {
    TmbStyleGaussianStateSpace model(n);
    std::vector<double> random(static_cast<std::size_t>(n), 0.0);
    quadra::LaplaceObjectiveOptions objective_options;
    objective_options.include_constant_m = true;
    objective_options.newton_m.gradient_tolerance_m = 1.0e-9;
    objective_options.newton_m.step_tolerance_m = 1.0e-11;

    std::vector<double> reference_gradient;
    double reference_objective = std::nan("");
    double serial_setup_ms = std::nan("");
    double serial_warm_ms = std::nan("");
    double serial_hdot_ms = std::nan("");
    for (int requested_workers : {1, 2, 0}) {
      quadra::laplace::ExactLaplaceGradientEngineOptions engine_options;
      engine_options.hdot_workers = requested_workers;
      std::unique_ptr<
          quadra::stats::ExactLaplaceEvaluator<TmbStyleGaussianStateSpace>>
          evaluator;
      const double setup_ms = milliseconds([&]() {
        evaluator.reset(new quadra::stats::ExactLaplaceEvaluator<
                        TmbStyleGaussianStateSpace>(
            model, fixed, random, model.parameters(), objective_options,
            engine_options));
      });

      quadra::stats::ExactLaplaceResult result;
      const int repetitions = n <= 100 ? 30 : (n <= 300 ? 15 : 8);
      const double total_ms = milliseconds([&]() {
        for (int repetition = 0; repetition < repetitions; ++repetition) {
          result = evaluator->evaluate(fixed);
        }
      });
      const double warm_ms = total_ms / static_cast<double>(repetitions);

      const auto diagnostics = quadra::laplace::analyze_hessian_structure(
          result.objective.hessian_random_m);
      if (requested_workers == 1) {
        reference_gradient = result.gradient;
        reference_objective = result.objective.laplace_objective_m;
        serial_setup_ms = setup_ms;
        serial_warm_ms = warm_ms;
        serial_hdot_ms = result.timings.hdot_ms;
      }
      const double objective_difference =
          std::abs(result.objective.laplace_objective_m - reference_objective);
      const double gradient_difference =
          max_abs_difference(result.gradient, reference_gradient);
      const bool success =
          result.success && result.active_directions.size() == 3 &&
          objective_difference < 1.0e-10 && gradient_difference < 1.0e-10 &&
          diagnostics.is_tridiagonal() && diagnostics.bandwidth == 1;

      std::cout << "gaussian_ar1_state_space," << n << "," << requested_workers
                << "," << evaluator->hdot_worker_count() << ","
                << evaluator->hdot_shared_topology_owner_count() << ","
                << evaluator->hdot_shared_operation_owner_count() << ","
                << evaluator->hdot_operation_count() << ","
                << (evaluator->hdot_worker_count() - 1) *
                       evaluator->hdot_operation_count() *
                       sizeof(had::ADOperation)
                << ","
                << (evaluator->hdot_worker_count() - 1) *
                       evaluator->hdot_operation_count() * 2 *
                       sizeof(had::VertexId)
                << "," << result.active_directions.size() << "," << setup_ms
                << "," << warm_ms << "," << result.timings.objective_ms << ","
                << result.timings.factorization_ms << ","
                << result.timings.mode_sensitivity_ms << ","
                << result.timings.hdot_ms << "," << result.timings.trace_ms
                << "," << result.timings.total_ms << ","
                << serial_warm_ms / warm_ms << ","
                << serial_hdot_ms / result.timings.hdot_ms << ","
                << setup_ms / serial_setup_ms << "," << diagnostics.nnz << ","
                << diagnostics.fill_ratio << "," << diagnostics.bandwidth << ","
                << quadra::laplace::ToString(diagnostics.structure) << ","
                << quadra::laplace::ToString(diagnostics.backend) << ","
                << quadra::laplace::ToString(diagnostics.solver) << ","
                << diagnostics.complexity << ","
                << (diagnostics.supports_symbolic_reuse ? 1 : 0) << ","
                << evaluator->objective_tape_rebuild_count() << ","
                << evaluator->hdot_tape_rebuild_count() << ","
                << objective_difference << "," << gradient_difference << ","
                << (success ? 1 : 0) << "\n";
      if (!success) {
        return 1;
      }
    }
  }
  return 0;
}
