#pragma once

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "laplace_objective_cached.hpp"

namespace quadra {

struct LaplaceEvaluatorOptions {
    CachedLaplaceObjectiveOptions objective_m;
    bool warm_start_random_m = true;
};

template <class Model>
class LaplaceEvaluator {
public:
    LaplaceEvaluator(
        Model& model,
        const ParameterSet& parameters,
        std::vector<double> random_initial,
        LaplaceEvaluatorOptions options = LaplaceEvaluatorOptions()
    )
        : model_m(model),
          partition_m(partition_parameters(parameters)),
          random_start_m(std::move(random_initial)),
          options_m(options)
    {
        if (random_start_m.size() != partition_m.random_indices_m.size()) {
            throw std::invalid_argument(
                "LaplaceEvaluator: random_initial has incorrect length."
            );
        }
    }

    LaplaceObjectiveResult evaluate(const std::vector<double>& fixed) {
        if (fixed.size() != partition_m.fixed_indices_m.size()) {
            throw std::invalid_argument(
                "LaplaceEvaluator::evaluate fixed vector has incorrect length."
            );
        }

        LaplaceObjectiveResult result =
            evaluate_laplace_objective_cached(
                model_m,
                fixed,
                random_start_m,
                partition_m,
                cache_state_m,
                options_m.objective_m
            );

        last_result_m = result;
        has_last_result_m = true;
        evaluations_m += 1;

        if (options_m.warm_start_random_m && result.converged_m) {
            random_start_m = result.u_hat_m;
        }

        return result;
    }

    const ParameterPartition& partition() const {
        return partition_m;
    }

    const std::vector<double>& random_start() const {
        return random_start_m;
    }

    void set_random_start(const std::vector<double>& random_start) {
        if (random_start.size() != partition_m.random_indices_m.size()) {
            throw std::invalid_argument(
                "LaplaceEvaluator::set_random_start length mismatch."
            );
        }

        random_start_m = random_start;
    }

    bool has_last_result() const {
        return has_last_result_m;
    }

    const LaplaceObjectiveResult& last_result() const {
        if (!has_last_result_m) {
            throw std::runtime_error(
                "LaplaceEvaluator::last_result called before evaluate."
            );
        }

        return last_result_m;
    }

    CachedLaplaceObjectiveState& cache_state() {
        return cache_state_m;
    }

    const CachedLaplaceObjectiveState& cache_state() const {
        return cache_state_m;
    }

    int evaluations() const {
        return evaluations_m;
    }

private:
    Model& model_m;
    ParameterPartition partition_m;
    std::vector<double> random_start_m;

    LaplaceEvaluatorOptions options_m;
    CachedLaplaceObjectiveState cache_state_m;

    LaplaceObjectiveResult last_result_m;
    bool has_last_result_m = false;

    int evaluations_m = 0;
};

} // namespace quadra
