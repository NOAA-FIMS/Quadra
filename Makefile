UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
QUADRA_TIME_CMD := $(QUADRA_TIME_CMD)
else
QUADRA_TIME_CMD := /usr/bin/time -v
endif

include config/quadra_includes.mk

CXX ?= $(QUADRA_CXX)
CXXFLAGS ?= $(QUADRA_CXXFLAGS) $(QUADRA_INCLUDE_FLAGS)

TESTS = \
	tests/test_curvature_depends_on_theta \
	tests/test_poisson_random_effect \
	tests/test_ar1_random_walk \
	tests/test_hdot_validation \
	tests/test_scalar_generic_model \
	tests/test_parameter_transforms

EXAMPLES = \
	examples/fisheries_random_year_effects \
	examples/fisheries_age_selectivity_random_walk \
	examples/fisheries_index_cpue_laplace \
	examples/normal_mean_ad \
	examples/transformed_parameters

all: $(TESTS) $(EXAMPLES)

tests/%: tests/%.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

examples: $(EXAMPLES)

examples/%: examples/%.cpp

	$(CXX) $(CXXFLAGS) -o $@ $<

run-examples: examples

	@for ex in $(EXAMPLES); do \
		echo ""; \
		echo "======================================"; \
		echo "Running $$ex"; \
		echo "======================================"; \
		$(QUADRA_TIME_CMD) ./$$ex; \
	done

clean-examples:

	rm -f $(EXAMPLES)

validate-hdot:
	$(CXX) $(CXXFLAGS) -DQUADRA_VALIDATE_HDOT -o tests/test_hdot_validation tests/test_hdot_validation.cpp
	./tests/test_hdot_validation
	./tests/test_scalar_generic_model
	./tests/test_parameter_transforms
	./tests/test_transform_jacobian
	./tests/test_parameter_partition
	./tests/test_random_effect_objective
	./tests/test_random_effect_hessian
	./tests/test_random_effect_newton
	./tests/test_laplace_objective
	./tests/test_laplace_fixed_gradient
	./tests/test_laplace_optimizer
	./tests/test_laplace_lbfgs_optimizer
	./tests/test_laplace_exact_lbfgs_optimizer
	./tests/test_laplace_exact_objective_gradient
	./tests/test_laplace_exact_objective_lbfgs
	./tests/test_joint_only_exact_gradient
	./tests/test_had_quadra_reusable_tape_probe
	./tests/test_had_quadra_value_mutation_behavior
	./tests/test_had_quadra_adjoint_reset_probe
	./tests/test_had_quadra_forward_replay_math_ops
	./tests/test_reusable_tape_random_intercept_gradient
	./tests/test_multi_random_intercept_laplace
	./tests/test_correlated_random_intercept_hessian
	./tests/test_sparse_factorization_cache
	./tests/test_laplace_objective_cached
	./tests/test_laplace_evaluator
	./tests/test_laplace_profile
	./tests/test_laplace_exact_gradient
	./tests/test_laplace_implicit_derivatives

run-tests: $(TESTS) test-big-laplace-convergence-contract
	./tests/test_curvature_depends_on_theta
	./tests/test_poisson_random_effect
	./tests/test_ar1_random_walk
	./tests/test_hdot_validation

clean:
	rm -f $(TESTS) $(EXAMPLES)

.PHONY: all run-tests validate-hdot clean


tests/test_parameter_transforms: tests/test_parameter_transforms.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


examples/transformed_parameters: examples/transformed_parameters.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


tests/test_transform_jacobian: tests/test_transform_jacobian.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


examples/transform_jacobian_demo: examples/transform_jacobian_demo.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


tests/test_parameter_partition: tests/test_parameter_partition.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


examples/parameter_partition_demo: examples/parameter_partition_demo.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


tests/test_random_effect_objective: tests/test_random_effect_objective.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


examples/random_effect_objective_demo: examples/random_effect_objective_demo.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


tests/test_random_effect_hessian: tests/test_random_effect_hessian.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


examples/random_effect_hessian_demo: examples/random_effect_hessian_demo.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


tests/test_random_effect_newton: tests/test_random_effect_newton.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


examples/random_effect_newton_demo: examples/random_effect_newton_demo.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


tests/test_laplace_objective: tests/test_laplace_objective.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


examples/laplace_objective_demo: examples/laplace_objective_demo.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


tests/test_laplace_fixed_gradient: tests/test_laplace_fixed_gradient.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


examples/laplace_fixed_gradient_demo: examples/laplace_fixed_gradient_demo.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


tests/test_laplace_optimizer: tests/test_laplace_optimizer.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


examples/laplace_optimizer_demo: examples/laplace_optimizer_demo.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


tests/test_laplace_exact_gradient: tests/test_laplace_exact_gradient.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


examples/laplace_exact_gradient_demo: examples/laplace_exact_gradient_demo.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

tests/test_laplace_implicit_derivatives: tests/test_laplace_implicit_derivatives.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

examples/laplace_implicit_derivatives_demo: examples/laplace_implicit_derivatives_demo.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<	


.PHONY: laplace-stack-check
laplace-stack-check: \
	tests/test_random_effect_objective \
	tests/test_random_effect_hessian \
	tests/test_random_effect_newton \
	tests/test_laplace_objective \
	tests/test_laplace_fixed_gradient \
	tests/test_laplace_optimizer \
	tests/test_laplace_exact_gradient \
	tests/test_laplace_implicit_derivatives
	./tests/test_random_effect_objective
	./tests/test_random_effect_hessian
	./tests/test_random_effect_newton
	./tests/test_laplace_objective
	./tests/test_laplace_fixed_gradient
	./tests/test_laplace_optimizer
	./tests/test_laplace_exact_gradient
	./tests/test_laplace_implicit_derivatives


tests/test_laplace_lbfgs_optimizer: tests/test_laplace_lbfgs_optimizer.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


examples/laplace_lbfgs_optimizer_demo: examples/laplace_lbfgs_optimizer_demo.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


benchmarks/random_intercept_scaling: benchmarks/random_intercept_scaling.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


.PHONY: benchmark-laplace
benchmark-laplace: benchmarks/random_intercept_scaling
	./benchmarks/random_intercept_scaling


tests/test_laplace_exact_lbfgs_optimizer: tests/test_laplace_exact_lbfgs_optimizer.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


benchmarks/random_intercept_scaling_exact: benchmarks/random_intercept_scaling_exact.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


.PHONY: benchmark-laplace-exact
benchmark-laplace-exact: benchmarks/random_intercept_scaling_exact
	./benchmarks/random_intercept_scaling_exact


tests/test_laplace_profile: tests/test_laplace_profile.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


benchmarks/laplace_profile_random_intercept: benchmarks/laplace_profile_random_intercept.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


.PHONY: benchmark-laplace-profile
benchmark-laplace-profile: benchmarks/laplace_profile_random_intercept
	./benchmarks/laplace_profile_random_intercept


tests/test_laplace_exact_objective_gradient: tests/test_laplace_exact_objective_gradient.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


tests/test_laplace_exact_objective_lbfgs: tests/test_laplace_exact_objective_lbfgs.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


benchmarks/random_intercept_scaling_combined: benchmarks/random_intercept_scaling_combined.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


.PHONY: benchmark-laplace-combined
benchmark-laplace-combined: benchmarks/random_intercept_scaling_combined
	./benchmarks/random_intercept_scaling_combined


benchmarks/ad_cost_breakdown: benchmarks/ad_cost_breakdown.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

.PHONY: benchmark-ad-cost
benchmark-ad-cost: benchmarks/ad_cost_breakdown
	./benchmarks/ad_cost_breakdown


tests/test_joint_only_exact_gradient: tests/test_joint_only_exact_gradient.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


benchmarks/random_intercept_joint_only_scaling: benchmarks/random_intercept_joint_only_scaling.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


.PHONY: benchmark-joint-only
benchmark-joint-only: benchmarks/random_intercept_joint_only_scaling
	./benchmarks/random_intercept_joint_only_scaling


tests/test_had_quadra_reusable_tape_probe: tests/test_had_quadra_reusable_tape_probe.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


tests/test_had_quadra_value_mutation_behavior: tests/test_had_quadra_value_mutation_behavior.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


tests/test_had_quadra_adjoint_reset_probe: tests/test_had_quadra_adjoint_reset_probe.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

tests/test_had_quadra_zero_adjoints: tests/test_had_quadra_zero_adjoints.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

tests/test_had_quadra_forward_replay_probe: tests/test_had_quadra_forward_replay_probe.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

tests/test_forward_replay_design_prototype: tests/test_forward_replay_design_prototype.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

tests/test_had_quadra_forward_replay: tests/test_had_quadra_forward_replay.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


tests/test_had_quadra_forward_replay_math_ops: tests/test_had_quadra_forward_replay_math_ops.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


tests/test_reusable_tape_random_intercept_gradient: tests/test_reusable_tape_random_intercept_gradient.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


benchmarks/reusable_tape_random_intercept: benchmarks/reusable_tape_random_intercept.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


.PHONY: benchmark-reusable-tape
benchmark-reusable-tape: benchmarks/reusable_tape_random_intercept
	./benchmarks/reusable_tape_random_intercept


tests/test_multi_random_intercept_laplace: tests/test_multi_random_intercept_laplace.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


benchmarks/multi_random_intercept_scaling: benchmarks/multi_random_intercept_scaling.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


.PHONY: benchmark-multi-random
benchmark-multi-random: benchmarks/multi_random_intercept_scaling
	./benchmarks/multi_random_intercept_scaling


tests/test_correlated_random_intercept_hessian: tests/test_correlated_random_intercept_hessian.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


benchmarks/correlated_random_intercept_scaling: benchmarks/correlated_random_intercept_scaling.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


.PHONY: benchmark-correlated-random
benchmark-correlated-random: benchmarks/correlated_random_intercept_scaling
	./benchmarks/correlated_random_intercept_scaling


tests/test_sparse_factorization_cache: tests/test_sparse_factorization_cache.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


benchmarks/sparse_factorization_cache_benchmark: benchmarks/sparse_factorization_cache_benchmark.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


.PHONY: benchmark-sparse-factorization-cache
benchmark-sparse-factorization-cache: benchmarks/sparse_factorization_cache_benchmark
	./benchmarks/sparse_factorization_cache_benchmark


tests/test_laplace_objective_cached: tests/test_laplace_objective_cached.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


benchmarks/cached_laplace_objective_benchmark: benchmarks/cached_laplace_objective_benchmark.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


.PHONY: benchmark-cached-laplace-objective
benchmark-cached-laplace-objective: benchmarks/cached_laplace_objective_benchmark
	./benchmarks/cached_laplace_objective_benchmark


tests/test_laplace_evaluator: tests/test_laplace_evaluator.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


benchmarks/laplace_evaluator_benchmark: benchmarks/laplace_evaluator_benchmark.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<


.PHONY: benchmark-laplace-evaluator
benchmark-laplace-evaluator: benchmarks/laplace_evaluator_benchmark
	./benchmarks/laplace_evaluator_benchmark

.PHONY: benchmark-arena-pool
benchmark-arena-pool: benchmarks/arena_pool_benchmark
	./benchmarks/arena_pool_benchmark

benchmarks/arena_pool_benchmark: benchmarks/arena_pool_benchmark.cpp core/memory/arena_pool.hpp
	$(CXX) $(CXXFLAGS) -std=c++17 -O3 -I. -o benchmarks/arena_pool_benchmark benchmarks/arena_pool_benchmark.cpp

.PHONY: benchmark-had-quadra-workspace
benchmark-had-quadra-workspace: benchmarks/had_quadra_workspace_benchmark
	./benchmarks/had_quadra_workspace_benchmark

benchmarks/had_quadra_workspace_benchmark: benchmarks/had_quadra_workspace_benchmark.cpp core/autodiff/had_quadra_workspace.hpp core/memory/arena_pool.hpp
	$(CXX) $(CXXFLAGS) -std=c++17 -O3 -I. -o benchmarks/had_quadra_workspace_benchmark benchmarks/had_quadra_workspace_benchmark.cpp

.PHONY: benchmark-laplace-evaluator-workspace
benchmark-laplace-evaluator-workspace: benchmarks/laplace_evaluator_workspace_benchmark
	./benchmarks/laplace_evaluator_workspace_benchmark

benchmarks/laplace_evaluator_workspace_benchmark: benchmarks/laplace_evaluator_workspace_benchmark.cpp core/autodiff/had_quadra_workspace.hpp core/memory/arena_pool.hpp
	$(CXX) $(CXXFLAGS) -std=c++17 -O3 -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o benchmarks/laplace_evaluator_workspace_benchmark benchmarks/laplace_evaluator_workspace_benchmark.cpp

.PHONY: benchmark-had-quadra-graph-reserve
benchmark-had-quadra-graph-reserve: benchmarks/had_quadra_graph_reserve_benchmark
	./benchmarks/had_quadra_graph_reserve_benchmark

benchmarks/had_quadra_graph_reserve_benchmark: benchmarks/had_quadra_graph_reserve_benchmark.cpp core/had_quadra.hpp core/autodiff/had_quadra_graph_reserve.hpp
	$(CXX) $(CXXFLAGS) -std=c++17 -O3 -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o benchmarks/had_quadra_graph_reserve_benchmark benchmarks/had_quadra_graph_reserve_benchmark.cpp

.PHONY: benchmark-btree-accumulation
benchmark-btree-accumulation: benchmarks/btree_accumulation_benchmark
	./benchmarks/btree_accumulation_benchmark

benchmarks/btree_accumulation_benchmark: benchmarks/btree_accumulation_benchmark.cpp
	$(CXX) $(CXXFLAGS) -std=c++17 -O3 -I. -o benchmarks/btree_accumulation_benchmark benchmarks/btree_accumulation_benchmark.cpp

.PHONY: benchmark-sparse-accumulator
benchmark-sparse-accumulator: benchmarks/sparse_accumulator_benchmark
	./benchmarks/sparse_accumulator_benchmark

benchmarks/sparse_accumulator_benchmark: benchmarks/sparse_accumulator_benchmark.cpp core/sparse/sparse_accumulator.hpp
	$(CXX) $(CXXFLAGS) -std=c++17 -O3 -I. -o benchmarks/sparse_accumulator_benchmark benchmarks/sparse_accumulator_benchmark.cpp

.PHONY: benchmark-hessian-slot-accumulator
benchmark-hessian-slot-accumulator: benchmarks/hessian_slot_accumulator_benchmark
	./benchmarks/hessian_slot_accumulator_benchmark

benchmarks/hessian_slot_accumulator_benchmark: benchmarks/hessian_slot_accumulator_benchmark.cpp core/sparse/hessian_slot_accumulator.hpp
	$(CXX) $(CXXFLAGS) -std=c++17 -O3 -I. -o benchmarks/hessian_slot_accumulator_benchmark benchmarks/hessian_slot_accumulator_benchmark.cpp

.PHONY: benchmark-hessian-slot-hotloop
benchmark-hessian-slot-hotloop: benchmarks/hessian_slot_accumulator_hotloop_benchmark
	./benchmarks/hessian_slot_accumulator_hotloop_benchmark

benchmarks/hessian_slot_accumulator_hotloop_benchmark: benchmarks/hessian_slot_accumulator_hotloop_benchmark.cpp core/sparse/hessian_slot_accumulator.hpp
	$(CXX) $(CXXFLAGS) -std=c++17 -O3 -I. -o benchmarks/hessian_slot_accumulator_hotloop_benchmark benchmarks/hessian_slot_accumulator_hotloop_benchmark.cpp

.PHONY: benchmark-hessian-slot-tape-workspace
benchmark-hessian-slot-tape-workspace: benchmarks/hessian_slot_tape_workspace_benchmark
	./benchmarks/hessian_slot_tape_workspace_benchmark

benchmarks/hessian_slot_tape_workspace_benchmark: benchmarks/hessian_slot_tape_workspace_benchmark.cpp core/laplace/hessian_slot_tape_workspace.hpp core/sparse/hessian_slot_accumulator.hpp
	$(CXX) $(CXXFLAGS) -std=c++17 -O3 -I. -o benchmarks/hessian_slot_tape_workspace_benchmark benchmarks/hessian_slot_tape_workspace_benchmark.cpp

# ---- Quadra examples -------------------------------------------------------

CXX ?= clang++
CXXFLAGS ?= -std=c++17 -O3 \
	-I. \
	-I./external/eigen \
	-I./external/had \
	-I./external/LBFGSpp/include

EXAMPLES := $(patsubst examples/%.cpp,examples/%,$(wildcard examples/*.cpp))
BIG_EXAMPLES := $(patsubst examples/big/%.cpp,examples/big/%,$(wildcard examples/big/*.cpp))

examples: $(EXAMPLES)

big-examples: $(BIG_EXAMPLES)

examples/%: examples/%.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

examples/big/%: examples/big/%.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

run-examples: examples
	@for ex in $(EXAMPLES); do \
		echo ""; \
		echo "======================================"; \
		echo "Running $$ex"; \
		echo "======================================"; \
		$(QUADRA_TIME_CMD) ./$$ex; \
	done

run-big-examples: big-examples
	@for ex in $(BIG_EXAMPLES); do \
		echo ""; \
		echo "======================================"; \
		echo "Running $$ex"; \
		echo "======================================"; \
		$(QUADRA_TIME_CMD) ./$$ex; \
	done

clean-examples:
	rm -f $(EXAMPLES)

clean-big-examples:
	rm -f $(BIG_EXAMPLES)

.PHONY: test-big-laplace-convergence-contract
test-big-laplace-convergence-contract:
	./tests/run_big_laplace_convergence_contract.sh

.PHONY: test-contracts
test-contracts: test-big-laplace-convergence-contract

test-fixed-effect-covariance: tests/test_fixed_effect_covariance
	./tests/test_fixed_effect_covariance

tests/test_fixed_effect_covariance: tests/test_fixed_effect_covariance.cpp core/inference/fixed_effect_covariance.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ tests/test_fixed_effect_covariance.cpp

test-fixed-effect-report: tests/test_fixed_effect_report
	./tests/test_fixed_effect_report

tests/test_fixed_effect_report: tests/test_fixed_effect_report.cpp core/inference/fixed_effect_report.hpp core/inference/fixed_effect_covariance.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ tests/test_fixed_effect_report.cpp

test-delta-method: tests/test_delta_method
	./tests/test_delta_method

tests/test_delta_method: tests/test_delta_method.cpp core/inference/delta_method.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ tests/test_delta_method.cpp

test-big-catch-at-age-derived: tests/test_big_catch_at_age_derived
	./tests/test_big_catch_at_age_derived

tests/test_big_catch_at_age_derived: tests/test_big_catch_at_age_derived.cpp examples/big/catch_at_age_derived.hpp examples/big/catch_at_age_shared.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ tests/test_big_catch_at_age_derived.cpp

test-ad-delta-method: tests/test_ad_delta_method
	./tests/test_ad_delta_method

tests/test_ad_delta_method: tests/test_ad_delta_method.cpp core/inference/ad_delta_method.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ tests/test_ad_delta_method.cpp

test-ad-delta-method-matches-fd: tests/test_ad_delta_method_matches_fd
	./tests/test_ad_delta_method_matches_fd

tests/test_ad_delta_method_matches_fd: tests/test_ad_delta_method_matches_fd.cpp core/inference/ad_delta_method.hpp core/inference/delta_method.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ tests/test_ad_delta_method_matches_fd.cpp

test-ad-delta-method-vector: tests/test_ad_delta_method_vector
	./tests/test_ad_delta_method_vector

tests/test_ad_delta_method_vector: tests/test_ad_delta_method_vector.cpp core/inference/ad_delta_method_vector.hpp core/inference/ad_delta_method.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ tests/test_ad_delta_method_vector.cpp

test-ift-dense-validation: tests/test_ift_dense_validation
	./tests/test_ift_dense_validation

tests/test_ift_dense_validation: tests/test_ift_dense_validation.cpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ tests/test_ift_dense_validation.cpp

test-ift-vector-validation: tests/test_ift_vector_validation
	./tests/test_ift_vector_validation

tests/test_ift_vector_validation: tests/test_ift_vector_validation.cpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -o $@ tests/test_ift_vector_validation.cpp

test-quadra-ift-blocks: tests/test_quadra_ift_blocks
	./tests/test_quadra_ift_blocks

tests/test_quadra_ift_blocks: tests/test_quadra_ift_blocks.cpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ tests/test_quadra_ift_blocks.cpp

test-laplace-mode-sensitivity: tests/test_laplace_mode_sensitivity
	./tests/test_laplace_mode_sensitivity

tests/test_laplace_mode_sensitivity: tests/test_laplace_mode_sensitivity.cpp core/laplace/laplace_mode_sensitivity.hpp core/inference/ift_mode_sensitivity.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ tests/test_laplace_mode_sensitivity.cpp

test-laplace-derived-gradient: tests/test_laplace_derived_gradient
	./tests/test_laplace_derived_gradient

tests/test_laplace_derived_gradient: tests/test_laplace_derived_gradient.cpp core/laplace/laplace_derived_gradient.hpp core/laplace/laplace_mode_sensitivity.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ tests/test_laplace_derived_gradient.cpp

test-laplace-implicit-derived-gradient: tests/test_laplace_implicit_derived_gradient
	./tests/test_laplace_implicit_derived_gradient

tests/test_laplace_implicit-derived-gradient: tests/test_laplace_implicit_derived_gradient.cpp core/laplace/laplace_implicit_derivatives.hpp core/laplace/laplace_derived_gradient.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ tests/test_laplace_implicit_derived_gradient.cpp

test-laplace-profiled-derived-gradient: tests/test_laplace_profiled_derived_gradient
	./tests/test_laplace_profiled_derived_gradient

tests/test_laplace_profiled_derived_gradient: tests/test_laplace_profiled_derived_gradient.cpp core/laplace/laplace_profiled_derived_gradient.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ tests/test_laplace_profiled_derived_gradient.cpp

test-laplace-profiled-ad-gradient: tests/test_laplace_profiled_ad_gradient
	./tests/test_laplace_profiled_ad_gradient

tests/test_laplace_profiled_ad_gradient: tests/test_laplace_profiled_ad_gradient.cpp core/laplace/laplace_profiled_ad_gradient.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ tests/test_laplace_profiled_ad_gradient.cpp

test-ad-mixed-hessian-matches-fd: tests/test_ad_mixed_hessian_matches_fd
	./tests/test_ad_mixed_hessian_matches_fd

tests/test_ad_mixed_hessian_matches_fd: tests/test_ad_mixed_hessian_matches_fd.cpp examples/big/catch_at_age_shared.hpp core/laplace/laplace_implicit_derivatives.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ tests/test_ad_mixed_hessian_matches_fd.cpp

test-sparse-ldlt-factorization-cache: tests/test_sparse_ldlt_factorization_cache
	./tests/test_sparse_ldlt_factorization_cache

tests/test_sparse_ldlt_factorization_cache: tests/test_sparse_ldlt_factorization_cache.cpp core/laplace/sparse_factorization_cache.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ tests/test_sparse_ldlt_factorization_cache.cpp

test-laplace-profiled-delta-method: tests/test_laplace_profiled_delta_method
	./tests/test_laplace_profiled_delta_method

tests/test_laplace_profiled_delta_method: tests/test_laplace_profiled_delta_method.cpp core/laplace/laplace_profiled_delta_method.hpp core/laplace/laplace_profiled_derived_gradient.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ tests/test_laplace_profiled_delta_method.cpp

test-laplace-profiled-delta-method-vector: tests/test_laplace_profiled_delta_method_vector
	./tests/test_laplace_profiled_delta_method_vector

tests/test_laplace_profiled_delta_method_vector: tests/test_laplace_profiled_delta_method_vector.cpp core/laplace/laplace_profiled_delta_method_vector.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ tests/test_laplace_profiled_delta_method_vector.cpp

test-laplace-implicit-workspace: tests/test_laplace_implicit_workspace
	./tests/test_laplace_implicit_workspace

tests/test_laplace_implicit_workspace: tests/test_laplace_implicit_workspace.cpp core/laplace/laplace_implicit_workspace.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ tests/test_laplace_implicit_workspace.cpp

test-laplace-profiled-derived-report: tests/test_laplace_profiled_derived_report
	./tests/test_laplace_profiled_derived_report

tests/test_laplace_profiled_derived_report: tests/test_laplace_profiled_derived_report.cpp core/laplace/laplace_profiled_derived_report.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ tests/test_laplace_profiled_derived_report.cpp

test-report-serialization: tests/test_report_serialization
	./tests/test_report_serialization

tests/test_report_serialization: tests/test_report_serialization.cpp core/inference/report_serialization.hpp core/laplace/laplace_profiled_derived_report.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ tests/test_report_serialization.cpp

run-simple-random-intercept-example: examples/simple/random_intercept_model
	./examples/simple/random_intercept_model

examples/simple/random_intercept_model: examples/simple/random_intercept_model.cpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ examples/simple/random_intercept_model.cpp

test-public-api-headers: tests/test_public_api_headers
	./tests/test_public_api_headers

tests/test_public_api_headers: tests/test_public_api_headers.cpp include/quadra/quadra.hpp include/quadra/laplace.hpp include/quadra/inference.hpp include/quadra/io.hpp include/quadra/model.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ tests/test_public_api_headers.cpp

benchmark-random-intercept: benchmarks/random_intercept/benchmark_random_intercept
	./benchmarks/random_intercept/benchmark_random_intercept

benchmarks/random_intercept/benchmark_random_intercept: benchmarks/random_intercept/benchmark_random_intercept.cpp include/quadra/quadra.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ benchmarks/random_intercept/benchmark_random_intercept.cpp

benchmark-state-space: benchmarks/state_space/benchmark_state_space
	./benchmarks/state_space/benchmark_state_space

benchmarks/state_space/benchmark_state_space: benchmarks/state_space/benchmark_state_space.cpp include/quadra/quadra.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ benchmarks/state_space/benchmark_state_space.cpp

benchmark-quadra-tmb-random-intercept-quadra: benchmarks/comparisons/tmb_random_intercept/quadra_random_intercept_compare
	./benchmarks/comparisons/tmb_random_intercept/quadra_random_intercept_compare

benchmarks/comparisons/tmb_random_intercept/quadra_random_intercept_compare: benchmarks/comparisons/tmb_random_intercept/quadra_random_intercept_compare.cpp include/quadra/quadra.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ benchmarks/comparisons/tmb_random_intercept/quadra_random_intercept_compare.cpp

benchmark-quadra-tmb-random-intercept-tmb:
	./benchmarks/comparisons/tmb_random_intercept/run_tmb_random_intercept.R

benchmark-normalize-random-intercept:
	python3 benchmarks/analysis/normalize_random_intercept_results.py

benchmark-plot-random-intercept:
	Rscript benchmarks/analysis/plot_random_intercept_scaling.R

benchmark-normalize-rss:
	python3 benchmarks/analysis/parse_rss_logs.py

benchmark-normalize-all: benchmark-normalize-random-intercept benchmark-normalize-state-space benchmark-normalize-rss

benchmark-quadra-tmb-state-space-quadra: benchmarks/comparisons/tmb_state_space/quadra_state_space_compare
	./benchmarks/comparisons/tmb_state_space/quadra_state_space_compare

benchmarks/comparisons/tmb_state_space/quadra_state_space_compare: benchmarks/comparisons/tmb_state_space/quadra_state_space_compare.cpp include/quadra/quadra.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ benchmarks/comparisons/tmb_state_space/quadra_state_space_compare.cpp

benchmark-quadra-tmb-state-space-tmb:
	./benchmarks/comparisons/tmb_state_space/run_tmb_state_space.R

benchmark-plot-state-space-structure:
	Rscript benchmarks/analysis/plot_state_space_structure.R

benchmark-plot-random-intercept-comparison:
	Rscript benchmarks/analysis/plot_random_intercept_comparison.R

benchmark-normalize-state-space:
	python3 benchmarks/analysis/normalize_state_space_results.py

benchmark-plot-state-space-comparison:
	Rscript benchmarks/analysis/plot_state_space_comparison.R

benchmark-summarize:
	python3 benchmarks/analysis/summarize_benchmarks.py

benchmark-exact-laplace-gradient-state-space: benchmarks/exact_laplace_gradient/state_space_exact_gradient_benchmark
	./benchmarks/exact_laplace_gradient/state_space_exact_gradient_benchmark

benchmarks/exact_laplace_gradient/state_space_exact_gradient_benchmark: benchmarks/exact_laplace_gradient/state_space_exact_gradient_benchmark.cpp include/quadra/quadra.hpp core/laplace/laplace_exact_objective_gradient.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ benchmarks/exact_laplace_gradient/state_space_exact_gradient_benchmark.cpp

benchmark-exact-gradient-reuse: benchmarks/exact_laplace_gradient/exact_gradient_reuse_benchmark
	./benchmarks/exact_laplace_gradient/exact_gradient_reuse_benchmark

benchmarks/exact_laplace_gradient/exact_gradient_reuse_benchmark: benchmarks/exact_laplace_gradient/exact_gradient_reuse_benchmark.cpp include/quadra/quadra.hpp core/laplace/laplace_exact_objective_gradient.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ benchmarks/exact_laplace_gradient/exact_gradient_reuse_benchmark.cpp

benchmark-analyze-exact-gradient-reuse:
	python3 benchmarks/analysis/analyze_exact_gradient_reuse.py

benchmark-plot-exact-gradient-reuse:
	Rscript benchmarks/analysis/plot_exact_gradient_reuse.R

benchmark-plot-exact-gradient-scaling:
	Rscript benchmarks/analysis/plot_exact_gradient_scaling.R

benchmark-normalize-exact-gradient-comparison:
	python3 benchmarks/analysis/normalize_exact_gradient_comparison.py

benchmark-plot-exact-gradient-tmb-comparison:
	Rscript benchmarks/analysis/plot_exact_gradient_tmb_comparison.R

benchmark-factorization-reuse: benchmarks/exact_laplace_gradient/factorization_reuse_benchmark
	./benchmarks/exact_laplace_gradient/factorization_reuse_benchmark

benchmarks/exact_laplace_gradient/factorization_reuse_benchmark: benchmarks/exact_laplace_gradient/factorization_reuse_benchmark.cpp include/quadra/quadra.hpp core/laplace/sparse_factorization_cache.hpp
	$(CXX) $(CXXFLAGS) -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o $@ benchmarks/exact_laplace_gradient/factorization_reuse_benchmark.cpp

benchmark-analyze-factorization-reuse:
	python3 benchmarks/analysis/analyze_factorization_reuse.py

benchmark-plot-factorization-reuse:
	Rscript benchmarks/analysis/plot_factorization_reuse.R

benchmark-plot-report:
	Rscript benchmarks/analysis/build_benchmark_plot_report.R

