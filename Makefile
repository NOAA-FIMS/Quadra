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

examples/%: examples/%.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

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

run-tests: $(TESTS)
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
