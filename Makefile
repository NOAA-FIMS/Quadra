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
