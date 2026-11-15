.DEFAULT_GOAL := prepare
CXX ?= c++
PYTHON ?= python3
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic
QUADRA_ROOT := ../../..
DEPFLAGS := -MMD -MP
-include $(wildcard build/*.d build/framework/*.d)
CPPFLAGS += -Iinclude -Ibuild -I$(QUADRA_ROOT)

.NOTPARALLEL:
.PHONY: prepare audit audit-biology audit-population audit-likelihood test
prepare:
	$(PYTHON) scripts/prepare_model.py
	$(PYTHON) scripts/prepare_biology.py
	$(PYTHON) scripts/prepare_population.py
	$(PYTHON) scripts/prepare_likelihood.py

build/audit_dynamics: src/audit_dynamics.cpp include/swordfish/data.hpp include/swordfish/dynamics.hpp include/swordfish/biology.hpp build/swordfish_biology_input.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEPFLAGS) $< -o $@

build/swordfish_input.hpp: scripts/prepare_model.py scripts/swordfish_inputs.py data/raw/swo2023_v004.dat data/raw/swo2023_v007.ctl data/raw/Report.sso
	$(PYTHON) scripts/prepare_model.py

audit: prepare build/audit_dynamics
	./build/audit_dynamics build/ss3_replay.txt build/dynamics_audit.csv build/biology_fitted_parameters.txt

build/test_dynamics: tests/test_dynamics.cpp include/swordfish/dynamics.hpp include/swordfish/data.hpp
	$(CXX) $(CPPFLAGS) -I$(QUADRA_ROOT)/external/eigen $(CXXFLAGS) $(DEPFLAGS) $< -o $@

build/swordfish_biology_input.hpp: scripts/prepare_biology.py scripts/prepare_model.py scripts/swordfish_inputs.py data/raw/swo2023_v004.dat data/raw/swo2023_v007.ctl data/raw/Report.sso
	$(PYTHON) scripts/prepare_biology.py

build/audit_biology: src/audit_biology.cpp include/swordfish/biology.hpp include/swordfish/data.hpp include/swordfish/dynamics.hpp build/swordfish_biology_input.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEPFLAGS) $< -o $@

build/test_biology: tests/test_biology.cpp include/swordfish/biology.hpp include/swordfish/data.hpp include/swordfish/dynamics.hpp build/swordfish_biology_input.hpp
	$(CXX) $(CPPFLAGS) -I$(QUADRA_ROOT)/external/eigen $(CXXFLAGS) $(DEPFLAGS) $< -o $@

audit-biology: prepare build/audit_biology
	./build/audit_biology build/ss3_biology_fixtures.txt build/biology_audit.csv

test: audit audit-biology audit-population audit-likelihood build/test_dynamics build/test_biology build/test_population build/test_input_factory build/test_likelihood build/test_model_gradient build/test_benchmarks build/test_projection build/test_delta_method build/test_csv build/test_control_start build/test_fixed_biology
	./build/test_dynamics
	./build/test_biology
	./build/test_population
	./build/test_input_factory
	./build/test_likelihood
	./build/test_model_gradient
	./build/test_benchmarks
	./build/test_projection
	./build/test_delta_method
	./build/test_csv
	./build/test_control_start
	./build/test_fixed_biology
	$(PYTHON) -m unittest discover -s tests -p 'test_*.py'

build/swordfish_population_input.hpp: scripts/prepare_population.py scripts/prepare_model.py scripts/swordfish_inputs.py data/raw/swo2023_v004.dat data/raw/swo2023_v007.ctl data/raw/Report.sso
	$(PYTHON) scripts/prepare_population.py

build/audit_population: src/audit_population.cpp include/swordfish/population.hpp include/swordfish/biology.hpp include/swordfish/dynamics.hpp include/swordfish/data.hpp build/swordfish_population_input.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEPFLAGS) $< -o $@

audit-population: prepare build/audit_population
	./build/audit_population build/biology_fitted_parameters.txt build/population_fitted_deviations.txt build/initialization_targets.txt build/ss3_replay.txt build/population_audit.csv build/population_trajectory.csv

build/test_population: tests/test_population.cpp include/swordfish/population.hpp include/swordfish/biology.hpp include/swordfish/dynamics.hpp include/swordfish/data.hpp build/swordfish_biology_input.hpp
	$(CXX) $(CPPFLAGS) -I$(QUADRA_ROOT)/external/eigen $(CXXFLAGS) $(DEPFLAGS) $< -o $@

build/test_input_factory: tests/test_input_factory.cpp include/swordfish/data.hpp build/swordfish_input.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEPFLAGS) $< -o $@

build/audit_likelihood: src/audit_likelihood.cpp include/swordfish/likelihood.hpp include/swordfish/population.hpp include/swordfish/biology.hpp include/swordfish/dynamics.hpp include/swordfish/data.hpp build/swordfish_input.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEPFLAGS) $< -o $@

audit-likelihood: prepare build/audit_likelihood
	./build/audit_likelihood build/biology_fitted_parameters.txt build/population_fitted_deviations.txt build/likelihood_targets.txt build

build/swordfish_assessment: src/swordfish_assessment.cpp $(QUADRA_ROOT)/core/optimizer.hpp $(QUADRA_ROOT)/core/laplace.hpp $(QUADRA_ROOT)/core/autodiff/model_gradient.hpp include/swordfish/model.hpp include/swordfish/reference_io.hpp include/swordfish/report.hpp include/swordfish/likelihood.hpp include/swordfish/population.hpp include/swordfish/biology.hpp include/swordfish/dynamics.hpp include/swordfish/data.hpp build/swordfish_input.hpp
	$(CXX) $(CPPFLAGS) -I$(QUADRA_ROOT)/external/eigen $(CXXFLAGS) $(DEPFLAGS) $< -o $@

.PHONY: evaluate-objective objective-gradient fit-swordfish
evaluate-objective: prepare build/swordfish_assessment
	./build/swordfish_assessment --evaluate
objective-gradient: prepare build/swordfish_assessment
	./build/swordfish_assessment --gradient
fit-swordfish: prepare build/swordfish_assessment
	./build/swordfish_assessment --fit

build/test_likelihood: tests/test_likelihood.cpp include/swordfish/likelihood.hpp include/swordfish/population.hpp include/swordfish/biology.hpp include/swordfish/dynamics.hpp include/swordfish/data.hpp build/swordfish_biology_input.hpp
	$(CXX) $(CPPFLAGS) -I$(QUADRA_ROOT)/external/eigen $(CXXFLAGS) $(DEPFLAGS) $< -o $@

build/test_model_gradient: tests/test_model_gradient.cpp include/swordfish/model.hpp include/swordfish/reference_io.hpp include/swordfish/likelihood.hpp include/swordfish/population.hpp include/swordfish/biology.hpp include/swordfish/dynamics.hpp include/swordfish/data.hpp $(QUADRA_ROOT)/core/autodiff/model_gradient.hpp build/swordfish_input.hpp
	$(CXX) $(CPPFLAGS) -I$(QUADRA_ROOT)/external/eigen $(CXXFLAGS) $(DEPFLAGS) $< -o $@

.PHONY: fitted-report
fitted-report: build/swordfish_assessment
	./build/swordfish_assessment --report

build/test_dynamics build/test_biology build/test_population build/test_likelihood: $(QUADRA_ROOT)/core/autodiff/model_gradient.hpp

.PHONY: compare-fit
compare-fit: fitted-report
	$(PYTHON) scripts/compare_fit.py

build/swordfish_diagnostics: src/swordfish_diagnostics.cpp $(wildcard include/swordfish/*.hpp) $(QUADRA_ROOT)/core/optimizer.hpp $(QUADRA_ROOT)/core/autodiff/model_gradient.hpp build/swordfish_input.hpp
	$(CXX) $(CPPFLAGS) -I$(QUADRA_ROOT)/external/eigen $(CXXFLAGS) $(DEPFLAGS) $< -o $@
.PHONY: robustness curvature
robustness: build/swordfish_diagnostics build/coordinate_covariance.csv verify-covariance
	./build/swordfish_diagnostics starts
	$(PYTHON) scripts/compare_robustness.py
curvature: build/swordfish_diagnostics
	./build/swordfish_diagnostics curvature

build/summarize_curvature: src/summarize_curvature.cpp include/swordfish/csv.hpp
	$(CXX) $(CPPFLAGS) -I$(QUADRA_ROOT)/external/eigen $(CXXFLAGS) $(DEPFLAGS) $< -o $@
.PHONY: uncertainty benchmarks
uncertainty: build/summarize_curvature build/swordfish_diagnostics
	$(PYTHON) scripts/runtime.py begin
	./build/swordfish_diagnostics curvature
	./build/swordfish_diagnostics derived 0.0001
	./build/summarize_curvature
	cp build/derived_uncertainty.csv build/derived_uncertainty_step1e-4.csv
	./build/swordfish_diagnostics derived 0.00005
	./build/summarize_curvature
	$(PYTHON) scripts/compare_uncertainty.py
	$(PYTHON) scripts/runtime.py finish
benchmarks: build/swordfish_diagnostics
	./build/swordfish_diagnostics benchmarks
	$(PYTHON) scripts/compare_benchmarks.py

build/test_benchmarks: tests/test_benchmarks.cpp $(wildcard include/swordfish/*.hpp) build/swordfish_input.hpp
	$(CXX) $(CPPFLAGS) -I$(QUADRA_ROOT)/external/eigen $(CXXFLAGS) $(DEPFLAGS) $< -o $@
.PHONY: test-benchmarks
test-benchmarks: build/test_benchmarks benchmarks
	./build/test_benchmarks

build/coordinate_covariance.csv: build/quadra_fitted_parameters.csv
	$(MAKE) uncertainty

build/swordfish_projection: src/swordfish_projection.cpp $(wildcard include/swordfish/*.hpp) build/swordfish_input.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEPFLAGS) $< -o $@
.PHONY: project-swordfish
project-swordfish: build/swordfish_projection
	./build/swordfish_projection
	$(PYTHON) scripts/compare_projections.py

build/test_projection: tests/test_projection.cpp $(wildcard include/swordfish/*.hpp) build/swordfish_input.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEPFLAGS) $< -o $@
.PHONY: test-projection
test-projection: build/test_projection project-swordfish
	./build/test_projection

.PHONY: plot-projections
plot-projections: project-swordfish
	Rscript scripts/plot_projections.R

build/swordfish_projection_uncertainty: src/swordfish_projection_uncertainty.cpp $(wildcard include/swordfish/*.hpp) build/swordfish_input.hpp
	$(CXX) $(CPPFLAGS) -I$(QUADRA_ROOT)/external/eigen $(CXXFLAGS) $(DEPFLAGS) $< -o $@
.PHONY: projection-uncertainty
projection-uncertainty: build/swordfish_projection_uncertainty build/coordinate_covariance.csv verify-covariance project-swordfish
	./build/swordfish_projection_uncertainty
	$(PYTHON) scripts/check_projection_uncertainty.py

build/test_delta_method: tests/test_delta_method.cpp include/swordfish/delta_method.hpp
	$(CXX) $(CPPFLAGS) -I$(QUADRA_ROOT)/external/eigen $(CXXFLAGS) $(DEPFLAGS) $< -o $@
.PHONY: test-projection-uncertainty
test-projection-uncertainty: build/test_delta_method projection-uncertainty
	./build/test_delta_method

.PHONY: plot-projection-uncertainty
plot-projection-uncertainty:
	$(PYTHON) scripts/check_projection_uncertainty.py
	Rscript scripts/plot_projection_uncertainty.R

.PHONY: verify-covariance test-framework release
verify-covariance:
	$(PYTHON) scripts/runtime.py verify

FRAMEWORK_TESTS := test_third_order_comparisons test_fixed_only_optimizer test_scalar_generic_model test_parameter_transforms test_stats_laplace_optimizer
test-framework:
	@mkdir -p build/framework
	@set -e; for name in $(FRAMEWORK_TESTS); do $(CXX) -std=c++17 -O2 -I$(QUADRA_ROOT)/external/eigen $(QUADRA_ROOT)/tests/$$name.cpp -o build/framework/$$name; ./build/framework/$$name; done

build/test_csv: tests/test_csv.cpp include/swordfish/csv.hpp include/swordfish/scenarios.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEPFLAGS) $< -o $@

release:
	$(PYTHON) scripts/run_release.py

CPP_PROGRAMS := $(patsubst src/%.cpp,build/%,$(wildcard src/*.cpp)) $(patsubst tests/%.cpp,build/%,$(wildcard tests/*.cpp))
$(CPP_PROGRAMS): | build
build:
	mkdir -p build

.PHONY: test-portable
test-portable: build/test_dynamics build/test_delta_method build/test_csv
	./build/test_dynamics
	./build/test_delta_method
	./build/test_csv
	$(PYTHON) -m unittest discover -s tests -p 'test_runtime.py'

.PHONY: build-all
build-all: prepare $(CPP_PROGRAMS)

build/swordfish_control_start: src/swordfish_control_start.cpp $(wildcard include/swordfish/*.hpp) build/swordfish_control_input.hpp
	$(CXX) $(CPPFLAGS) -I$(QUADRA_ROOT)/external/eigen $(CXXFLAGS) $(DEPFLAGS) $< -o $@

build/test_control_start: tests/test_control_start.cpp $(wildcard include/swordfish/*.hpp) build/swordfish_control_input.hpp
	$(CXX) $(CPPFLAGS) -I$(QUADRA_ROOT)/external/eigen $(CXXFLAGS) $(DEPFLAGS) $< -o $@
.PHONY: test-control-start
test-control-start: build/test_control_start
	./build/test_control_start

build/swordfish_control_input.hpp: scripts/prepare_control_start.py scripts/prepare_model.py scripts/swordfish_inputs.py data/raw/swo2023_v004.dat data/raw/swo2023_v007.ctl | build
	$(PYTHON) scripts/prepare_control_start.py

.PHONY: fit-control-start
fit-control-start: build/swordfish_control_start
	$(PYTHON) scripts/run_control_start.py

SS_PAR ?= data/raw/ss.par
.PHONY: fit-ss-par
fit-ss-par:
	$(PYTHON) scripts/run_ss_par.py "$(SS_PAR)"

.PHONY: native-swordfish
native-swordfish: build/swordfish_native

build/swordfish_native: src/swordfish_native.cpp $(wildcard include/swordfish/*.hpp) | build
	$(CXX) $(CPPFLAGS) -I$(QUADRA_ROOT)/external/eigen $(CXXFLAGS) $(DEPFLAGS) $< -o $@

.PHONY: test-fixed-biology
test-fixed-biology: build/test_fixed_biology
	./build/test_fixed_biology

build/test_fixed_biology: tests/test_fixed_biology.cpp $(wildcard include/swordfish/*.hpp) build/swordfish_input.hpp
	$(CXX) $(CPPFLAGS) -I$(QUADRA_ROOT)/external/eigen $(CXXFLAGS) $(DEPFLAGS) $< -o $@

.PHONY: test-optimizer-scaling
test-optimizer-scaling: build/test_optimizer_scaling
	./build/test_optimizer_scaling

build/test_optimizer_scaling: tests/test_optimizer_scaling.cpp include/swordfish/optimizer_scaling.hpp include/swordfish/delta_method.hpp include/swordfish/csv.hpp | build
	$(CXX) $(CPPFLAGS) -I$(QUADRA_ROOT)/external/eigen $(CXXFLAGS) $(DEPFLAGS) $< -o $@

NATIVE_INPUT ?= build/native/swordfish-ss-par-start.json
PDF_PYTHON ?= $(PYTHON)
.PHONY: native-release
native-release:
	$(PYTHON) scripts/run_native_release.py --input "$(NATIVE_INPUT)" --pdf-python "$(PDF_PYTHON)"
