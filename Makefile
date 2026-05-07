CXX ?= clang++
CXXFLAGS ?= -std=c++17 -O3 -flto -Icore/eigen

CORE_SRC = core/autodiff/adgraph.cpp

TESTS = \
	tests/test_curvature_depends_on_theta \
	tests/test_poisson_random_effect \
	tests/test_ar1_random_walk \
	tests/test_hdot_validation \
	tests/test_fixed_covariance \
	tests/test_report_registry \
	tests/test_report_v2 
	

EXAMPLES = \
	examples/fisheries_random_year_effects \
	examples/fisheries_age_selectivity_random_walk \
	examples/fisheries_index_cpue_laplace \
	examples/fisheries_report_registry \
	examples/fisheries_report_v2 
	
	

all: $(TESTS) $(EXAMPLES)

tests/%: tests/%.cpp $(CORE_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $< $(CORE_SRC)

examples/%: examples/%.cpp $(CORE_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $< $(CORE_SRC)

validate-hdot:
	$(CXX) $(CXXFLAGS) -DQUADRA_VALIDATE_HDOT -o tests/test_hdot_validation tests/test_hdot_validation.cpp $(CORE_SRC)
	./tests/test_hdot_validation

run-tests: $(TESTS)
	./tests/test_curvature_depends_on_theta
	./tests/test_poisson_random_effect
	./tests/test_ar1_random_walk
	./tests/test_hdot_validation
	./tests/test_fixed_covariance
	./tests/test_report_registry

clean:
	rm -f $(TESTS) $(EXAMPLES)

.PHONY: all run-tests validate-hdot clean
