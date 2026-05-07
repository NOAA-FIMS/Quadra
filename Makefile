
CXX ?= clang++
CXXFLAGS ?= -std=c++17 -O3 -flto -Icore/eigen

TESTS = \
	tests/test_curvature_depends_on_theta \
	tests/test_poisson_random_effect \
	tests/test_ar1_random_walk \
	tests/test_hdot_validation

EXAMPLES = \
	examples/fisheries_random_year_effects \
	examples/fisheries_age_selectivity_random_walk \
	examples/fisheries_index_cpue_laplace

all: $(TESTS) $(EXAMPLES)

tests/%: tests/%.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

examples/%: examples/%.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

validate-hdot:
	$(CXX) $(CXXFLAGS) -DQUADRA_VALIDATE_HDOT -o tests/test_hdot_validation tests/test_hdot_validation.cpp
	./tests/test_hdot_validation

run-tests: $(TESTS)
	./tests/test_curvature_depends_on_theta
	./tests/test_poisson_random_effect
	./tests/test_ar1_random_walk
	./tests/test_hdot_validation

clean:
	rm -f $(TESTS) $(EXAMPLES)

.PHONY: all run-tests validate-hdot clean
