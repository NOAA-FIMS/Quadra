# Shared include paths for local builds.
#
# Keep dependency include paths centralized so tests, examples, and ad hoc
# compile commands do not drift from one another.

QUADRA_ROOT ?= .

QUADRA_INCLUDE_FLAGS := \
	-I$(QUADRA_ROOT) \
	-I$(QUADRA_ROOT)/external/eigen \
	-I$(QUADRA_ROOT)/external/had \
	-I$(QUADRA_ROOT)/external/LBFGSpp/include

QUADRA_CXX ?= c++
QUADRA_CXXFLAGS ?= -std=c++17 -O3 -flto

.PHONY: check-external-deps
check-external-deps:
	@test -f "$(QUADRA_ROOT)/external/eigen/Eigen/Core" || \
		(echo "Missing Eigen: expected $(QUADRA_ROOT)/external/eigen/Eigen/Core"; exit 1)
