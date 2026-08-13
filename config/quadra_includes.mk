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

# Rtools installs g++.exe and its assembler in the same target bin directory.
# On Windows, selecting that compiler without also exporting its directory on
# PATH lets GNU Make find g++ but prevents g++ from launching as.exe.
ifeq ($(OS),Windows_NT)
QUADRA_WINDOWS_CXX_CANDIDATES := \
	C:/rtools45/x86_64-w64-mingw32.static.posix/bin/g++.exe \
	C:/rtools44/x86_64-w64-mingw32.static.posix/bin/g++.exe \
	C:/rtools43/x86_64-w64-mingw32.static.posix/bin/g++.exe \
	C:/msys64/ucrt64/bin/g++.exe \
	C:/msys64/mingw64/bin/g++.exe
QUADRA_DETECTED_WINDOWS_CXX := \
	$(firstword $(foreach compiler,$(QUADRA_WINDOWS_CXX_CANDIDATES),\
		$(if $(wildcard $(compiler)),$(compiler))))
ifneq ($(strip $(QUADRA_DETECTED_WINDOWS_CXX)),)
QUADRA_CXX ?= $(QUADRA_DETECTED_WINDOWS_CXX)
QUADRA_TOOLCHAIN_BIN := $(patsubst %/,%,$(dir $(QUADRA_DETECTED_WINDOWS_CXX)))
export PATH := $(QUADRA_TOOLCHAIN_BIN):$(PATH)
endif
endif

QUADRA_CXX ?= c++
QUADRA_CXXFLAGS ?= -std=c++17 -O3 -flto

.PHONY: check-external-deps
check-external-deps:
	@test -f "$(QUADRA_ROOT)/external/eigen/Eigen/Core" || \
		(echo "Missing Eigen: expected $(QUADRA_ROOT)/external/eigen/Eigen/Core"; exit 1)
