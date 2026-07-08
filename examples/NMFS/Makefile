# Standalone NMFS examples Makefile
#
# Intended location:
#   examples/NMFS/Makefile
#
# Usage from repository root:
#   make -C examples/NMFS help
#   make -C examples/NMFS red_snapper
#   make -C examples/NMFS opakapaka
#   make -C examples/NMFS all
#
# Usage from examples/NMFS:
#   make help
#   make red_snapper
#   make opakapaka
#
# Notes:
# - This Makefile intentionally does not depend on the top-level Makefile.
# - It compiles the Quadra examples directly with clang++/g++.
# - If an example runner script exists, targets prefer that script.
# - If no known runner script/source is found, the target prints searched paths.

SHELL := /usr/bin/env bash

CXX ?= clang++
STD ?= c++17
OPT ?= -O3
WARN ?= -Wall -Wextra
LTO ?=

# Repo root relative to examples/NMFS.
REPO_ROOT := ../..

# Common include paths used by Quadra examples.
INCLUDES := \
	-I$(REPO_ROOT) \
	-I$(REPO_ROOT)/core \
	-I$(REPO_ROOT)/external/eigen \
	-I$(REPO_ROOT)/external \
	-I$(REPO_ROOT)/external/LBFGSpp/include

CXXFLAGS ?= -std=$(STD) $(OPT) $(WARN) $(LTO) $(INCLUDES)

BIN_DIR := build

.PHONY: help all red_snapper opakapaka clean print-config

help:
	@echo "Standalone NMFS examples Makefile"
	@echo
	@echo "Targets:"
	@echo "  make red_snapper   Build/run the SEFSC red snapper Quadra example"
	@echo "  make opakapaka     Build/run the PIFSC opakapaka Quadra example"
	@echo "  make all           Build/run both examples"
	@echo "  make clean         Remove local examples/NMFS/build outputs"
	@echo "  make print-config  Show compiler and include settings"
	@echo
	@echo "Overrides:"
	@echo "  CXX=clang++|g++ OPT=-O0|-O2|-O3 LTO=-flto"

print-config:
	@echo "CXX      = $(CXX)"
	@echo "CXXFLAGS = $(CXXFLAGS)"
	@echo "REPO_ROOT= $(REPO_ROOT)"
	@echo "BIN_DIR  = $(BIN_DIR)"

all: red_snapper opakapaka

$(BIN_DIR):
	@mkdir -p "$(BIN_DIR)"

red_snapper: $(BIN_DIR)
	@set -euo pipefail; \
	echo "== NMFS red snapper =="; \
	if [[ -x "$(REPO_ROOT)/run_sefsc_red_snapper.sh" ]]; then \
	  echo "using $(REPO_ROOT)/run_sefsc_red_snapper.sh"; \
	  cd "$(REPO_ROOT)" && ./run_sefsc_red_snapper.sh; \
	elif [[ -x "$(REPO_ROOT)/run_red_snapper.sh" ]]; then \
	  echo "using $(REPO_ROOT)/run_red_snapper.sh"; \
	  cd "$(REPO_ROOT)" && ./run_red_snapper.sh; \
	elif [[ -x "$(REPO_ROOT)/examples/sefsc_red_snapper/run_sefsc_red_snapper.sh" ]]; then \
	  echo "using $(REPO_ROOT)/examples/sefsc_red_snapper/run_sefsc_red_snapper.sh"; \
	  cd "$(REPO_ROOT)" && ./examples/sefsc_red_snapper/run_sefsc_red_snapper.sh; \
	elif [[ -f "$(REPO_ROOT)/examples/sefsc_red_snapper/quadra/sefsc_red_snapper.cpp" ]]; then \
	  src="$(REPO_ROOT)/examples/sefsc_red_snapper/quadra/sefsc_red_snapper.cpp"; \
	  bin="$(BIN_DIR)/sefsc_red_snapper"; \
	  echo "building $$src -> $$bin"; \
	  $(CXX) $(CXXFLAGS) "$$src" -o "$$bin"; \
	  echo "running $$bin"; \
	  "$$bin"; \
	elif [[ -f "$(REPO_ROOT)/examples/SEFSC/red_snapper/quadra/red_snapper.cpp" ]]; then \
	  src="$(REPO_ROOT)/examples/SEFSC/red_snapper/quadra/red_snapper.cpp"; \
	  bin="$(BIN_DIR)/sefsc_red_snapper"; \
	  echo "building $$src -> $$bin"; \
	  $(CXX) $(CXXFLAGS) "$$src" -o "$$bin"; \
	  echo "running $$bin"; \
	  "$$bin"; \
	else \
	  echo "ERROR: could not find a known red snapper runner or source."; \
	  echo "Searched:"; \
	  echo "  $(REPO_ROOT)/run_sefsc_red_snapper.sh"; \
	  echo "  $(REPO_ROOT)/run_red_snapper.sh"; \
	  echo "  $(REPO_ROOT)/examples/sefsc_red_snapper/run_sefsc_red_snapper.sh"; \
	  echo "  $(REPO_ROOT)/examples/sefsc_red_snapper/quadra/sefsc_red_snapper.cpp"; \
	  echo "  $(REPO_ROOT)/examples/SEFSC/red_snapper/quadra/red_snapper.cpp"; \
	  exit 1; \
	fi

opakapaka: $(BIN_DIR)
	@set -euo pipefail; \
	echo "== NMFS opakapaka =="; \
	if [[ -x "$(REPO_ROOT)/run_opakapaka.sh" ]]; then \
	  echo "using $(REPO_ROOT)/run_opakapaka.sh"; \
	  cd "$(REPO_ROOT)" && ./run_opakapaka.sh; \
	elif [[ -x "$(REPO_ROOT)/examples/NMFS/pifsc_opakapaka/run_opakapaka.sh" ]]; then \
	  echo "using $(REPO_ROOT)/examples/NMFS/pifsc_opakapaka/run_opakapaka.sh"; \
	  cd "$(REPO_ROOT)" && ./examples/NMFS/pifsc_opakapaka/run_opakapaka.sh; \
	elif [[ -x "$(REPO_ROOT)/examples/NMFS/pifsc_opakapaka/quadra/run_opakapaka.sh" ]]; then \
	  echo "using $(REPO_ROOT)/examples/NMFS/pifsc_opakapaka/quadra/run_opakapaka.sh"; \
	  cd "$(REPO_ROOT)" && ./examples/NMFS/pifsc_opakapaka/quadra/run_opakapaka.sh; \
	elif [[ -f "$(REPO_ROOT)/examples/NMFS/pifsc_opakapaka/quadra/opakapaka.cpp" ]]; then \
	  src="$(REPO_ROOT)/examples/NMFS/pifsc_opakapaka/quadra/opakapaka.cpp"; \
	  bin="$(BIN_DIR)/pifsc_opakapaka"; \
	  echo "building $$src -> $$bin"; \
	  $(CXX) $(CXXFLAGS) "$$src" -o "$$bin"; \
	  echo "running $$bin"; \
	  "$$bin"; \
	elif [[ -f "$(REPO_ROOT)/examples/NMFS/pifsc_opakapaka/quadra/pifsc_opakapaka.cpp" ]]; then \
	  src="$(REPO_ROOT)/examples/NMFS/pifsc_opakapaka/quadra/pifsc_opakapaka.cpp"; \
	  bin="$(BIN_DIR)/pifsc_opakapaka"; \
	  echo "building $$src -> $$bin"; \
	  $(CXX) $(CXXFLAGS) "$$src" -o "$$bin"; \
	  echo "running $$bin"; \
	  "$$bin"; \
	else \
	  echo "ERROR: could not find a known opakapaka runner or source."; \
	  echo "Searched:"; \
	  echo "  $(REPO_ROOT)/run_opakapaka.sh"; \
	  echo "  $(REPO_ROOT)/examples/NMFS/pifsc_opakapaka/run_opakapaka.sh"; \
	  echo "  $(REPO_ROOT)/examples/NMFS/pifsc_opakapaka/quadra/run_opakapaka.sh"; \
	  echo "  $(REPO_ROOT)/examples/NMFS/pifsc_opakapaka/quadra/opakapaka.cpp"; \
	  echo "  $(REPO_ROOT)/examples/NMFS/pifsc_opakapaka/quadra/pifsc_opakapaka.cpp"; \
	  exit 1; \
	fi

clean:
	@rm -rf "$(BIN_DIR)"
	@echo "removed $(BIN_DIR)"
