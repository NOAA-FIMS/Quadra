# Quadra top-level Makefile
#
# Intended location:
#   ./Makefile
#
# Usage:
#   make help
#   make nmfs
#   make red_snapper
#   make opakapaka
#   make pollock

SHELL := /usr/bin/env bash

CXX ?= clang++
STD ?= c++17
OPT ?= -O3
WARN ?= -Wall -Wextra
LTO ?=

CXXFLAGS ?= -std=$(STD) $(OPT) $(WARN) $(LTO) \
	-I. \
	-Icore \
	-Iexternal/eigen \
	-Iexternal \
	-Iexternal/LBFGSpp/include

BUILD_DIR ?= build/examples

.PHONY: help all nmfs red_snapper opakapaka pollock clean print-config

help:
	@echo "Quadra top-level Makefile"
	@echo
	@echo "Targets:"
	@echo "  make nmfs         Run maintained NMFS examples"
	@echo "  make red_snapper  Run SEFSC red snapper Quadra fit"
	@echo "  make opakapaka    Build/run PIFSC opakapaka"
	@echo "  make pollock      Run AFSC walleye pollock"
	@echo "  make all          Alias for nmfs"
	@echo "  make clean        Remove local example build outputs"
	@echo "  make print-config Show compiler settings"

all: nmfs

nmfs: red_snapper opakapaka pollock

print-config:
	@echo "CXX       = $(CXX)"
	@echo "CXXFLAGS  = $(CXXFLAGS)"
	@echo "BUILD_DIR = $(BUILD_DIR)"

$(BUILD_DIR):
	@mkdir -p "$(BUILD_DIR)"

red_snapper:
	@./examples/NMFS/sefsc_red_snapper/run_red_snapper_quadra_fit.sh

opakapaka: $(BUILD_DIR)
	@echo "== Building PIFSC opakapaka =="
	@$(CXX) $(CXXFLAGS) \
		examples/NMFS/pifsc_opakapaka/quadra/opakapaka.cpp \
		examples/NMFS/pifsc_opakapaka/quadra/opakapaka_adgraph_global.cpp \
		-o "$(BUILD_DIR)/opakapaka"
	@echo "== Running PIFSC opakapaka =="
	@"$(BUILD_DIR)/opakapaka"

pollock:
	@./examples/NMFS/afsc_walleye_pollock/run_walleye_pollock_example.sh

clean:
	@rm -rf "$(BUILD_DIR)"
	@echo "removed $(BUILD_DIR)"
