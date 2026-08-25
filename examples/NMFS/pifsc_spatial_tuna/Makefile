CXX := c++
QUADRA_ROOT ?= ../../..
BUILD_DIR ?= build
CONFIG ?= config/tuna_assessment.conf
HDOT_WORKERS ?= 0
COMPREHENSIVE_MULTISTART ?= 4
RUN_OUTPUT_DIR ?= $(BUILD_DIR)/assessment_outputs/data
REPORT_OUTPUT_DIR ?= $(BUILD_DIR)/assessment_outputs/report
RUN_NUTS ?= 0
RSCRIPT ?= Rscript

CPPFLAGS += -Iinclude -I$(QUADRA_ROOT) -I$(QUADRA_ROOT)/external/eigen
CXXFLAGS ?= -O2 -Wall -Wextra
CXXFLAGS += -std=c++17
EXAMPLE_SRC := src/advanced_tuna_spatial_assessment_example.cpp
EXAMPLE_BIN := $(BUILD_DIR)/advanced_tuna_spatial_assessment_example
EXAMPLE_DEP := $(EXAMPLE_BIN).d
NATIVE_FLOW_SRC := src/train_dependency_free_transport_flow.cpp
NATIVE_FLOW_BIN := $(BUILD_DIR)/train_native_transport_flow
NATIVE_FLOW_DEP := $(NATIVE_FLOW_BIN).d
WORKFLOW_SRC := src/run_full_workflow.cpp
WORKFLOW_BIN := $(BUILD_DIR)/run_tuna_full_workflow
WORKFLOW_DEP := $(WORKFLOW_BIN).d
SMOKE_SRC := tests/catch_conditioning_smoke.cpp
SMOKE_BIN := $(BUILD_DIR)/catch_conditioning_smoke
SMOKE_DEP := $(SMOKE_BIN).d
RECOVERY_STAMP := $(BUILD_DIR)/.simulation_recovery_passed

.PHONY: all advanced-tuna-example full-tuna-workflow run-full-tuna-workflow \
	run-advanced-tuna-example \
	fit-advanced-tuna sample-advanced-tuna run-advanced-tuna-nuts \
	run-comprehensive-analysis test test-fast \
	test-simulation-recovery force-simulation-recovery \
	check-simulation-recovery sensitivity-grid audit-opal-data \
	check-posterior-assessment assessment-report train-transport-flow \
	spatial-pulse-map \
	train-native-transport-flow \
	check-transport-gradient check-transport-reproducibility \
	test-catch-conditioning clean

all: advanced-tuna-example full-tuna-workflow

advanced-tuna-example: $(EXAMPLE_BIN)

full-tuna-workflow: $(EXAMPLE_BIN) $(WORKFLOW_BIN)

$(WORKFLOW_BIN): $(WORKFLOW_SRC) | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP $< -o $@ \
		$(LDFLAGS) $(LDLIBS)

run-full-tuna-workflow: $(EXAMPLE_BIN) $(WORKFLOW_BIN)
	$(WORKFLOW_BIN) --config $(CONFIG) --output-dir $(RUN_OUTPUT_DIR) \
		--report-dir $(REPORT_OUTPUT_DIR) --rscript "$(RSCRIPT)"

$(BUILD_DIR):
	mkdir -p $@

$(EXAMPLE_BIN): $(EXAMPLE_SRC) | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP $< -o $@ \
		$(LDFLAGS) $(LDLIBS)

run-advanced-tuna-example: $(EXAMPLE_BIN)
	QUADRA_TUNA_BASELINE_ONLY=1 \
		QUADRA_TUNA_OUTPUT_DIR=$(RUN_OUTPUT_DIR) \
		QUADRA_TUNA_RUN_NUTS=$(RUN_NUTS) \
		$< --config $(CONFIG)
	"$(RSCRIPT)" scripts/build_assessment_report.R $(RUN_OUTPUT_DIR) $(REPORT_OUTPUT_DIR)

fit-advanced-tuna: RUN_NUTS=0
fit-advanced-tuna: run-advanced-tuna-example

sample-advanced-tuna: $(EXAMPLE_BIN)
	QUADRA_TUNA_BASELINE_ONLY=1 \
		QUADRA_TUNA_LOAD_FIT_CHECKPOINT=1 \
		QUADRA_TUNA_OUTPUT_DIR=$(RUN_OUTPUT_DIR) \
		QUADRA_TUNA_RUN_NUTS=1 \
		$< --config $(CONFIG)
	"$(RSCRIPT)" tests/check_posterior_assessment.R $(RUN_OUTPUT_DIR) $(REPORT_OUTPUT_DIR)
	"$(RSCRIPT)" scripts/build_assessment_report.R $(RUN_OUTPUT_DIR) $(REPORT_OUTPUT_DIR)

run-advanced-tuna-nuts: sample-advanced-tuna

# Explicit opt-in to the expensive composition grid, retrospective peels,
# sensitivity scenarios, and simulation-estimation loop.
run-comprehensive-analysis: $(EXAMPLE_BIN)
	QUADRA_TUNA_MULTISTART=$(COMPREHENSIVE_MULTISTART) $< --config $(CONFIG)

# `test` is cached after a successful full recovery run. Use
# `force-simulation-recovery` when an unconditional rerun is desired.
test: test-fast test-simulation-recovery

test-fast: test-catch-conditioning

test-simulation-recovery: $(RECOVERY_STAMP)

$(RECOVERY_STAMP): $(EXAMPLE_BIN) tests/check_simulation_recovery.sh | $(BUILD_DIR)
	QUADRA_TUNA_BASELINE_ONLY=1 QUADRA_TUNA_MULTISTART=1 \
		QUADRA_TUNA_MAX_PHASE_ITERATIONS=100 \
		QUADRA_TUNA_OUTPUT_DIR=$(RUN_OUTPUT_DIR) \
		QUADRA_TUNA_LBFGS_PRINT_EVERY=0 $(EXAMPLE_BIN) --config $(CONFIG)
	bash tests/check_simulation_recovery.sh $(RUN_OUTPUT_DIR)
	touch $@

force-simulation-recovery: $(EXAMPLE_BIN)
	QUADRA_TUNA_BASELINE_ONLY=1 QUADRA_TUNA_MULTISTART=1 \
		QUADRA_TUNA_MAX_PHASE_ITERATIONS=100 \
		QUADRA_TUNA_OUTPUT_DIR=$(RUN_OUTPUT_DIR) \
		QUADRA_TUNA_LBFGS_PRINT_EVERY=0 $(EXAMPLE_BIN) --config $(CONFIG)
	bash tests/check_simulation_recovery.sh $(RUN_OUTPUT_DIR)
	touch $(RECOVERY_STAMP)

check-simulation-recovery:
	bash tests/check_simulation_recovery.sh $(RUN_OUTPUT_DIR)

check-posterior-assessment:
	"$(RSCRIPT)" tests/check_posterior_assessment.R $(RUN_OUTPUT_DIR) $(REPORT_OUTPUT_DIR)

check-transport-gradient: $(NATIVE_FLOW_BIN)
	$(NATIVE_FLOW_BIN) \
		--draws $(BUILD_DIR)/assessment_outputs/data/posterior_draws.csv \
		--output $(BUILD_DIR)/transport_flow/gradient_check_smoke.qflow \
		--geometry-cache $(RUN_OUTPUT_DIR)/marginal_whitening_cache.csv \
		--fit-metadata $(RUN_OUTPUT_DIR)/fit_checkpoint_metadata.csv \
		--layers 3 --hidden 8 --epochs 2 --patience 2 --seed 20260823

check-transport-reproducibility: $(NATIVE_FLOW_BIN)
	$(NATIVE_FLOW_BIN) --draws \
		$(BUILD_DIR)/assessment_outputs/data/posterior_draws.csv \
		--output $(BUILD_DIR)/transport_flow/reproducibility_a.qflow \
		--geometry-cache $(RUN_OUTPUT_DIR)/marginal_whitening_cache.csv \
		--fit-metadata $(RUN_OUTPUT_DIR)/fit_checkpoint_metadata.csv \
		--layers 3 --hidden 8 --epochs 2 --patience 2 --seed 20260823
	$(NATIVE_FLOW_BIN) --draws \
		$(BUILD_DIR)/assessment_outputs/data/posterior_draws.csv \
		--output $(BUILD_DIR)/transport_flow/reproducibility_b.qflow \
		--geometry-cache $(RUN_OUTPUT_DIR)/marginal_whitening_cache.csv \
		--fit-metadata $(RUN_OUTPUT_DIR)/fit_checkpoint_metadata.csv \
		--layers 3 --hidden 8 --epochs 2 --patience 2 --seed 20260823
	cmp $(BUILD_DIR)/transport_flow/reproducibility_a.qflow \
		$(BUILD_DIR)/transport_flow/reproducibility_b.qflow

sensitivity-grid: $(EXAMPLE_BIN)
	QUADRA_TUNA_HDOT_WORKERS=$(HDOT_WORKERS) \
		bash tests/run_sensitivity_grid.sh $(EXAMPLE_BIN) $(BUILD_DIR)/sensitivity_grid

audit-opal-data: | $(BUILD_DIR)
	"$(RSCRIPT)" scripts/audit_opal_wcpo_bet.R data/opal_raw $(BUILD_DIR)/opal_data_audit.csv

assessment-report:
	"$(RSCRIPT)" scripts/build_assessment_report.R $(RUN_OUTPUT_DIR) $(REPORT_OUTPUT_DIR)

spatial-pulse-map:
	"$(RSCRIPT)" scripts/build_spatial_pulse_map.R $(RUN_OUTPUT_DIR) $(REPORT_OUTPUT_DIR)

train-transport-flow: train-native-transport-flow

train-native-transport-flow: $(NATIVE_FLOW_BIN)
	@for seed in 20260823 20260824 20260825 20260826; do \
		$(NATIVE_FLOW_BIN) \
			--draws $(BUILD_DIR)/assessment_outputs/data/posterior_draws.csv \
			--output $(BUILD_DIR)/transport_flow/native_realnvp_$${seed}.qflow \
			--geometry-cache $(RUN_OUTPUT_DIR)/marginal_whitening_cache.csv \
			--fit-metadata $(RUN_OUTPUT_DIR)/fit_checkpoint_metadata.csv \
			--layers 8 --hidden 64 --epochs 1000 --patience 120 \
			--seed $${seed}; \
	done

$(NATIVE_FLOW_BIN): $(NATIVE_FLOW_SRC) include/tuna/dependency_free_transport_flow.hpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP $< -o $@ \
		$(LDFLAGS) $(LDLIBS)

test-catch-conditioning: $(SMOKE_BIN)
	$<

$(SMOKE_BIN): $(SMOKE_SRC) | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP $< -o $@ $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(EXAMPLE_BIN) $(EXAMPLE_DEP) $(SMOKE_BIN) $(SMOKE_DEP) \
		$(NATIVE_FLOW_BIN) $(NATIVE_FLOW_DEP) $(WORKFLOW_BIN) \
		$(WORKFLOW_DEP) $(RECOVERY_STAMP)

-include $(EXAMPLE_DEP) $(SMOKE_DEP) $(NATIVE_FLOW_DEP) $(WORKFLOW_DEP)
