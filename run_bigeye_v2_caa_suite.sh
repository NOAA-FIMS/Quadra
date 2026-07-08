#!/usr/bin/env bash
set -euo pipefail

echo "== Bigeye v2 CAA suite =="

./run_bigeye_v2_01_life_history_caa_check.sh
./run_bigeye_v2_level02_selectivity_caa_check.sh
./run_bigeye_v2_03_fleet_package_caa_check.sh
./run_bigeye_v2_04_population_recruitment_caa_check.sh
./run_bigeye_v2_05_population_survival_caa_check.sh
./run_bigeye_v2_06_population_aging_caa_check.sh
./run_bigeye_v2_07_population_package_caa_check.sh
./run_bigeye_v2_08_population_package_orchestrator_caa_check.sh
./run_bigeye_v2_09_fleet_package_orchestrator_caa_check.sh
./run_bigeye_v2_10_assessment_cycle_caa_check.sh
./run_bigeye_v2_11_observation_biomass_index_caa_check.sh
./run_bigeye_v2_12_observation_agecomp_caa_check.sh
./run_bigeye_v2_13_observation_package_caa_check.sh
./run_bigeye_v2_14_likelihood_package_caa_check.sh
./run_bigeye_v2_15_assessment_cycle_objective_caa_check.sh
./run_bigeye_v2_16_movement_package_caa_check.sh
./run_bigeye_v2_package_interface_unification_check.sh

echo
echo "PASSED: Bigeye v2 CAA suite"
