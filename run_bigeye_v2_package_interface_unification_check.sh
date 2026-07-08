#!/usr/bin/env bash
set -euo pipefail

./run_bigeye_v2_08_population_package_orchestrator_caa_check.sh
./run_bigeye_v2_09_fleet_package_orchestrator_caa_check.sh
./run_bigeye_v2_13_observation_package_caa_check.sh
./run_bigeye_v2_14_likelihood_package_caa_check.sh
./run_bigeye_v2_10_assessment_cycle_caa_check.sh

echo
echo "PASSED: Bigeye v2 CAA package interface unification check"
