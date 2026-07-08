#!/usr/bin/env bash
set -euo pipefail

echo "== Bigeye v2 regression suite =="

./run_bigeye_v2_level00_life_history_check.sh
./run_bigeye_v2_level01_population_check.sh
./run_bigeye_v2_level02_selectivity_check.sh
./run_bigeye_v2_level03_mortality_check.sh
./run_bigeye_v2_level04_catch_check.sh
./run_bigeye_v2_level05_catch_likelihood_check.sh
./run_bigeye_v2_level06_index_likelihood_check.sh
./run_bigeye_v2_level07_agecomp_likelihood_check.sh
./run_bigeye_v2_level08_joint_objective_check.sh
./run_bigeye_v2_level09_optimizer_check.sh

echo
echo "PASSED: Bigeye v2 full regression suite"
