#!/usr/bin/env bash
set -euo pipefail

echo "== Bigeye v2 CAA build =="

echo
echo "== Validate manifest =="
./validate_bigeye_v2_caa_manifest.sh

echo
echo "== Generate catalog =="
./generate_bigeye_v2_caa_catalog.sh

echo
echo "== Generate registry =="
./generate_bigeye_v2_caa_registry.sh

echo
echo "== Generate operation catalog =="
./generate_bigeye_v2_caa_operation_catalog.sh

echo
echo "== Generate IR =="
./generate_bigeye_v2_caa_ir.sh

echo
echo "== Generate operation graph =="
./generate_bigeye_v2_caa_operation_graph.sh

echo
echo "== Inspect operation graph =="
./inspect_bigeye_v2_caa_operation_graph.sh

echo
echo "== Validate operation graph =="
./validate_bigeye_v2_caa_operation_graph.sh

echo
echo "== Generate IR AssessmentCycle =="
./generate_bigeye_v2_caa_assessment_cycle_from_ir.sh

echo
echo "== Inspect execution plan =="
./inspect_bigeye_v2_caa_execution_plan.sh

echo
echo "== Inspect architecture =="
./inspect_bigeye_v2_caa_architecture.sh

echo
echo "== Run CAA suite =="
./run_bigeye_v2_caa_suite.sh

echo
echo "PASSED: Bigeye v2 CAA build"
