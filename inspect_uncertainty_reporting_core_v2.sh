#!/usr/bin/env bash
set -euo pipefail

echo "== Eigen includes =="
grep -n "#include <Eigen/" core/uncertainty/reporting.hpp

echo
echo "== Run core reporting test =="
./run_uncertainty_reporting_test.sh
