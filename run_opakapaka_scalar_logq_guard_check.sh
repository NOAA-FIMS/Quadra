#!/usr/bin/env bash
set -euo pipefail

./inspect_opakapaka_scalar_logq_guard.sh

echo
echo "== Remove stale scalar-only uncertainty files before run =="
rm -f \
  examples/NMFS/pifsc_opakapaka/outputs/uncertainty_summary.csv \
  examples/NMFS/pifsc_opakapaka/outputs/covariance_matrix.csv \
  examples/NMFS/pifsc_opakapaka/outputs/correlation_matrix.csv \
  examples/NMFS/pifsc_opakapaka/outputs/standard_errors.csv \
  examples/NMFS/pifsc_opakapaka/outputs/confidence_intervals.csv

echo
echo "== O3 build after scalar logq uncertainty guard =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_opakapaka/quadra/opakapaka.cpp \
  examples/NMFS/pifsc_opakapaka/quadra/opakapaka_adgraph_global.cpp \
  -o build/examples/pifsc_opakapaka_scalar_logq_guard_check

echo
echo "== Run Opakapaka guard check binary =="
./build/examples/pifsc_opakapaka_scalar_logq_guard_check

echo
echo "== Check scalar-only uncertainty output files =="
for f in \
  examples/NMFS/pifsc_opakapaka/outputs/uncertainty_summary.csv \
  examples/NMFS/pifsc_opakapaka/outputs/covariance_matrix.csv \
  examples/NMFS/pifsc_opakapaka/outputs/correlation_matrix.csv \
  examples/NMFS/pifsc_opakapaka/outputs/standard_errors.csv \
  examples/NMFS/pifsc_opakapaka/outputs/confidence_intervals.csv
do
  if [[ -f "$f" ]]; then
    echo "ERROR: scalar-only file was written for multi-fixed fit: $f"
    exit 1
  else
    echo "OK: not written for multi-fixed fit: $f"
  fi
done
