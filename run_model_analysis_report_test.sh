#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/tests

set -x
c++ -std=c++17 -O2 -DNDEBUG -g \
  -Iexternal/Eigen \
  -I. \
  tests/test_model_analysis_report.cpp \
  -o build/tests/test_model_analysis_report

./build/tests/test_model_analysis_report
