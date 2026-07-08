#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/tests

clang++ -std=c++17 -g -I"external/eigen/" \
  tests/test_laplace_structure_report.cpp \
  -o build/tests/test_laplace_structure_report

build/tests/test_laplace_structure_report
