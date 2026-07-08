#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/tests

c++ -std=c++17 -O3 \
  -I. -Iexternal/eigen \
  tests/laplace/test_laplace_logdet_theta_dependent_regression.cpp \
  -o build/tests/test_laplace_logdet_theta_dependent_regression

./build/tests/test_laplace_logdet_theta_dependent_regression
