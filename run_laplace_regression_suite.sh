#!/usr/bin/env bash
set -euo pipefail

./run_laplace_regression_test.sh
./run_laplace_theta_dependent_regression_test.sh

echo
echo "PASSED: full Laplace regression suite"
