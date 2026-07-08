#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/level08_joint_objective_check/bigeye_v2_level08_joint_objective_check.cpp \
  -o build/examples/bigeye_v2_level08_joint_objective_check

./build/examples/bigeye_v2_level08_joint_objective_check
