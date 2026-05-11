#!/usr/bin/env bash
set -euo pipefail

root="${1:-.}"

missing=0

check_file() {
  local path="$1"
  local label="$2"

  if [[ ! -f "${root}/${path}" ]]; then
    echo "Missing ${label}: expected ${root}/${path}"
    missing=1
  fi
}

check_dir() {
  local path="$1"
  local label="$2"

  if [[ ! -d "${root}/${path}" ]]; then
    echo "Missing ${label}: expected ${root}/${path}"
    missing=1
  fi
}

check_file "external/eigen/Eigen/Core" "Eigen Core header"
check_file "external/eigen/Eigen/Dense" "Eigen Dense header"
check_dir "external/had" "had dependency directory"
check_dir "external/LBFGSpp" "LBFGSpp dependency directory"
