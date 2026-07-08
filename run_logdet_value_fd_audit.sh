#!/usr/bin/env bash
set -euo pipefail

CXXFLAGS="-DQUADRA_AUDIT_LOGDET_VALUE_FD ${CXXFLAGS:-}" \
  ./run_bigeye_level21_age_based_m_check.sh > level21_logdet_value_fd_audit.log 2>&1

grep -m1 -A8 "QUADRA_AUDIT_LOGDET_VALUE_FD" level21_logdet_value_fd_audit.log
