#!/usr/bin/env bash
set -euo pipefail

log="level21_logdet_component_audit.log"

CXXFLAGS="${CXXFLAGS:-} -DQUADRA_AUDIT_LOGDET_COMPONENTS" \
  ./run_bigeye_level21_age_based_m_check.sh > "$log" 2>&1

echo "wrote: $log"
echo
grep -m1 -A8 "QUADRA_AUDIT_LOGDET_COMPONENTS" "$log"
