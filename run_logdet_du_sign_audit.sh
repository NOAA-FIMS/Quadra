#!/usr/bin/env bash
set -euo pipefail

grep -n "QUADRA_AUDIT_LOGDET_DU_SIGN\|minus_du_logdet_grad\|plus_du_logdet_grad" -A8 -B8 core/laplace.hpp

CXXFLAGS="${CXXFLAGS:-} -DQUADRA_AUDIT_LOGDET_DU_SIGN" \
  ./run_bigeye_level21_age_based_m_check.sh
