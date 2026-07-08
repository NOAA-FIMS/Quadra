#!/usr/bin/env bash
set -euo pipefail
OBJ="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/objective/bigeye_quadra_objective.hpp"
BAK="$OBJ.before_level21_fixed_m_offsets_diagnostic"

if [[ ! -f "$BAK" ]]; then
  echo "ERROR: backup not found: $BAK"
  exit 1
fi

cp "$BAK" "$OBJ"
echo "Restored Level 21 objective from:"
echo "  $BAK"
