#!/usr/bin/env bash
set -euo pipefail

output_dir="${1:-build/assessment_outputs}"
summary="${output_dir}/acceptance_summary.json"
recovery="${output_dir}/parameter_recovery.csv"
residuals="${output_dir}/residual_summary.csv"
audit="${output_dir}/laplace_gradient_audit.csv"

for required in "$summary" "$recovery" "$residuals" "$audit"; do
    if [[ ! -s "$required" ]]; then
        echo "simulation recovery check: missing or empty $required" >&2
        exit 1
    fi
done

grep -q '"converged": true' "$summary" || {
    echo "simulation recovery check: fit did not satisfy convergence rule" >&2
    exit 1
}

awk -F, '
    NR == 1 { next }
    $1 == "log_r0" && $5 > 0.15 { exit 1 }
    $1 ~ /^retention50_fleet_/ && $5 > 0.50 { exit 1 }
    $1 ~ /^retention_slope_fleet_/ && $5 > 0.30 { exit 1 }
    END { if (NR < 2) exit 1 }
' "$recovery" || {
    echo "simulation recovery check: key parameter recovery threshold failed" >&2
    exit 1
}

awk -F, '
    NR == 1 { next }
    $5 < 0.85 || $5 > 1.15 { exit 1 }
    END { if (NR < 2) exit 1 }
' "$residuals" || {
    echo "simulation recovery check: residual SDNR outside [0.85, 1.15]" >&2
    exit 1
}

awk -F, '
    NR == 1 { next }
    $7 != 1 || $8 != 1 { exit 1 }
    {
        difference = ($3 + $5) - $4
        if (difference < 0) difference = -difference
        scale = 1 + ($4 < 0 ? -$4 : $4)
        if (difference / scale > 1e-4) exit 1
    }
    END { if (NR < 2) exit 1 }
' "$audit" || {
    echo "simulation recovery check: Laplace finite-difference audit failed" >&2
    exit 1
}

echo "simulation recovery check: PASS"
