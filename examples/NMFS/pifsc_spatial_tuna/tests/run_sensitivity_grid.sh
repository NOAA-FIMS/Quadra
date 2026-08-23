#!/usr/bin/env bash
set -euo pipefail

executable="${1:-build/advanced_tuna_spatial_assessment_example}"
grid_dir="${2:-build/sensitivity_grid}"
max_iterations="${QUADRA_TUNA_SENSITIVITY_MAX_ITERATIONS:-100}"
scenario_filter="${QUADRA_TUNA_SENSITIVITY_SCENARIOS:-all}"

mkdir -p "$grid_dir"
comparison="$grid_dir/sensitivity_comparison.csv"
echo "scenario,converged,iterations,gradient_norm,nll,depletion_terminal,observation_sd_multiplier,composition_effective_n,retention_prior_sd_multiplier,anchor_fleet,max_key_parameter_error,min_sdnr,max_sdnr" > "$comparison"

run_scenario() {
    scenario="$1"
    shift
    if [[ "$scenario_filter" != "all" && ",$scenario_filter," != *",$scenario,"* ]]; then
        return
    fi

    scenario_dir="$grid_dir/$scenario"
    mkdir -p "$scenario_dir"
    echo "sensitivity grid: running $scenario"
    env \
        QUADRA_TUNA_BASELINE_ONLY=1 \
        QUADRA_TUNA_MULTISTART=1 \
        QUADRA_TUNA_MAX_PHASE_ITERATIONS="$max_iterations" \
        QUADRA_TUNA_OUTPUT_DIR="$scenario_dir" \
        "$@" "$executable"

    run_values="$(awk -F, 'NR == 2 { print $0 }' "$scenario_dir/sensitivity_run_summary.csv")"
    max_error="$(awk -F, '
        NR > 1 && ($1 == "log_r0" || $1 ~ /^retention50_fleet_/ ||
                   $1 ~ /^retention_slope_fleet_/) {
            if ($5 > maximum) maximum = $5
        }
        END { print maximum + 0 }
    ' "$scenario_dir/parameter_recovery.csv")"
    sdnr_range="$(awk -F, '
        NR == 2 { minimum = $5; maximum = $5 }
        NR > 1 { if ($5 < minimum) minimum = $5; if ($5 > maximum) maximum = $5 }
        END { print minimum "," maximum }
    ' "$scenario_dir/residual_summary.csv")"
    echo "$scenario,$run_values,$max_error,$sdnr_range" >> "$comparison"
}

run_scenario baseline
run_scenario observation_tight QUADRA_TUNA_OBSERVATION_SD_MULTIPLIER=0.7
run_scenario observation_relaxed QUADRA_TUNA_OBSERVATION_SD_MULTIPLIER=1.4
run_scenario composition_neff_5 QUADRA_TUNA_COMPOSITION_EFFECTIVE_N=5
run_scenario composition_neff_25 QUADRA_TUNA_COMPOSITION_EFFECTIVE_N=25
run_scenario retention_prior_tight QUADRA_TUNA_RETENTION_PRIOR_SD_MULTIPLIER=0.5
run_scenario retention_prior_relaxed QUADRA_TUNA_RETENTION_PRIOR_SD_MULTIPLIER=2.0
run_scenario anchor_fleet_2 QUADRA_TUNA_ANCHOR_FLEET=2

echo "sensitivity grid: wrote $comparison"
