#!/usr/bin/env bash
set -euo pipefail

echo "== Core uncertainty reporting API =="
grep -n "QUADRA_UNCERTAINTY_REPORTING_CORE\\|covariance_to_correlation_matrix\\|lognormal_delta_covariance\\|diagnose_covariance_matrix\\|correlation_decay_summary\\|ProjectionEnvelopeRow\\|summarize_samples" \
  core/uncertainty/reporting.hpp

echo
echo "== Run core reporting test =="
./run_uncertainty_reporting_test.sh
