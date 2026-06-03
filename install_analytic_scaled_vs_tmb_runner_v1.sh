#!/usr/bin/env bash
set -euo pipefail

# install_analytic_scaled_vs_tmb_runner_v1.sh
#
# Creates a clean comparison runner using the analytic scaled latent-state
# tridiagonal benchmark instead of the finite-difference scaled benchmark.

quadra_src="examples/state_space_surplus_production/benchmark_latent_tridiagonal_analytic_scaled.cpp"

if [[ ! -f "$quadra_src" ]]; then
  echo "ERROR: missing $quadra_src"
  echo
  echo "Available scaled benchmarks:"
  find examples/state_space_surplus_production -maxdepth 1 -name '*scaled*.cpp' -print
  exit 1
fi

cat > run_quadra_analytic_vs_tmb_scaled_fixed_theta_benchmark.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-10}"
LENGTHS="${2:-25,50,100,250,500,1000}"

mkdir -p build/examples

echo "== Quadra scaled analytic latent-state tridiagonal =="
set -x
c++ -std=c++17 -O2 -DNDEBUG -g \
  -Iexternal/Eigen \
  -Iexternal/LBFGSpp/include \
  -Iexamples/state_space_surplus_production \
  -Iexamples/surplus_production \
  examples/state_space_surplus_production/benchmark_latent_tridiagonal_analytic_scaled.cpp \
  -o build/examples/benchmark_latent_tridiagonal_analytic_scaled

./build/examples/benchmark_latent_tridiagonal_analytic_scaled "$REPS" "$LENGTHS"
set +x

echo
echo "== TMB scaled AD/Laplace =="

if [[ -x ./run_tmb_scaled_state_space_surplus_benchmark.sh ]]; then
  ./run_tmb_scaled_state_space_surplus_benchmark.sh "$REPS" "$LENGTHS"
elif [[ -x ./run_quadra_vs_tmb_scaled_fixed_theta_benchmark.sh ]]; then
  ./run_quadra_vs_tmb_scaled_fixed_theta_benchmark.sh "$REPS" "$LENGTHS" \
    | sed -n '/== TMB scaled AD\\/Laplace ==/,$p'
else
  echo "Could not find a TMB scaled benchmark runner."
  echo "Known candidates:"
  find . -maxdepth 2 -name '*tmb*scaled*' -o -name '*scaled*tmb*'
  exit 1
fi
EOF

chmod +x run_quadra_analytic_vs_tmb_scaled_fixed_theta_benchmark.sh

cat <<'EOF'

Installed analytic scaled Quadra vs TMB runner.

Run:
  ./run_quadra_analytic_vs_tmb_scaled_fixed_theta_benchmark.sh 10 25,50,100,250,500,1000

This avoids the finite-difference scaled path and uses:
  examples/state_space_surplus_production/benchmark_latent_tridiagonal_analytic_scaled.cpp

EOF
