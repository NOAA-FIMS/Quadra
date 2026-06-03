#!/usr/bin/env bash
set -euo pipefail

mkdir -p examples/age_structured_recruitment

src="examples/age_structured_recruitment/benchmark_age_structured_no_plus_rolling_window.cpp"
dst="examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band.cpp"

if [[ ! -f "$src" ]]; then
  echo "ERROR: missing $src"
  echo "Run install_age_structured_no_plus_rolling_window_benchmark_v1.sh first."
  exit 1
fi

cp "$src" "$dst"

python3 - <<'PYEOF'
from pathlib import Path

p = Path("examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band.cpp")
s = p.read_text()

s = s.replace("#include <unordered_map>\n", "")
s = s.replace("#include <utility>\n", "")

start = s.find("static inline long long pack_key")
if start >= 0:
    end = s.find("\n\nEvalAll eval_all_rolling", start)
    if end >= 0:
        s = s[:start] + s[end+2:]

s = s.replace("EvalAll eval_all_rolling", "EvalAll eval_all_flat_band")
s = s.replace("eval_all_rolling", "eval_all_flat_band")

old = '''  // Store Hessian in a small-band map. This avoids a dense n x n matrix.
  std::unordered_map<long long, double> Hmap;
  Hmap.reserve(static_cast<std::size_t>(n * (2 * A + 1)));

  auto add_H = [&](int i, int j, double value) {
    if (i > j) std::swap(i, j);
    Hmap[pack_key(i, j)] += value;
  };
'''

new = '''  // Direct flat band storage. For this no-plus model, max distance is A.
  const int bw = A;
  std::vector<double> Hband(static_cast<std::size_t>((2 * bw + 1) * n), 0.0);

  auto add_H = [&](int i, int j, double value) {
    if (i < 0 || j < 0 || i >= n || j >= n) return;
    const int d = j - i;
    if (std::abs(d) > bw) {
      throw std::runtime_error("flat-band Hessian write outside bandwidth");
    }
    const int band = d + bw;
    Hband[static_cast<std::size_t>(band * n + i)] += value;
  };
'''

if old not in s:
    raise SystemExit("Could not find Hmap block")

s = s.replace(old, new, 1)

old2 = '''  for (const auto& kv : Hmap) {
    const long long key = kv.first;
    const int i = static_cast<int>(key >> 32);
    const int j = static_cast<int>(key & 0xffffffffu);
    const double v = kv.second;

    if (std::abs(v) > 1e-12) {
      triplets.emplace_back(i, j, v);
      if (i != j) triplets.emplace_back(j, i, v);
    }
  }
'''

new2 = '''  for (int i = 0; i < n; ++i) {
    for (int d = -bw; d <= bw; ++d) {
      const int j = i + d;
      if (j < 0 || j >= n) continue;

      const double v = Hband[static_cast<std::size_t>((d + bw) * n + i)];
      if (std::abs(v) > 1e-12) {
        triplets.emplace_back(i, j, v);
      }
    }
  }
'''

if old2 not in s:
    raise SystemExit("Could not find Hmap triplet conversion block")

s = s.replace(old2, new2, 1)

s = s.replace(
    "Quadra no-plus age-structured rolling-window analytic Laplace benchmark",
    "Quadra no-plus age-structured flat-band analytic Laplace benchmark"
)

p.write_text(s)
PYEOF

cat > run_quadra_age_structured_no_plus_flat_band_benchmark.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-10}"
LENGTHS="${2:-25,50,100,250,500,1000}"
AGES="${3:-10}"

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band.cpp \
  -o build/examples/benchmark_age_structured_no_plus_flat_band

./build/examples/benchmark_age_structured_no_plus_flat_band "$REPS" "$LENGTHS" "$AGES"
EOF

chmod +x run_quadra_age_structured_no_plus_flat_band_benchmark.sh

cat > run_quadra_flat_band_vs_tmb_age_structured_no_plus_benchmark.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-10}"
LENGTHS="${2:-25,50,100,250,500,1000}"
AGES="${3:-10}"

echo "== Quadra no-plus flat-band analytic Laplace =="
./run_quadra_age_structured_no_plus_flat_band_benchmark.sh "$REPS" "$LENGTHS" "$AGES"

echo
echo "== TMB no-plus age-structured AD/Laplace =="
./run_tmb_age_structured_no_plus_benchmark.sh "$REPS" "$LENGTHS" "$AGES"
EOF

chmod +x run_quadra_flat_band_vs_tmb_age_structured_no_plus_benchmark.sh

cat <<'EOF'

Installed no-plus age-structured flat-band benchmark.

Run:
  ./run_quadra_flat_band_vs_tmb_age_structured_no_plus_benchmark.sh 10 25,50,100,250,500,1000 10

Quick:
  ./run_quadra_age_structured_no_plus_flat_band_benchmark.sh 3 25,50 10

EOF
