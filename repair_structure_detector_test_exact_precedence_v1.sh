#!/usr/bin/env bash
set -euo pipefail

target="tests/test_structure_detector_registry.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/test_structure_detector_registry.cpp.exact_precedence_expectation.$(date +%Y%m%d_%H%M%S).bak"

python3 - "$target" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()

old = '''void test_small_matrix_dense_preference() {
  StructureDetectorOptions opts;
  opts.prefer_dense_for_small_matrices = true;
  opts.dense_size_cutoff = 16;
  opts.banded_width_cutoff = 8;
  opts.dense_fill_ratio = 0.40;

  StructureDetector det(opts);

  std::vector<Eigen::Triplet<double>> t;
  t.emplace_back(0, 0, 2.0);
  t.emplace_back(1, 1, 3.0);
  t.emplace_back(2, 2, 4.0);

  auto H = make_sparse(3, t);
  const auto rec = det.Analyze(H);

  expect_true(rec.backend == LaplaceBackendKind::DenseLDLT,
              "small matrix dense preference");
}'''

new = '''void test_small_matrix_exact_structure_precedence() {
  StructureDetectorOptions opts;
  opts.prefer_dense_for_small_matrices = true;
  opts.dense_size_cutoff = 16;
  opts.banded_width_cutoff = 8;
  opts.dense_fill_ratio = 0.40;

  StructureDetector det(opts);

  std::vector<Eigen::Triplet<double>> t;
  t.emplace_back(0, 0, 2.0);
  t.emplace_back(1, 1, 3.0);
  t.emplace_back(2, 2, 4.0);

  auto H = make_sparse(3, t);
  const auto rec = det.Analyze(H);

  expect_true(rec.backend == LaplaceBackendKind::Diagonal,
              "small exact diagonal structure takes precedence over dense preference");
}'''

if old not in s:
    raise SystemExit("Could not find old small matrix dense preference test")

s = s.replace(old, new, 1)
s = s.replace("test_small_matrix_dense_preference();",
              "test_small_matrix_exact_structure_precedence();")

p.write_text(s)
PY

cat <<'EOF'

Updated structure detector test expectations for exact-structure precedence.

Run:
  ./run_structure_detector_registry_test.sh
  ./run_laplace_backend_factory_test.sh

Then rebuild/rerun the state-space example.

EOF
