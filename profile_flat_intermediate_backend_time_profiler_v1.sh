#!/usr/bin/env bash
set -euo pipefail

# profile_flat_intermediate_backend_time_profiler_v1.sh
#
# Profiles the current flat-intermediate exact-gradient reuse path on macOS.
#
# Produces:
#   build benchmark
#   normal timing run
#   .trace package
#   exported XML
#   quick hotspot grep summary
#
# Usage:
#   ./profile_flat_intermediate_backend_time_profiler_v1.sh 10
#
# Optional env:
#   PROFILE_M=500   # benchmark currently profiles all cases unless --reuse-only mode exists
#   PROFILE_SECONDS unused; xctrace stops when target exits.

REPS="${1:-10}"
BENCH="build/benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse"
SRC="benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp"

STAMP="$(date +%Y%m%d_%H%M%S)"
OUT_PREFIX="quadra_flat_intermediate_profile_${STAMP}"
TRACE="${OUT_PREFIX}.trace"
XML="${OUT_PREFIX}.xml"
TIMING_LOG="${OUT_PREFIX}_timing.log"
GREP_LOG="${OUT_PREFIX}_hotspots.txt"

if [[ ! -f "$SRC" ]]; then
  echo "ERROR: missing $SRC"
  exit 1
fi

mkdir -p build/benchmarks profiles

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O3 -DNDEBUG -g}"
EIGEN_INCLUDE=""
if [[ -d external/Eigen ]]; then
  EIGEN_INCLUDE="-Iexternal/Eigen"
elif [[ -d core/eigen ]]; then
  EIGEN_INCLUDE="-Icore/eigen"
fi

LBFGS_INCLUDE=""
if [[ -d external/LBFGSpp/include ]]; then
  LBFGS_INCLUDE="-Iexternal/LBFGSpp/include"
elif [[ -d external/LBFGSpp ]]; then
  LBFGS_INCLUDE="-Iexternal/LBFGSpp"
fi

echo "Building benchmark..."
set -x
"${CXX}" ${CXXFLAGS} ${EIGEN_INCLUDE} ${LBFGS_INCLUDE} -I. \
  "$SRC" \
  -o "$BENCH"
set +x

echo
echo "Running normal timing benchmark..."
"./$BENCH" "$REPS" 2>&1 | tee "$TIMING_LOG"

echo
echo "Recording Time Profiler trace..."
rm -rf "$TRACE"

# Prefer reuse-only if benchmark supports it; harmless if it ignores the flag.
# The output will show whether mode = reuse-only profiling was recognized.
set -x
xcrun xctrace record \
  --template "Time Profiler" \
  --output "$TRACE" \
  --launch -- "./$BENCH" "$REPS" --reuse-only
set +x

echo
echo "Exporting trace to XML..."
set -x
xcrun xctrace export \
  --input "$TRACE" \
  --output "$XML" \
  --xpath '/trace-toc/run[@number="1"]/data/table[@schema="time-sample"]'
set +x

echo
echo "Generating hotspot grep summary..."
{
  echo "Hotspot grep summary for $XML"
  echo

  echo "== Core reverse / flat backend symbols =="
  grep -E \
    "PropagateAdjointDirectionalBatch|AddFlatIntermediate|TryGetFlatIntermediate|BatchDirectionalFlatAccumulator|IntermediateEdgeSlotRegistry|PushEdgeDotBatchValue|BatchDirection|ComputeBatchActiveDirectionMasks|BTree::|EnsureBatch" \
    "$XML" | head -250 || true

  echo
  echo "== C++ library / allocator / memory movement =="
  grep -E \
    "memset|memcpy|malloc|free|operator new|allocator|vector|unordered_map|_platform_mem" \
    "$XML" | head -150 || true

  echo
  echo "== Eigen / sparse factorization symbols =="
  grep -E \
    "Eigen::|SimplicialLDLT|SimplicialLLT|solve|factor|chol|ldlt" \
    "$XML" | head -150 || true
} | tee "$GREP_LOG"

mkdir -p profiles
mv "$TRACE" profiles/
mv "$XML" profiles/
mv "$TIMING_LOG" profiles/
mv "$GREP_LOG" profiles/

cat <<EOF

Profile complete.

Files:
  profiles/$TRACE
  profiles/$XML
  profiles/$TIMING_LOG
  profiles/$GREP_LOG

Next:
  paste the top 100-200 lines from:
    profiles/$GREP_LOG

Useful manual checks:
  grep -E "PropagateAdjointDirectionalBatch|BatchDirectionalFlatAccumulator|BTree::|memset|vector|unordered_map" profiles/$XML | head -200

EOF
