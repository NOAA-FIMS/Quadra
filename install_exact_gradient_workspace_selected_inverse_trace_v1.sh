#!/usr/bin/env bash
set -euo pipefail

header="core/laplace/exact_gradient_workspace.hpp"
test_file="tests/test_exact_gradient_workspace.cpp"

mkdir -p .quadra_patch_backups

if [[ ! -f "$header" ]]; then
  echo "ERROR: missing $header"
  exit 1
fi

if [[ ! -f "$test_file" ]]; then
  echo "ERROR: missing $test_file"
  exit 1
fi

cp "$header" ".quadra_patch_backups/exact_gradient_workspace.hpp.selected_inverse_trace.$(date +%Y%m%d_%H%M%S).bak"
cp "$test_file" ".quadra_patch_backups/test_exact_gradient_workspace.cpp.selected_inverse_trace.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/quadra_selected_inverse_trace_patch.py <<'PYEOF'
from pathlib import Path

header = Path("core/laplace/exact_gradient_workspace.hpp")
s = header.read_text()

if "TraceTermsSelectedInverse" not in s:
    marker = "    HadGraphWorkspace& HadWorkspace() { return had_workspace_; }"
    idx = s.find(marker)
    if idx < 0:
        raise SystemExit("Could not find HadWorkspace marker")

    method = """
    template <class SelectedInverseAccessor>
    Eigen::VectorXd TraceTermsSelectedInverse(
        SelectedInverseAccessor&& selected_inverse,
        const std::vector<SparseHdotPatternEntry>& pattern) {
        RequireBuilt();

        Eigen::VectorXd traces =
            Eigen::VectorXd::Zero(static_cast<int>(n_directions_));

        had_workspace_.Activate();

        for (std::size_t k = 0; k < n_directions_; ++k) {
            double trace = 0.0;

            for (const auto& entry : pattern) {
                CheckRandomIndex(entry.row);
                CheckRandomIndex(entry.col);

                const double hdot =
                    had::GetAdjointDotBatch(
                        (*random_effects_)[static_cast<std::size_t>(entry.row)],
                        (*random_effects_)[static_cast<std::size_t>(entry.col)],
                        static_cast<int>(k));

                const double hinv = selected_inverse(entry.row, entry.col);

                if (entry.row == entry.col) {
                    trace += hinv * hdot;
                } else {
                    trace += 2.0 * hinv * hdot;
                }
            }

            traces[static_cast<int>(k)] = trace;
        }

        return traces;
    }

"""
    s = s[:idx] + method + s[idx:]
    header.write_text(s)

test = Path("tests/test_exact_gradient_workspace.cpp")
t = test.read_text()

if "TraceTermsSelectedInverse" not in t:
    old = """    const Eigen::VectorXd traces = workspace.TraceTerms(Hinv, pattern);

    if (traces.size() != 2) {
        throw std::runtime_error("TraceTerms returned wrong number of directions");
    }

    if (!std::isfinite(traces[0]) || !std::isfinite(traces[1])) {
        throw std::runtime_error("TraceTerms returned non-finite values");
    }"""

    new = """    const Eigen::VectorXd traces = workspace.TraceTerms(Hinv, pattern);

    const Eigen::VectorXd selected_inverse_traces =
        workspace.TraceTermsSelectedInverse(
            [&](int row, int col) {
                return Hinv(row, col);
            },
            pattern);

    if (traces.size() != 2) {
        throw std::runtime_error("TraceTerms returned wrong number of directions");
    }

    if (selected_inverse_traces.size() != traces.size()) {
        throw std::runtime_error(
            "TraceTermsSelectedInverse returned wrong number of directions");
    }

    if ((traces - selected_inverse_traces).cwiseAbs().maxCoeff() > 1.0e-12) {
        throw std::runtime_error(
            "TraceTermsSelectedInverse does not match dense TraceTerms");
    }

    if (!std::isfinite(traces[0]) || !std::isfinite(traces[1])) {
        throw std::runtime_error("TraceTerms returned non-finite values");
    }"""

    if old not in t:
        raise SystemExit("Could not find TraceTerms test block")

    t = t.replace(old, new, 1)
    test.write_text(t)
PYEOF

python3 /tmp/quadra_selected_inverse_trace_patch.py

cat <<'EOF'

Installed ExactGradientWorkspace selected-inverse trace API.

Patched:
  core/laplace/exact_gradient_workspace.hpp
  tests/test_exact_gradient_workspace.cpp

Run:
  ./run_exact_gradient_workspace_test.sh
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

EOF
