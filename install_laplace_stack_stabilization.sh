#!/usr/bin/env bash
set -euo pipefail

# Run from the Quadra repository root.
#
# Stabilization patch:
# - adds Makefile targets for implicit derivative test/demo
# - wires test_laplace_implicit_derivatives into run-tests if possible
# - adds a focused laplace-stack-check target

python3 - <<'PY'
from pathlib import Path

p = Path("Makefile")
if not p.exists():
    raise SystemExit("Makefile not found")

s = p.read_text()

if "tests/test_laplace_implicit_derivatives: tests/test_laplace_implicit_derivatives.cpp" not in s:
    s += """

tests/test_laplace_implicit_derivatives: tests/test_laplace_implicit_derivatives.cpp
\t$(CXX) $(CXXFLAGS) -o $@ $<
"""

if "examples/laplace_implicit_derivatives_demo: examples/laplace_implicit_derivatives_demo.cpp" not in s:
    s += """

examples/laplace_implicit_derivatives_demo: examples/laplace_implicit_derivatives_demo.cpp
\t$(CXX) $(CXXFLAGS) -o $@ $<
"""

if "./tests/test_laplace_implicit_derivatives" not in s and "run-tests:" in s:
    markers = [
        "\t./tests/test_laplace_exact_gradient",
        "\t./tests/test_laplace_optimizer",
        "\t./tests/test_laplace_fixed_gradient",
        "\t./tests/test_laplace_objective",
    ]

    for marker in markers:
        if marker in s:
            s = s.replace(marker, marker + "\n\t./tests/test_laplace_implicit_derivatives", 1)
            break

if "laplace-stack-check:" not in s:
    s += """

.PHONY: laplace-stack-check
laplace-stack-check: \\
\ttests/test_random_effect_objective \\
\ttests/test_random_effect_hessian \\
\ttests/test_random_effect_newton \\
\ttests/test_laplace_objective \\
\ttests/test_laplace_fixed_gradient \\
\ttests/test_laplace_optimizer \\
\ttests/test_laplace_exact_gradient \\
\ttests/test_laplace_implicit_derivatives
\t./tests/test_random_effect_objective
\t./tests/test_random_effect_hessian
\t./tests/test_random_effect_newton
\t./tests/test_laplace_objective
\t./tests/test_laplace_fixed_gradient
\t./tests/test_laplace_optimizer
\t./tests/test_laplace_exact_gradient
\t./tests/test_laplace_implicit_derivatives
"""

p.write_text(s.rstrip() + "\n")
print("Updated Makefile with implicit derivative targets and laplace-stack-check.")
PY

echo
echo "Verify with:"
echo "  make tests/test_laplace_implicit_derivatives && ./tests/test_laplace_implicit_derivatives"
echo "  make laplace-stack-check"
echo "  make run-tests"
