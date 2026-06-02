#!/usr/bin/env bash
set -euo pipefail

mkdir -p core/laplace tests .quadra_patch_backups

backup_file() {
  local f="$1"
  if [[ -f "$f" ]]; then
    local b=".quadra_patch_backups/$(echo "$f" | tr '/ ' '__').$(date +%Y%m%d_%H%M%S).bak"
    cp "$f" "$b"
    echo "Backed up $f -> $b"
  fi
}

backup_file core/laplace/laplace_exact_gradient_evaluator.hpp
backup_file tests/test_laplace_exact_gradient_evaluator.cpp
backup_file run_laplace_exact_gradient_evaluator_test.sh

if [[ ! -f core/laplace/exact_gradient_workspace.hpp ]]; then
  echo "ERROR: missing core/laplace/exact_gradient_workspace.hpp"
  exit 1
fi

cat > core/laplace/laplace_exact_gradient_evaluator.hpp <<'EOF'
#pragma once

#include <Eigen/Dense>

#include <cstddef>
#include <utility>
#include <vector>

#include "exact_gradient_workspace.hpp"

namespace quadra {
namespace laplace {

struct LaplaceExactGradientEvaluation {
    double objective = 0.0;
    Eigen::VectorXd gradient;
    Eigen::VectorXd trace_terms;
};

// Callback-driven evaluator scaffold over ExactGradientWorkspace.
//
// v1 intentionally avoids imposing a model concept. Existing evaluators can
// provide lambdas for graph construction, direction seeding, selected inverse
// access, and Hdot pattern extraction.
class LaplaceExactGradientEvaluator {
public:
    LaplaceExactGradientEvaluator() = default;

    LaplaceExactGradientEvaluator(const LaplaceExactGradientEvaluator&) = delete;
    LaplaceExactGradientEvaluator& operator=(const LaplaceExactGradientEvaluator&) = delete;

    LaplaceExactGradientEvaluator(LaplaceExactGradientEvaluator&&) = default;
    LaplaceExactGradientEvaluator& operator=(LaplaceExactGradientEvaluator&&) = default;

    template <class Builder,
              class DirectionProvider,
              class SelectedInverseAccessor>
    LaplaceExactGradientEvaluation Evaluate(
        Builder&& builder,
        std::vector<had::AReal>* fixed_effects,
        std::vector<had::AReal>* random_effects,
        std::size_t n_directions,
        DirectionProvider&& direction_provider,
        SelectedInverseAccessor&& selected_inverse,
        const std::vector<SparseHdotPatternEntry>& hdot_pattern,
        double joint_objective,
        double logdet_huu,
        const Eigen::VectorXd& joint_envelope_gradient) {
        workspace_.Build(
            std::forward<Builder>(builder),
            fixed_effects,
            random_effects);

        workspace_.PropagateBaseAdjoint();

        workspace_.SeedTotalDirections(
            n_directions,
            std::forward<DirectionProvider>(direction_provider));

        workspace_.PropagateDirectionalBatch();

        const ExactGradientEvaluation assembled =
            workspace_.AssembleExactGradient(
                joint_objective,
                logdet_huu,
                joint_envelope_gradient,
                std::forward<SelectedInverseAccessor>(selected_inverse),
                hdot_pattern);

        LaplaceExactGradientEvaluation out;
        out.objective = assembled.objective;
        out.gradient = assembled.gradient;
        out.trace_terms = assembled.trace_terms;
        return out;
    }

    ExactGradientWorkspace& Workspace() { return workspace_; }
    const ExactGradientWorkspace& Workspace() const { return workspace_; }

private:
    ExactGradientWorkspace workspace_;
};

}  // namespace laplace
}  // namespace quadra
EOF

cat > tests/test_laplace_exact_gradient_evaluator.cpp <<'EOF'
#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../core/laplace/laplace_exact_gradient_evaluator.hpp"

DECLARE_ADGRAPH()

namespace {

had::AReal objective(const std::vector<had::AReal>& theta,
                     const std::vector<had::AReal>& u) {
    const had::AReal& a = theta[0];
    const had::AReal& b = theta[1];

    had::AReal f(0.0);

    for (std::size_t i = 0; i < u.size(); ++i) {
        f = f + 0.5 * (u[i] - a) * (u[i] - a) + exp(b) * u[i] * u[i];
    }

    for (std::size_t i = 1; i < u.size(); ++i) {
        const had::AReal diff = u[i] - u[i - 1];
        f = f + 0.25 * diff * diff;
    }

    return f;
}

}  // namespace

int main() {
    using quadra::laplace::LaplaceExactGradientEvaluator;
    using quadra::laplace::MakeTridiagonalHdotPattern;

    std::vector<had::AReal> theta(2);
    std::vector<had::AReal> u(4);

    LaplaceExactGradientEvaluator evaluator;

    const Eigen::MatrixXd Hinv = Eigen::MatrixXd::Identity(4, 4);
    const auto pattern = MakeTridiagonalHdotPattern(4);
    const Eigen::VectorXd joint_gradient = Eigen::VectorXd::Zero(2);

    const auto result =
        evaluator.Evaluate(
            [&]() {
                theta[0] = had::AReal(0.2);
                theta[1] = had::AReal(-0.5);

                for (std::size_t i = 0; i < u.size(); ++i) {
                    u[i] = had::AReal(0.1 * static_cast<double>(i + 1));
                }

                return objective(theta, u);
            },
            &theta,
            &u,
            2,
            [](std::size_t k,
               Eigen::VectorXd& theta_direction,
               Eigen::VectorXd& random_direction) {
                theta_direction = Eigen::VectorXd::Zero(2);
                random_direction = Eigen::VectorXd::Zero(4);

                theta_direction[static_cast<int>(k)] = 1.0;

                for (int i = 0; i < random_direction.size(); ++i) {
                    random_direction[i] =
                        0.01 * static_cast<double>((i + 1) * (k + 1));
                }
            },
            [&](int row, int col) {
                return Hinv(row, col);
            },
            pattern,
            12.5,
            1.25,
            joint_gradient);

    if (std::abs(result.objective - 13.125) > 1.0e-12) {
        throw std::runtime_error("wrong Laplace exact-gradient objective");
    }

    if (result.gradient.size() != 2 || result.trace_terms.size() != 2) {
        throw std::runtime_error("wrong Laplace exact-gradient vector sizes");
    }

    if (!std::isfinite(result.gradient[0]) ||
        !std::isfinite(result.gradient[1])) {
        throw std::runtime_error("non-finite Laplace exact-gradient result");
    }

    std::cout << "laplace exact-gradient evaluator tests passed\n";
    std::cout << "objective = " << result.objective << "\n";
    std::cout << "gradient = " << result.gradient.transpose() << "\n";
    std::cout << "trace_terms = " << result.trace_terms.transpose() << "\n";

    return 0;
}
EOF

cat > run_laplace_exact_gradient_evaluator_test.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2}"

EIGEN_INCLUDE=""
if [[ -d external/Eigen ]]; then
  EIGEN_INCLUDE="-Iexternal/Eigen"
elif [[ -d core/eigen ]]; then
  EIGEN_INCLUDE="-Icore/eigen"
fi

mkdir -p build/tests

set -x
"${CXX}" ${CXXFLAGS} ${EIGEN_INCLUDE} -I.   tests/test_laplace_exact_gradient_evaluator.cpp   -o build/tests/test_laplace_exact_gradient_evaluator

./build/tests/test_laplace_exact_gradient_evaluator
EOF

chmod +x run_laplace_exact_gradient_evaluator_test.sh

cat <<'EOF'

Installed LaplaceExactGradientEvaluator v1.

Files added:
  core/laplace/laplace_exact_gradient_evaluator.hpp
  tests/test_laplace_exact_gradient_evaluator.cpp
  run_laplace_exact_gradient_evaluator_test.sh

Run:
  ./run_laplace_exact_gradient_evaluator_test.sh
  ./run_exact_gradient_workspace_test.sh

EOF
