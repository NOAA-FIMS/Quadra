#!/usr/bin/env bash
set -euo pipefail

# install_exact_gradient_workspace_v1.sh
#
# Adds first production-facing exact-gradient workspace scaffold:
#
#   core/laplace/exact_gradient_workspace.hpp
#
# This wraps:
#   quadra::HadGraphWorkspace
#   direction seeding
#   batched directional propagation
#   sparse-pattern Hdot extraction
#
# The class is intentionally generic enough to support model-specific builders
# while keeping RW1 benchmark validation straightforward.

mkdir -p core/laplace tests .quadra_patch_backups

backup_file() {
  local f="$1"
  if [[ -f "$f" ]]; then
    local b=".quadra_patch_backups/$(echo "$f" | tr '/ ' '__').$(date +%Y%m%d_%H%M%S).bak"
    cp "$f" "$b"
    echo "Backed up $f -> $b"
  fi
}

backup_file core/laplace/exact_gradient_workspace.hpp
backup_file tests/test_exact_gradient_workspace.cpp
backup_file run_exact_gradient_workspace_test.sh

if [[ ! -f core/had_graph_workspace.hpp ]]; then
  echo "ERROR: missing core/had_graph_workspace.hpp"
  echo "Install HadGraphWorkspace first."
  exit 1
fi

cat > core/laplace/exact_gradient_workspace.hpp <<'EOF'
#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../had_graph_workspace.hpp"

namespace quadra {
namespace laplace {

struct SparseHdotPatternEntry {
    int row = 0;
    int col = 0;

    SparseHdotPatternEntry() = default;
    SparseHdotPatternEntry(int row_, int col_) : row(row_), col(col_) {}
};

// Production-facing scaffold for exact Laplace gradient Hdot reuse.
//
// Responsibilities:
//   - own/reuse a HAD graph via HadGraphWorkspace
//   - seed batched total derivative directions
//   - run batched directional reverse propagation
//   - extract Hdot values over a sparse pattern
//
// Non-responsibilities in v1:
//   - solving uhat
//   - factorization ownership
//   - logdet/objective assembly
//
// Those stay with higher-level Laplace evaluators.
class ExactGradientWorkspace {
public:
    ExactGradientWorkspace() = default;

    ExactGradientWorkspace(const ExactGradientWorkspace&) = delete;
    ExactGradientWorkspace& operator=(const ExactGradientWorkspace&) = delete;

    ExactGradientWorkspace(ExactGradientWorkspace&&) = default;
    ExactGradientWorkspace& operator=(ExactGradientWorkspace&&) = default;

    template <class Builder>
    had::AReal Build(Builder&& builder,
                     std::vector<had::AReal>* fixed_effects,
                     std::vector<had::AReal>* random_effects) {
        if (fixed_effects == nullptr || random_effects == nullptr) {
            throw std::invalid_argument(
                "ExactGradientWorkspace::Build requires non-null variable handles.");
        }

        fixed_effects_ = fixed_effects;
        random_effects_ = random_effects;

        output_ = had_workspace_.Build(std::forward<Builder>(builder));
        built_ = true;

        return output_;
    }

    void PropagateBaseAdjoint() {
        RequireBuilt();
        had_workspace_.PropagateAdjoint(output_.varId);
    }

    void ResizeDirectionalBatch(std::size_t n_directions) {
        RequireBuilt();
        n_directions_ = n_directions;
        had_workspace_.ResizeDirectionalBatch(n_directions);
    }

    template <class DirectionProvider>
    void SeedTotalDirections(std::size_t n_directions,
                             DirectionProvider&& direction_provider) {
        RequireBuilt();
        ResizeDirectionalBatch(n_directions);

        const std::size_t n_fixed = fixed_effects_->size();
        const std::size_t n_random = random_effects_->size();

        had_workspace_.Activate();

        for (std::size_t k = 0; k < n_directions; ++k) {
            Eigen::VectorXd theta_direction;
            Eigen::VectorXd random_direction;

            direction_provider(k, theta_direction, random_direction);

            if (theta_direction.size() != static_cast<int>(n_fixed)) {
                throw std::invalid_argument(
                    "ExactGradientWorkspace::SeedTotalDirections theta direction size mismatch.");
            }
            if (random_direction.size() != static_cast<int>(n_random)) {
                throw std::invalid_argument(
                    "ExactGradientWorkspace::SeedTotalDirections random direction size mismatch.");
            }

            for (std::size_t j = 0; j < n_fixed; ++j) {
                had::SetARealDotBatch((*fixed_effects_)[j],
                                      static_cast<int>(k),
                                      theta_direction[static_cast<int>(j)]);
            }

            for (std::size_t i = 0; i < n_random; ++i) {
                had::SetARealDotBatch((*random_effects_)[i],
                                      static_cast<int>(k),
                                      random_direction[static_cast<int>(i)]);
            }
        }
    }

    void PropagateDirectionalBatch() {
        RequireBuilt();
        had_workspace_.PropagateAdjointDirectionalBatch();
    }

    Eigen::MatrixXd ExtractHdotDense(
        std::size_t direction_index,
        const std::vector<SparseHdotPatternEntry>& pattern) {
        RequireBuilt();

        if (direction_index >= n_directions_) {
            throw std::out_of_range(
                "ExactGradientWorkspace::ExtractHdotDense direction_index out of range.");
        }

        const int n_random = static_cast<int>(random_effects_->size());
        Eigen::MatrixXd out = Eigen::MatrixXd::Zero(n_random, n_random);

        had_workspace_.Activate();

        for (const auto& entry : pattern) {
            CheckRandomIndex(entry.row);
            CheckRandomIndex(entry.col);

            const double value =
                had::GetAdjointDotBatch(
                    (*random_effects_)[static_cast<std::size_t>(entry.row)],
                    (*random_effects_)[static_cast<std::size_t>(entry.col)],
                    static_cast<int>(direction_index));

            out(entry.row, entry.col) = value;
            out(entry.col, entry.row) = value;
        }

        return out;
    }

    std::vector<Eigen::Triplet<double>> ExtractHdotTriplets(
        std::size_t direction_index,
        const std::vector<SparseHdotPatternEntry>& pattern) {
        RequireBuilt();

        if (direction_index >= n_directions_) {
            throw std::out_of_range(
                "ExactGradientWorkspace::ExtractHdotTriplets direction_index out of range.");
        }

        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(pattern.size() * 2);

        had_workspace_.Activate();

        for (const auto& entry : pattern) {
            CheckRandomIndex(entry.row);
            CheckRandomIndex(entry.col);

            const double value =
                had::GetAdjointDotBatch(
                    (*random_effects_)[static_cast<std::size_t>(entry.row)],
                    (*random_effects_)[static_cast<std::size_t>(entry.col)],
                    static_cast<int>(direction_index));

            triplets.emplace_back(entry.row, entry.col, value);

            if (entry.row != entry.col) {
                triplets.emplace_back(entry.col, entry.row, value);
            }
        }

        return triplets;
    }

    Eigen::VectorXd TraceTerms(
        const Eigen::MatrixXd& Hinv,
        const std::vector<SparseHdotPatternEntry>& pattern) {
        RequireBuilt();

        const int n_random = static_cast<int>(random_effects_->size());
        if (Hinv.rows() != n_random || Hinv.cols() != n_random) {
            throw std::invalid_argument(
                "ExactGradientWorkspace::TraceTerms Hinv dimension mismatch.");
        }

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

                if (entry.row == entry.col) {
                    trace += Hinv(entry.row, entry.col) * hdot;
                } else {
                    trace += 2.0 * Hinv(entry.row, entry.col) * hdot;
                }
            }

            traces[static_cast<int>(k)] = trace;
        }

        return traces;
    }

    HadGraphWorkspace& HadWorkspace() { return had_workspace_; }
    const HadGraphWorkspace& HadWorkspace() const { return had_workspace_; }

    std::size_t DirectionCount() const { return n_directions_; }

private:
    void RequireBuilt() const {
        if (!built_) {
            throw std::logic_error(
                "ExactGradientWorkspace used before Build.");
        }
    }

    void CheckRandomIndex(int index) const {
        if (index < 0 ||
            index >= static_cast<int>(random_effects_->size())) {
            throw std::out_of_range(
                "ExactGradientWorkspace random-effect index out of range.");
        }
    }

    HadGraphWorkspace had_workspace_;
    had::AReal output_;
    std::vector<had::AReal>* fixed_effects_ = nullptr;
    std::vector<had::AReal>* random_effects_ = nullptr;
    std::size_t n_directions_ = 0;
    bool built_ = false;
};

inline std::vector<SparseHdotPatternEntry>
MakeTridiagonalHdotPattern(int n) {
    if (n < 0) {
        throw std::invalid_argument(
            "MakeTridiagonalHdotPattern requires nonnegative n.");
    }

    std::vector<SparseHdotPatternEntry> pattern;
    pattern.reserve(static_cast<std::size_t>(2 * n));

    for (int i = 0; i < n; ++i) {
        pattern.emplace_back(i, i);
        if (i > 0) {
            pattern.emplace_back(i, i - 1);
        }
    }

    return pattern;
}

}  // namespace laplace
}  // namespace quadra
EOF

cat > tests/test_exact_gradient_workspace.cpp <<'EOF'
#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../core/laplace/exact_gradient_workspace.hpp"

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
    using quadra::laplace::ExactGradientWorkspace;
    using quadra::laplace::MakeTridiagonalHdotPattern;

    std::vector<had::AReal> theta(2);
    std::vector<had::AReal> u(4);

    ExactGradientWorkspace workspace;

    had::AReal f = workspace.Build(
        [&]() {
            theta[0] = had::AReal(0.2);
            theta[1] = had::AReal(-0.5);

            for (std::size_t i = 0; i < u.size(); ++i) {
                u[i] = had::AReal(0.1 * static_cast<double>(i + 1));
            }

            return objective(theta, u);
        },
        &theta,
        &u);

    (void)f;

    workspace.PropagateBaseAdjoint();

    workspace.SeedTotalDirections(
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
        });

    workspace.PropagateDirectionalBatch();

    const auto pattern = MakeTridiagonalHdotPattern(4);

    Eigen::MatrixXd Hinv = Eigen::MatrixXd::Identity(4, 4);

    const Eigen::VectorXd traces = workspace.TraceTerms(Hinv, pattern);

    if (traces.size() != 2) {
        throw std::runtime_error("TraceTerms returned wrong number of directions");
    }

    if (!std::isfinite(traces[0]) || !std::isfinite(traces[1])) {
        throw std::runtime_error("TraceTerms returned non-finite values");
    }

    const Eigen::MatrixXd hdot0 = workspace.ExtractHdotDense(0, pattern);

    if (hdot0.rows() != 4 || hdot0.cols() != 4) {
        throw std::runtime_error("ExtractHdotDense returned wrong dimensions");
    }

    const auto triplets = workspace.ExtractHdotTriplets(1, pattern);

    if (triplets.empty()) {
        throw std::runtime_error("ExtractHdotTriplets returned empty result");
    }

    std::cout << "exact gradient workspace tests passed\n";
    std::cout << "vertices = "
              << workspace.HadWorkspace().VertexCount() << "\n";
    std::cout << "trace 0 = " << traces[0] << "\n";
    std::cout << "trace 1 = " << traces[1] << "\n";
    std::cout << "triplets = " << triplets.size() << "\n";

    return 0;
}
EOF

cat > run_exact_gradient_workspace_test.sh <<'EOF'
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
"${CXX}" ${CXXFLAGS} ${EIGEN_INCLUDE} -I. \
  tests/test_exact_gradient_workspace.cpp \
  -o build/tests/test_exact_gradient_workspace

./build/tests/test_exact_gradient_workspace
EOF

chmod +x run_exact_gradient_workspace_test.sh

cat <<'EOF'

Installed ExactGradientWorkspace v1.

Files added:
  core/laplace/exact_gradient_workspace.hpp
  tests/test_exact_gradient_workspace.cpp
  run_exact_gradient_workspace_test.sh

Run:
  ./run_had_graph_workspace_test.sh
  ./run_exact_gradient_workspace_test.sh
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

EOF
