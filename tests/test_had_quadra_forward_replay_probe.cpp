#include <cmath>
#include <iostream>
#include <vector>

#include "../core/autodiff.hpp"

DECLARE_ADGRAPH();

// This is an EXPECTED-FAIL diagnostic.
//
// It documents the behavior we eventually want from reusable tapes:
//
//   mutate independent values
//   replay forward pass
//   rerun reverse pass
//
// If this unexpectedly passes, then had_quadra already has enough replay
// machinery hidden internally.

int main() {

    std::cout << "had_quadra forward replay probe\n";

    quadra::TapeContext tape;
    quadra::ADScope scope(tape.graph);

    quadra::AD x = 2.0;
    quadra::AD y = x * x + 1.0;

    scope.backward(y);

    double y1 = quadra::value_of(y);

    Eigen::VectorXd g1 =
        quadra::extract_gradient(std::vector<quadra::AD>{x});

    std::cout << "initial y = " << y1 << "\n";
    std::cout << "initial grad = " << g1[0] << "\n";

    // ------------------------------------------------------------
    // mutate independent value
    // ------------------------------------------------------------

    x.val = 3.0;

    // ------------------------------------------------------------
    // zero adjoints
    // ------------------------------------------------------------

    scope.zero_adjoints();

    // ------------------------------------------------------------
    // NO FORWARD REPLAY EXISTS YET
    // ------------------------------------------------------------

    scope.backward(y);

    double y2 = quadra::value_of(y);

    Eigen::VectorXd g2 =
        quadra::extract_gradient(std::vector<quadra::AD>{x});

    std::cout << "mutated y = " << y2 << "\n";
    std::cout << "mutated grad = " << g2[0] << "\n";

    const bool value_ok =
        std::abs(y2 - 10.0) <= 1e-12;

    const bool grad_ok =
        std::abs(g2[0] - 6.0) <= 1e-12;

    std::cout << "forward replay value correct = "
              << value_ok << "\n";

    std::cout << "forward replay gradient correct = "
              << grad_ok << "\n";

    if (!value_ok || !grad_ok) {

        std::cout << "\n";
        std::cout << "EXPECTED RESULT:\n";
        std::cout << "Forward replay does not yet exist.\n";
        std::cout << "Graph primal values remain stale after x.val mutation.\n";
        std::cout << "\n";
        std::cout << "Next architectural requirement:\n";
        std::cout << "  replayable graph operations/opcodes\n";
        std::cout << "  or explicit forward recomputation support.\n";

        return 0;
    }

    std::cout << "\n";
    std::cout << "UNEXPECTED SUCCESS:\n";
    std::cout << "had_quadra appears to support forward replay already.\n";

    return 0;
}
