#include <cmath>
#include <iostream>
#include <vector>

#include "../core/autodiff.hpp"

DECLARE_ADGRAPH();

int main()
{
    std::cout << "had_quadra production forward replay test\n";

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

    quadra::set_value(x, 3.0);

    scope.forward();
    scope.zero_adjoints();
    scope.backward(y);

    double y2 = quadra::value_of(y);

    Eigen::VectorXd g2 =
        quadra::extract_gradient(std::vector<quadra::AD>{x});

    std::cout << "replayed y = " << y2 << "\n";
    std::cout << "replayed grad = " << g2[0] << "\n";

    if (std::abs(y1 - 5.0) > 1e-12)
    {
        std::cerr << "FAIL: initial value mismatch\n";
        return 1;
    }

    if (std::abs(g1[0] - 4.0) > 1e-12)
    {
        std::cerr << "FAIL: initial gradient mismatch\n";
        return 1;
    }

    if (std::abs(y2 - 10.0) > 1e-12)
    {
        std::cerr << "FAIL: replayed value mismatch\n";
        return 1;
    }

    if (std::abs(g2[0] - 6.0) > 1e-12)
    {
        std::cerr << "FAIL: replayed gradient mismatch\n";
        return 1;
    }

    std::cout << "PASS\n";
    return 0;
}
