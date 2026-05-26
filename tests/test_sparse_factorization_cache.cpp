#include <cmath>
#include <iostream>
#include <vector>

#include "../core/laplace/sparse_factorization_cache.hpp"

Eigen::SparseMatrix<double> make_tridiagonal_spd(
    int n,
    double diag_shift
) {
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<size_t>(3 * n - 2));

    for (int i = 0; i < n; ++i) {
        triplets.emplace_back(i, i, 2.0 + diag_shift);

        if (i > 0) {
            triplets.emplace_back(i, i - 1, -0.25);
        }

        if (i + 1 < n) {
            triplets.emplace_back(i, i + 1, -0.25);
        }
    }

    Eigen::SparseMatrix<double> H(n, n);
    H.setFromTriplets(triplets.begin(), triplets.end());
    H.makeCompressed();

    return H;
}

int main() {
    std::cout << "Testing SparseFactorizationCache\n";

    const int n = 10;

    auto H1 = make_tridiagonal_spd(n, 0.0);
    auto H2 = make_tridiagonal_spd(n, 0.5);

    quadra::SparseFactorizationCache cache;

    cache.analyze_pattern(H1);
    cache.factorize(H1);

    Eigen::VectorXd b = Eigen::VectorXd::Ones(n);
    Eigen::VectorXd x1 = cache.solve(b);

    const double ld1 = cache.logdet();

    cache.factorize(H2);

    Eigen::VectorXd x2 = cache.solve(b);

    const double ld2 = cache.logdet();

    std::cout << "analyzed = " << cache.analyzed() << "\n";
    std::cout << "factorized = " << cache.factorized() << "\n";
    std::cout << "nnz = " << cache.nonzeros() << "\n";
    std::cout << "logdet H1 = " << ld1 << "\n";
    std::cout << "logdet H2 = " << ld2 << "\n";
    std::cout << "x1 norm = " << x1.norm() << "\n";
    std::cout << "x2 norm = " << x2.norm() << "\n";

    if (!cache.analyzed() || !cache.factorized()) {
        std::cerr << "FAIL: cache state invalid\n";
        return 1;
    }

    if (cache.nonzeros() != 3 * n - 2) {
        std::cerr << "FAIL: nnz mismatch\n";
        return 1;
    }

    if (!(ld2 > ld1)) {
        std::cerr << "FAIL: expected larger logdet after diagonal shift\n";
        return 1;
    }

    if (!std::isfinite(x1.norm()) || !std::isfinite(x2.norm())) {
        std::cerr << "FAIL: solve produced nonfinite values\n";
        return 1;
    }

    std::cout << "PASS\n";
    return 0;
}
