#include <chrono>
#include <iomanip>
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

template <typename F>
double time_ms(F&& f) {
    const auto start = std::chrono::high_resolution_clock::now();
    f();
    const auto end = std::chrono::high_resolution_clock::now();

    return std::chrono::duration<double, std::milli>(end - start).count();
}

int main() {
    std::cout << "\nQuadra sparse factorization cache benchmark\n\n";

    std::cout
        << std::setw(8)  << "n"
        << std::setw(10) << "evals"
        << std::setw(16) << "compute ms"
        << std::setw(20) << "cached fact ms"
        << std::setw(14) << "speedup"
        << std::setw(14) << "nnz"
        << std::setw(18) << "last logdet"
        << "\n";

    std::cout << std::string(100, '-') << "\n";

    const int n_evals = 200;

    for (int n : std::vector<int>{50, 100, 250, 500, 1000, 2500}) {
        std::vector<Eigen::SparseMatrix<double>> matrices;
        matrices.reserve(n_evals);

        for (int i = 0; i < n_evals; ++i) {
            const double shift =
                0.001 * static_cast<double>(i);
            matrices.push_back(make_tridiagonal_spd(n, shift));
        }

        double last_logdet_compute = 0.0;

        double compute_ms = time_ms([&]() {
            for (const auto& H : matrices) {
                quadra::SparseFactorizationCache cache;
                cache.compute(H);
                last_logdet_compute = cache.logdet();
            }
        });

        double last_logdet_cached = 0.0;

        double cached_ms = time_ms([&]() {
            quadra::SparseFactorizationCache cache;
            cache.analyze_pattern(matrices.front());

            for (const auto& H : matrices) {
                cache.factorize(H);
                last_logdet_cached = cache.logdet();
            }
        });

        const double speedup = compute_ms / cached_ms;

        std::cout
            << std::setw(8)  << n
            << std::setw(10) << n_evals
            << std::setw(16) << std::fixed << std::setprecision(3) << compute_ms
            << std::setw(20) << cached_ms
            << std::setw(14) << speedup
            << std::setw(14) << matrices.back().nonZeros()
            << std::setw(18) << last_logdet_cached
            << "\n";

        if (std::abs(last_logdet_compute - last_logdet_cached) > 1e-8) {
            std::cerr << "logdet mismatch for n = " << n << "\n";
            return 1;
        }
    }

    return 0;
}
