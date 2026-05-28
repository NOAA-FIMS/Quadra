#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../core/laplace/sparse_trace_contraction.hpp"

namespace {

Eigen::MatrixXd make_spd_tridiagonal(int n) {
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(n, n);
    for (int i = 0; i < n; ++i) {
        H(i, i) = 4.0;
        if (i > 0) {
            H(i, i - 1) = -1.0;
            H(i - 1, i) = -1.0;
        }
    }
    return H;
}

Eigen::SparseMatrix<double> make_sparse_hdot(int n) {
    std::vector<Eigen::Triplet<double>> triplets;

    for (int i = 0; i < n; ++i) {
        triplets.emplace_back(i, i, 0.5 + 0.1 * i);
        if (i > 0) {
            triplets.emplace_back(i, i - 1, -0.25);
            triplets.emplace_back(i - 1, i, -0.25);
        }
    }

    Eigen::SparseMatrix<double> Hdot(n, n);
    Hdot.setFromTriplets(triplets.begin(), triplets.end());
    Hdot.makeCompressed();
    return Hdot;
}

void run_test() {
    for (int n : {5, 10, 25}) {
        const Eigen::MatrixXd H = make_spd_tridiagonal(n);
        const Eigen::SparseMatrix<double> Hdot = make_sparse_hdot(n);

        Eigen::LDLT<Eigen::MatrixXd> ldlt(H);
        if (ldlt.info() != Eigen::Success) {
            throw std::runtime_error("LDLT failed.");
        }

        const double dense =
            quadra::laplace::trace_hinv_hdot_dense_rhs(ldlt, Hdot);

        const double selected =
            quadra::laplace::trace_hinv_hdot_selected_inverse_columns(ldlt, Hdot);

        std::vector<Eigen::SparseMatrix<double>> matrices{Hdot};
        auto cols =
            quadra::laplace::needed_columns_from_sparse_matrices(matrices);

        quadra::laplace::SelectedInverseColumnTraceCache cache(ldlt, n, cols);
        const double cached = cache.trace(Hdot);

        if (std::abs(dense - selected) > 1.0e-10) {
            std::cerr << "n=" << n << " dense=" << dense
                      << " selected=" << selected << "\n";
            throw std::runtime_error("selected trace mismatch.");
        }

        if (std::abs(dense - cached) > 1.0e-10) {
            std::cerr << "n=" << n << " dense=" << dense
                      << " cached=" << cached << "\n";
            throw std::runtime_error("cached trace mismatch.");
        }
    }
}

}  // namespace

int main() {
    run_test();
    std::cout << "sparse trace contraction tests passed\n";
    return 0;
}
