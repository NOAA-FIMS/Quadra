#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../core/laplace/had_quadra_sparse_exact_hdot_provider.hpp"

DECLARE_ADGRAPH()

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
using Clock = std::chrono::steady_clock;

struct DirectionTiming {
    double ms = 0.0;
    int nnz = 0;
    double norm = 0.0;
};

struct BandedObjective {
    int m;

    template <class T>
    T operator()(const std::vector<T>& x) const {
        const T mu = x[0];
        const T log_sigma = x[1];
        const T log_lambda0 = x[2];
        const T log_lambda_rw = x[3];

        const T inv_sigma2 = exp(T(-2.0) * log_sigma);
        const T lambda0 = exp(log_lambda0);
        const T lambda_rw = exp(log_lambda_rw);

        T nll = T(0.0);

        for (int i = 0; i < m; ++i) {
            const double xd = static_cast<double>(i + 1);
            const T y = T(1.35 + 0.15 * std::sin(0.7 * xd)
                               + 0.03 * std::cos(1.3 * xd));
            const T u = x[4 + i];
            const T resid = y - mu - u;

            nll = nll
                + T(0.5) * resid * resid * inv_sigma2
                + log_sigma
                + T(0.5 * std::log(2.0 * kPi));

            nll = nll + T(0.5) * lambda0 * u * u;
        }

        for (int i = 1; i < m; ++i) {
            const T ui = x[4 + i];
            const T uim1 = x[4 + i - 1];
            const T diff = ui - uim1;
            nll = nll + T(0.5) * lambda_rw * diff * diff;
        }

        return nll;
    }
};

quadra::laplace::RandomHessianPattern tridiagonal_pattern(int m) {
    quadra::laplace::RandomHessianPattern pattern;
    pattern.reserve(static_cast<size_t>(2 * m - 1));

    for (int i = 0; i < m; ++i) {
        pattern.emplace_back(i, i);
        if (i > 0) {
            pattern.emplace_back(i, i - 1);
        }
    }

    return pattern;
}

Eigen::VectorXd make_y(int m) {
    Eigen::VectorXd y(m);
    for (int i = 0; i < m; ++i) {
        const double x = static_cast<double>(i + 1);
        y[i] = 1.35 + 0.15 * std::sin(0.7 * x)
                    + 0.03 * std::cos(1.3 * x);
    }
    return y;
}

Eigen::MatrixXd Huu_dense(const Eigen::VectorXd& theta, int m) {
    const double inv_sigma2 = std::exp(-2.0 * theta[1]);
    const double lambda0 = std::exp(theta[2]);
    const double lambda_rw = std::exp(theta[3]);

    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(m, m);

    for (int i = 0; i < m; ++i) {
        H(i, i) += inv_sigma2 + lambda0;
    }

    for (int i = 1; i < m; ++i) {
        H(i, i) += lambda_rw;
        H(i - 1, i - 1) += lambda_rw;
        H(i, i - 1) -= lambda_rw;
        H(i - 1, i) -= lambda_rw;
    }

    return H;
}

Eigen::VectorXd uhat_for(const Eigen::VectorXd& theta, int m) {
    const double mu = theta[0];
    const double inv_sigma2 = std::exp(-2.0 * theta[1]);
    const Eigen::VectorXd y = make_y(m);
    const Eigen::MatrixXd H = Huu_dense(theta, m);

    Eigen::VectorXd rhs(m);
    for (int i = 0; i < m; ++i) {
        rhs[i] = inv_sigma2 * (y[i] - mu);
    }

    Eigen::LDLT<Eigen::MatrixXd> ldlt(H);
    return ldlt.solve(rhs);
}

template <class Fn>
double mean_ms(Fn&& fn, int reps) {
    const auto start = Clock::now();

    for (int r = 0; r < reps; ++r) {
        volatile double sink = fn();
        (void)sink;
    }

    const auto end = Clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    return elapsed.count() / static_cast<double>(reps);
}

DirectionTiming time_direction(int m, int theta_index, int reps) {
    Eigen::VectorXd theta(4);
    theta << 1.1, std::log(0.45), std::log(0.75), std::log(1.35);

    const Eigen::VectorXd uhat = uhat_for(theta, m);

    BandedObjective objective{m};

    auto provider =
        quadra::laplace::make_had_quadra_sparse_exact_hdot_provider(
            objective,
            4,
            m,
            tridiagonal_pattern(m),
            0.0);

    Eigen::SparseMatrix<double> Hdot0 =
        provider.sparse(theta, uhat, theta_index);

    DirectionTiming out;
    out.nnz = Hdot0.nonZeros();
    out.norm = Hdot0.norm();

    out.ms = mean_ms([&]() {
        Eigen::SparseMatrix<double> Hdot =
            provider.sparse(theta, uhat, theta_index);
        return Hdot.norm() + static_cast<double>(Hdot.nonZeros());
    }, reps);

    return out;
}

}  // namespace

int main(int argc, char** argv) {
    int reps = 25;

    if (argc > 1) {
        reps = std::stoi(argv[1]);
    }

    if (reps <= 0) {
        throw std::invalid_argument("reps must be positive.");
    }

    const std::vector<int> dims = {10, 25, 50, 100, 250};
    const char* names[4] = {"mu", "log_sigma", "log_lambda0", "log_lambda_rw"};

    std::cout << "Per-direction sparse exact Hdot benchmark\n";
    std::cout << "reps per direction = " << reps << "\n\n";

    std::cout << std::setw(8) << "m"
              << std::setw(16) << "direction"
              << std::setw(18) << "Hdot ms"
              << std::setw(12) << "nnz"
              << std::setw(18) << "norm"
              << "\n";

    std::cout << std::scientific << std::setprecision(6);

    for (int m : dims) {
        for (int j = 0; j < 4; ++j) {
            const DirectionTiming r = time_direction(m, j, reps);

            std::cout << std::setw(8) << m
                      << std::setw(16) << names[j]
                      << std::setw(18) << r.ms
                      << std::setw(12) << r.nnz
                      << std::setw(18) << r.norm
                      << "\n";
        }
    }

    std::cout << "\nBenchmark complete.\n";
    return 0;
}
