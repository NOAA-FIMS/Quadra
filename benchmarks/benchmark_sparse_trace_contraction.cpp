#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../core/laplace/had_quadra_replay_reuse_sparse_hdot_provider.hpp"
#include "../core/laplace/sparse_trace_contraction.hpp"

DECLARE_ADGRAPH()

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
using Clock = std::chrono::steady_clock;

struct Result {
    double dense_rhs_ms = 0.0;
    double selected_cols_ms = 0.0;
    double cached_cols_ms = 0.0;
    double selected_speedup = 0.0;
    double cached_speedup = 0.0;
    double max_abs_diff = 0.0;
    int needed_cols = 0;
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
                + T(0.5 * std::log(2.0 * kPi))
                + T(0.5) * lambda0 * u * u;
        }

        for (int i = 1; i < m; ++i) {
            const T diff = x[4 + i] - x[4 + i - 1];
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

Eigen::MatrixXd Huu(const Eigen::VectorXd& theta, int m) {
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
    const Eigen::MatrixXd H = Huu(theta, m);

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

Result run_case(int m, int reps) {
    Eigen::VectorXd theta(4);
    theta << 1.1, std::log(0.45), std::log(0.75), std::log(1.35);

    const Eigen::VectorXd uhat = uhat_for(theta, m);
    const Eigen::MatrixXd H = Huu(theta, m);

    BandedObjective objective{m};

    auto provider =
        quadra::laplace::make_had_quadra_replay_reuse_sparse_hdot_provider(
            objective,
            4,
            m,
            tridiagonal_pattern(m),
            std::vector<int>{1, 2, 3},
            0.0);

    const auto Hdots = provider.compute_all_sparse(theta, uhat);

    Eigen::LDLT<Eigen::MatrixXd> ldlt(H);
    if (ldlt.info() != Eigen::Success) {
        throw std::runtime_error("LDLT failed.");
    }

    const auto needed =
        quadra::laplace::needed_columns_from_sparse_matrices(Hdots);

    quadra::laplace::SelectedInverseColumnTraceCache cache(ldlt, m, needed);

    auto dense_all = [&]() {
        double s = 0.0;
        for (const auto& Hdot : Hdots) {
            s += quadra::laplace::trace_hinv_hdot_dense_rhs(ldlt, Hdot);
        }
        return s;
    };

    auto selected_all = [&]() {
        double s = 0.0;
        for (const auto& Hdot : Hdots) {
            s += quadra::laplace::trace_hinv_hdot_selected_inverse_columns(ldlt, Hdot);
        }
        return s;
    };

    auto cached_all = [&]() {
        double s = 0.0;
        for (const auto& Hdot : Hdots) {
            s += cache.trace(Hdot);
        }
        return s;
    };

    const double dense_ref = dense_all();
    const double selected_ref = selected_all();
    const double cached_ref = cached_all();

    Result out;
    out.needed_cols = static_cast<int>(needed.size());
    out.max_abs_diff =
        std::max(std::abs(dense_ref - selected_ref),
                 std::abs(dense_ref - cached_ref));

    out.dense_rhs_ms = mean_ms(dense_all, reps);
    out.selected_cols_ms = mean_ms(selected_all, reps);
    out.cached_cols_ms = mean_ms(cached_all, reps);

    out.selected_speedup = out.dense_rhs_ms / out.selected_cols_ms;
    out.cached_speedup = out.dense_rhs_ms / out.cached_cols_ms;

    return out;
}

}  // namespace

int main(int argc, char** argv) {
    int reps = 20;
    if (argc > 1) {
        reps = std::stoi(argv[1]);
    }

    const std::vector<int> dims = {10, 25, 50, 100, 250};

    std::cout << "Sparse trace contraction benchmark\n";
    std::cout << "reps per case = " << reps << "\n\n";

    std::cout << std::setw(8) << "m"
              << std::setw(14) << "needed cols"
              << std::setw(18) << "dense RHS ms"
              << std::setw(20) << "selected ms"
              << std::setw(18) << "cached ms"
              << std::setw(18) << "selected spd"
              << std::setw(16) << "cached spd"
              << std::setw(16) << "max diff"
              << "\n";

    std::cout << std::scientific << std::setprecision(6);

    for (int m : dims) {
        const Result r = run_case(m, reps);

        std::cout << std::setw(8) << m
                  << std::setw(14) << r.needed_cols
                  << std::setw(18) << r.dense_rhs_ms
                  << std::setw(20) << r.selected_cols_ms
                  << std::setw(18) << r.cached_cols_ms
                  << std::setw(18) << r.selected_speedup
                  << std::setw(16) << r.cached_speedup
                  << std::setw(16) << r.max_abs_diff
                  << "\n";
    }

    std::cout << "\nBenchmark complete.\n";
    return 0;
}
