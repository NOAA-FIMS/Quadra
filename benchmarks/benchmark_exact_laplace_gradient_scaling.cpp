#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../core/laplace/full_exact_laplace_gradient_fd.hpp"
#include "../core/laplace/full_exact_laplace_gradient_hdot.hpp"
#include "../core/laplace/had_quadra_exact_edge_hdot_provider.hpp"
#include "../core/laplace/had_quadra_sparse_exact_hdot_provider.hpp"

DECLARE_ADGRAPH()

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
using Clock = std::chrono::steady_clock;

struct Result {
    double fd_ms = 0.0;
    double dense_ms = 0.0;
    double sparse_ms = 0.0;
    double dense_sparse_diff = 0.0;
    double fd_sparse_diff = 0.0;
    int pattern_nnz = 0;
};

struct GaussianObjective {
    int m;

    template <class T>
    T operator()(const std::vector<T>& x) const {
        const T mu = x[0];
        const T log_sigma = x[1];
        const T log_tau = x[2];
        const T sigma2 = exp(T(2.0) * log_sigma);
        const T tau2 = exp(T(2.0) * log_tau);
        T nll = T(0.0);

        for (int i = 0; i < m; ++i) {
            const double xd = static_cast<double>(i + 1);
            const T y = T(1.35 + 0.15 * std::sin(0.7 * xd)
                               + 0.03 * std::cos(1.3 * xd));
            const T u = x[3 + i];
            const T resid = y - mu - u;

            nll = nll
                + T(0.5) * resid * resid / sigma2
                + log_sigma
                + T(0.5 * std::log(2.0 * kPi))
                + T(0.5) * u * u / tau2
                + log_tau
                + T(0.5 * std::log(2.0 * kPi));
        }
        return nll;
    }
};

Eigen::VectorXd uhat_for(const Eigen::VectorXd& theta, int m) {
    const double mu = theta[0];
    const double sigma2 = std::exp(2.0 * theta[1]);
    const double tau2 = std::exp(2.0 * theta[2]);
    const double shrink = tau2 / (sigma2 + tau2);

    Eigen::VectorXd uhat(m);
    for (int i = 0; i < m; ++i) {
        const double x = static_cast<double>(i + 1);
        const double y = 1.35 + 0.15 * std::sin(0.7 * x)
                              + 0.03 * std::cos(1.3 * x);
        uhat[i] = shrink * (y - mu);
    }
    return uhat;
}

Eigen::VectorXd joint_grad(const Eigen::VectorXd& theta,
                           const Eigen::VectorXd& uhat) {
    const double mu = theta[0];
    const double sigma2 = std::exp(2.0 * theta[1]);
    const double tau2 = std::exp(2.0 * theta[2]);

    Eigen::VectorXd g(3);
    g.setZero();

    for (int i = 0; i < uhat.size(); ++i) {
        const double x = static_cast<double>(i + 1);
        const double y = 1.35 + 0.15 * std::sin(0.7 * x)
                              + 0.03 * std::cos(1.3 * x);
        const double resid = y - mu - uhat[i];

        g[0] += -resid / sigma2;
        g[1] += 1.0 - resid * resid / sigma2;
        g[2] += 1.0 - uhat[i] * uhat[i] / tau2;
    }
    return g;
}

Eigen::MatrixXd Huu(const Eigen::VectorXd& theta,
                    const Eigen::VectorXd& uhat) {
    const double sigma2 = std::exp(2.0 * theta[1]);
    const double tau2 = std::exp(2.0 * theta[2]);
    const double h = 1.0 / sigma2 + 1.0 / tau2;
    return h * Eigen::MatrixXd::Identity(uhat.size(), uhat.size());
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
    Eigen::VectorXd theta(3);
    theta << 1.1, std::log(0.45), std::log(0.8);

    const Eigen::VectorXd uhat = uhat_for(theta, m);
    const Eigen::VectorXd gj = joint_grad(theta, uhat);
    GaussianObjective objective{m};

    auto HthetaOnly = [&](const Eigen::VectorXd& th) {
        return Huu(th, uhat);
    };

    auto Hfn = [](const Eigen::VectorXd& th, const Eigen::VectorXd& uh) {
        return Huu(th, uh);
    };

    auto dense_provider =
        quadra::laplace::make_had_quadra_exact_edge_hdot_provider(
            objective, 3, m, 0.0);

    auto sparse_provider =
        quadra::laplace::make_had_quadra_sparse_exact_hdot_provider(
            objective,
            3,
            m,
            quadra::laplace::diagonal_random_hessian_pattern(m),
            0.0);

    quadra::laplace::FullExactLaplaceGradientFDOptions fd_options;
    fd_options.step = 1.0e-6;
    fd_options.relative_step = true;

    const Eigen::VectorXd g_fd =
        quadra::laplace::full_exact_laplace_gradient_fd(
            gj, HthetaOnly, theta, fd_options);

    const Eigen::VectorXd g_dense =
        quadra::laplace::full_exact_laplace_gradient_hdot(
            gj, Hfn, dense_provider, theta, uhat);

    const Eigen::VectorXd g_sparse =
        quadra::laplace::full_exact_laplace_gradient_hdot(
            gj, Hfn, sparse_provider, theta, uhat);

    Result out;
    out.pattern_nnz = m;
    out.dense_sparse_diff = (g_dense - g_sparse).cwiseAbs().maxCoeff();
    out.fd_sparse_diff = (g_fd - g_sparse).cwiseAbs().maxCoeff();

    out.fd_ms = mean_ms([&]() {
        return quadra::laplace::full_exact_laplace_gradient_fd(
            gj, HthetaOnly, theta, fd_options).sum();
    }, reps);

    out.dense_ms = mean_ms([&]() {
        return quadra::laplace::full_exact_laplace_gradient_hdot(
            gj, Hfn, dense_provider, theta, uhat).sum();
    }, reps);

    out.sparse_ms = mean_ms([&]() {
        return quadra::laplace::full_exact_laplace_gradient_hdot(
            gj, Hfn, sparse_provider, theta, uhat).sum();
    }, reps);

    return out;
}

}  // namespace

int main(int argc, char** argv) {
    int reps = 20;
    if (argc > 1) {
        reps = std::stoi(argv[1]);
    }
    if (reps <= 0) {
        throw std::invalid_argument("reps must be positive.");
    }

    const std::vector<int> dims = {5, 10, 25, 50, 100, 250};

    std::cout << "Exact Laplace gradient scaling benchmark\n";
    std::cout << "reps per case = " << reps << "\n\n";

    std::cout << std::setw(8) << "m"
              << std::setw(12) << "pattern"
              << std::setw(16) << "FD ms"
              << std::setw(18) << "dense ms"
              << std::setw(18) << "sparse ms"
              << std::setw(16) << "dense/sparse"
              << std::setw(16) << "FD/sparse"
              << std::setw(18) << "|dense-sparse|"
              << std::setw(18) << "|FD-sparse|"
              << "\n";

    std::cout << std::scientific << std::setprecision(6);

    for (int m : dims) {
        const Result r = run_case(m, reps);
        const double dense_ratio = r.dense_ms / r.sparse_ms;
        const double fd_ratio = r.fd_ms / r.sparse_ms;

        std::cout << std::setw(8) << m
                  << std::setw(12) << r.pattern_nnz
                  << std::setw(16) << r.fd_ms
                  << std::setw(18) << r.dense_ms
                  << std::setw(18) << r.sparse_ms
                  << std::setw(16) << dense_ratio
                  << std::setw(16) << fd_ratio
                  << std::setw(18) << r.dense_sparse_diff
                  << std::setw(18) << r.fd_sparse_diff
                  << "\n";
    }

    std::cout << "\nBenchmark complete.\n";
    return 0;
}
