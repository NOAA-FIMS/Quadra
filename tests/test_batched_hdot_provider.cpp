#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../core/laplace/batched_hdot_provider.hpp"
#include "../core/laplace/sparse_huu_factorization.hpp"

DECLARE_ADGRAPH()

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

struct SparseRw1Objective {
    int m;

    template <class T>
    T operator()(const std::vector<T>& x) const {
        const T mu = x[0];
        const T log_sigma = x[1];
        const T log_lambda0 = x[2];
        const T log_lambda_rw = x[3];
        const T log_beta = x[4];

        const T inv_sigma2 = exp(T(-2.0) * log_sigma);
        const T lambda0 = exp(log_lambda0);
        const T lambda_rw = exp(log_lambda_rw);
        const T beta = exp(log_beta);

        T nll = T(0.0);

        for (int i = 0; i < m; ++i) {
            const double xd = static_cast<double>(i + 1);
            const T y = T(0.6 + 0.10 * std::sin(0.21 * xd)
                              + 0.07 * std::cos(0.47 * xd));
            const T u = x[5 + i];
            const T resid = y - mu - u;

            nll = nll
                + T(0.5) * resid * resid * inv_sigma2
                + log_sigma
                + T(0.5 * std::log(2.0 * kPi))
                + T(0.5) * lambda0 * u * u
                + beta * exp(u);
        }

        for (int i = 1; i < m; ++i) {
            const T diff = x[5 + i] - x[5 + i - 1];
            nll = nll + T(0.5) * lambda_rw * diff * diff;
        }

        return nll;
    }
};

quadra::laplace::RandomHessianPattern tridiagonal_pattern(int m) {
    quadra::laplace::RandomHessianPattern pattern;
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
        y[i] = 0.6 + 0.10 * std::sin(0.21 * x)
                   + 0.07 * std::cos(0.47 * x);
    }
    return y;
}

Eigen::SparseMatrix<double> Huu_sparse_direct(const Eigen::VectorXd& theta,
                                              const Eigen::VectorXd& u) {
    const int m = static_cast<int>(u.size());
    const double inv_sigma2 = std::exp(-2.0 * theta[1]);
    const double lambda0 = std::exp(theta[2]);
    const double lambda_rw = std::exp(theta[3]);
    const double beta = std::exp(theta[4]);

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<size_t>(3 * m - 2));

    for (int i = 0; i < m; ++i) {
        double diag = inv_sigma2 + lambda0 + beta * std::exp(u[i]);
        if (i > 0) diag += lambda_rw;
        if (i + 1 < m) diag += lambda_rw;
        triplets.emplace_back(i, i, diag);
    }

    for (int i = 1; i < m; ++i) {
        triplets.emplace_back(i, i - 1, -lambda_rw);
        triplets.emplace_back(i - 1, i, -lambda_rw);
    }

    Eigen::SparseMatrix<double> H(m, m);
    H.setFromTriplets(triplets.begin(), triplets.end());
    H.makeCompressed();
    return H;
}

Eigen::VectorXd random_gradient(const Eigen::VectorXd& theta,
                                const Eigen::VectorXd& u) {
    const int m = static_cast<int>(u.size());
    const double mu = theta[0];
    const double inv_sigma2 = std::exp(-2.0 * theta[1]);
    const double lambda0 = std::exp(theta[2]);
    const double lambda_rw = std::exp(theta[3]);
    const double beta = std::exp(theta[4]);
    const Eigen::VectorXd y = make_y(m);

    Eigen::VectorXd g = Eigen::VectorXd::Zero(m);

    for (int i = 0; i < m; ++i) {
        const double resid = y[i] - mu - u[i];
        g[i] += -resid * inv_sigma2
              + lambda0 * u[i]
              + beta * std::exp(u[i]);
    }

    for (int i = 1; i < m; ++i) {
        const double diff = u[i] - u[i - 1];
        g[i] += lambda_rw * diff;
        g[i - 1] -= lambda_rw * diff;
    }

    return g;
}

Eigen::VectorXd solve_uhat(const Eigen::VectorXd& theta, int m) {
    Eigen::VectorXd u = Eigen::VectorXd::Zero(m);

    for (int iter = 0; iter < 80; ++iter) {
        const Eigen::VectorXd g = random_gradient(theta, u);
        Eigen::LDLT<Eigen::MatrixXd> ldlt(Eigen::MatrixXd(Huu_sparse_direct(theta, u)));
        const Eigen::VectorXd step = ldlt.solve(g);
        u -= step;
        if (step.lpNorm<Eigen::Infinity>() < 1.0e-12) {
            break;
        }
    }

    return u;
}

Eigen::VectorXd f_u_theta_column(const Eigen::VectorXd& theta,
                                 const Eigen::VectorXd& uhat,
                                 int theta_index) {
    const int m = static_cast<int>(uhat.size());
    const double inv_sigma2 = std::exp(-2.0 * theta[1]);
    const double lambda0 = std::exp(theta[2]);
    const double lambda_rw = std::exp(theta[3]);
    const double beta = std::exp(theta[4]);
    const Eigen::VectorXd y = make_y(m);

    Eigen::VectorXd col = Eigen::VectorXd::Zero(m);

    if (theta_index == 0) {
        col.array() = inv_sigma2;
        return col;
    }
    if (theta_index == 1) {
        for (int i = 0; i < m; ++i) {
            const double resid = y[i] - theta[0] - uhat[i];
            col[i] = 2.0 * resid * inv_sigma2;
        }
        return col;
    }
    if (theta_index == 2) {
        for (int i = 0; i < m; ++i) {
            col[i] = lambda0 * uhat[i];
        }
        return col;
    }
    if (theta_index == 3) {
        for (int i = 1; i < m; ++i) {
            const double diff = uhat[i] - uhat[i - 1];
            col[i] += lambda_rw * diff;
            col[i - 1] -= lambda_rw * diff;
        }
        return col;
    }
    if (theta_index == 4) {
        for (int i = 0; i < m; ++i) {
            col[i] = beta * std::exp(uhat[i]);
        }
        return col;
    }

    return col;
}

}  // namespace

int main() {
    const int m = 20;
    Eigen::VectorXd theta(5);
    theta << 0.55, std::log(0.65), std::log(0.55), std::log(0.90), std::log(0.25);

    const Eigen::VectorXd uhat = solve_uhat(theta, m);
    const Eigen::SparseMatrix<double> H = Huu_sparse_direct(theta, uhat);
    quadra::laplace::SparseHuuFactorization factor(H);

    auto direction_provider = [&](int theta_index) -> Eigen::VectorXd {
        return -factor.solve(f_u_theta_column(theta, uhat, theta_index));
    };

    SparseRw1Objective objective{m};

    auto batch =
        quadra::laplace::make_batched_hdot_provider(
            objective,
            direction_provider,
            5,
            m,
            tridiagonal_pattern(m),
            0.0);

    const auto result =
        batch.compute(theta, uhat, std::vector<int>{0, 1, 2, 3, 4});

    auto reference =
        quadra::laplace::make_had_quadra_replay_reuse_lazy_implicit_hdot_provider(
            objective,
            direction_provider,
            5,
            m,
            tridiagonal_pattern(m),
            std::vector<int>{0, 1, 2, 3, 4},
            0.0);

    const auto ref_hdots = reference.compute_all_sparse(theta, uhat);

    if (result.Hdots.size() != 5) {
        throw std::runtime_error("batched Hdot provider returned wrong Hdot count.");
    }

    for (int j = 0; j < 5; ++j) {
        Eigen::SparseMatrix<double> diff =
            result.Hdots[static_cast<size_t>(j)] -
            ref_hdots[static_cast<size_t>(j)];

        const double norm = diff.norm();
        std::cout << "direction " << j
                  << " diff norm = " << norm
                  << " nnz = " << result.direction_nnz[static_cast<size_t>(j)]
                  << "\n";

        if (norm > 1.0e-12) {
            throw std::runtime_error("batched Hdot differs from reference.");
        }
    }

    const auto full_vector =
        batch.compute_all_sparse(theta, uhat, std::vector<int>{1, 3});

    if (full_vector.size() != 5) {
        throw std::runtime_error("compute_all_sparse returned wrong vector length.");
    }

    if (full_vector[0].nonZeros() != 0 ||
        full_vector[2].nonZeros() != 0 ||
        full_vector[4].nonZeros() != 0) {
        throw std::runtime_error("inactive directions should be zero matrices.");
    }

    std::cout << "batched Hdot provider scaffold tests passed\n";
    return 0;
}
