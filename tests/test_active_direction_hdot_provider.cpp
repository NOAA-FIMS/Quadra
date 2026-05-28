#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../core/laplace/active_direction_hdot_provider.hpp"
#include "../core/laplace/full_exact_laplace_gradient_hdot.hpp"
#include "../core/laplace/had_quadra_sparse_exact_hdot_provider.hpp"

DECLARE_ADGRAPH()

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

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
            const T diff = x[4 + i] - x[4 + i - 1];
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

Eigen::MatrixXd Huu(const Eigen::VectorXd& theta, const Eigen::VectorXd& uhat) {
    const int m = static_cast<int>(uhat.size());
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

void run_test() {
    const int m = 20;

    Eigen::VectorXd theta(4);
    theta << 1.1, std::log(0.45), std::log(0.75), std::log(1.35);

    Eigen::VectorXd uhat = Eigen::VectorXd::Zero(m);

    BandedObjective objective{m};

    auto base =
        quadra::laplace::make_had_quadra_sparse_exact_hdot_provider(
            objective,
            4,
            m,
            tridiagonal_pattern(m),
            0.0);

    auto active =
        quadra::laplace::make_active_direction_hdot_provider(
            base,
            m,
            std::vector<int>{1, 2, 3});

    for (int j = 0; j < 4; ++j) {
        const Eigen::SparseMatrix<double> Hbase = base.sparse(theta, uhat, j);
        const Eigen::SparseMatrix<double> Hactive = active.sparse(theta, uhat, j);

        if (j == 0) {
            if (Hactive.nonZeros() != 0) {
                throw std::runtime_error("inactive direction should return zero Hdot.");
            }
        } else {
            const Eigen::MatrixXd diff =
                Eigen::MatrixXd(Hbase) - Eigen::MatrixXd(Hactive);

            if (diff.cwiseAbs().maxCoeff() > 1.0e-12) {
                throw std::runtime_error("active direction wrapper changed Hdot.");
            }
        }
    }
}

}  // namespace

int main() {
    run_test();
    std::cout << "active-direction Hdot provider tests passed\n";
    return 0;
}
