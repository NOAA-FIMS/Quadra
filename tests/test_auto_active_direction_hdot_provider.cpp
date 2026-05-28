#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../core/laplace/active_direction_hdot_provider.hpp"
#include "../core/laplace/auto_active_direction_hdot_provider.hpp"
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

Eigen::VectorXd make_y(int m) {
    Eigen::VectorXd y(m);
    for (int i = 0; i < m; ++i) {
        const double x = static_cast<double>(i + 1);
        y[i] = 1.35 + 0.15 * std::sin(0.7 * x)
                    + 0.03 * std::cos(1.3 * x);
    }
    return y;
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

Eigen::VectorXd uhat_for(const Eigen::VectorXd& theta, int m) {
    const double mu = theta[0];
    const double inv_sigma2 = std::exp(-2.0 * theta[1]);
    const Eigen::VectorXd y = make_y(m);
    const Eigen::MatrixXd H = Huu(theta, Eigen::VectorXd::Zero(m));

    Eigen::VectorXd rhs(m);
    for (int i = 0; i < m; ++i) {
        rhs[i] = inv_sigma2 * (y[i] - mu);
    }

    Eigen::LDLT<Eigen::MatrixXd> ldlt(H);
    return ldlt.solve(rhs);
}

Eigen::VectorXd joint_grad(const Eigen::VectorXd& theta,
                           const Eigen::VectorXd& uhat) {
    const int m = static_cast<int>(uhat.size());
    const double mu = theta[0];
    const double inv_sigma2 = std::exp(-2.0 * theta[1]);
    const double lambda0 = std::exp(theta[2]);
    const double lambda_rw = std::exp(theta[3]);
    const Eigen::VectorXd y = make_y(m);

    Eigen::VectorXd g(4);
    g.setZero();

    for (int i = 0; i < m; ++i) {
        const double resid = y[i] - mu - uhat[i];
        g[0] += -resid * inv_sigma2;
        g[1] += 1.0 - resid * resid * inv_sigma2;
        g[2] += 0.5 * lambda0 * uhat[i] * uhat[i];
    }

    for (int i = 1; i < m; ++i) {
        const double diff = uhat[i] - uhat[i - 1];
        g[3] += 0.5 * lambda_rw * diff * diff;
    }

    return g;
}

void run_test() {
    const int m = 30;

    Eigen::VectorXd theta(4);
    theta << 1.1, std::log(0.45), std::log(0.75), std::log(1.35);

    const Eigen::VectorXd uhat = uhat_for(theta, m);
    const Eigen::VectorXd gj = joint_grad(theta, uhat);

    BandedObjective objective{m};

    auto base =
        quadra::laplace::make_had_quadra_sparse_exact_hdot_provider(
            objective,
            4,
            m,
            tridiagonal_pattern(m),
            0.0);

    auto manual =
        quadra::laplace::make_active_direction_hdot_provider(
            base,
            m,
            std::vector<int>{1, 2, 3});

    auto auto_active =
        quadra::laplace::make_auto_active_direction_hdot_provider(
            base,
            theta,
            uhat,
            4,
            m,
            0.0);

    const auto& discovery = auto_active.discovery();

    if (discovery.active_directions != std::vector<int>({1, 2, 3})) {
        throw std::runtime_error(
            "auto active-direction discovery did not find {1,2,3}.");
    }

    if (discovery.hdot_nonzeros[0] != 0 ||
        discovery.hdot_norms[0] != 0.0) {
        throw std::runtime_error(
            "mu direction should have zero Hdot discovery metrics.");
    }

    auto Hfn = [](const Eigen::VectorXd& th, const Eigen::VectorXd& uh) {
        return Huu(th, uh);
    };

    const Eigen::VectorXd g_manual =
        quadra::laplace::full_exact_laplace_gradient_hdot(
            gj, Hfn, manual, theta, uhat);

    const Eigen::VectorXd g_auto =
        quadra::laplace::full_exact_laplace_gradient_hdot(
            gj, Hfn, auto_active, theta, uhat);

    const double max_abs = (g_manual - g_auto).cwiseAbs().maxCoeff();
    if (max_abs > 1.0e-12) {
        throw std::runtime_error(
            "auto active-direction gradient differs from manual active gradient.");
    }
}

}  // namespace

int main() {
    run_test();
    std::cout << "auto active-direction Hdot provider tests passed\n";
    return 0;
}
