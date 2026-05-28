#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../core/laplace/full_exact_laplace_gradient_hdot.hpp"
#include "../core/laplace/had_quadra_exact_edge_hdot_provider.hpp"
#include "../core/laplace/had_quadra_sparse_exact_hdot_provider.hpp"

DECLARE_ADGRAPH()

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

struct GaussianCombinedObjective {
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
                + T(0.5 * std::log(2.0 * kPi));

            nll = nll
                + T(0.5) * u * u / tau2
                + log_tau
                + T(0.5 * std::log(2.0 * kPi));
        }

        return nll;
    }
};

Eigen::VectorXd optimize_random_effects(const Eigen::VectorXd& theta, int m) {
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

Eigen::VectorXd joint_envelope_gradient(const Eigen::VectorXd& theta,
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

Eigen::MatrixXd hessian_uu(const Eigen::VectorXd& theta,
                           const Eigen::VectorXd& uhat) {
    const double sigma2 = std::exp(2.0 * theta[1]);
    const double tau2 = std::exp(2.0 * theta[2]);
    const double h = 1.0 / sigma2 + 1.0 / tau2;
    return h * Eigen::MatrixXd::Identity(uhat.size(), uhat.size());
}

double joint_nll_double(const Eigen::VectorXd& theta,
                        const Eigen::VectorXd& uhat) {
    GaussianCombinedObjective obj{static_cast<int>(uhat.size())};

    std::vector<double> x(static_cast<size_t>(3 + uhat.size()));
    x[0] = theta[0];
    x[1] = theta[1];
    x[2] = theta[2];

    for (int i = 0; i < uhat.size(); ++i) {
        x[static_cast<size_t>(3 + i)] = uhat[i];
    }

    return obj(x);
}

double laplace_objective(const Eigen::VectorXd& theta, int m) {
    const Eigen::VectorXd uhat = optimize_random_effects(theta, m);
    const double joint = joint_nll_double(theta, uhat);
    const Eigen::MatrixXd H = hessian_uu(theta, uhat);
    return joint + 0.5 * std::log(H.determinant());
}

Eigen::VectorXd central_fd_gradient(
    const std::function<double(const Eigen::VectorXd&)>& f,
    const Eigen::VectorXd& theta,
    double step = 1.0e-5) {
    Eigen::VectorXd g(theta.size());
    for (int j = 0; j < theta.size(); ++j) {
        const double h = step * std::max(1.0, std::abs(theta[j]));
        Eigen::VectorXd plus = theta;
        Eigen::VectorXd minus = theta;
        plus[j] += h;
        minus[j] -= h;
        g[j] = (f(plus) - f(minus)) / (2.0 * h);
    }
    return g;
}

void run_case(int m, const Eigen::VectorXd& theta) {
    const Eigen::VectorXd uhat = optimize_random_effects(theta, m);

    GaussianCombinedObjective objective{m};

    auto dense_exact =
        quadra::laplace::make_had_quadra_exact_edge_hdot_provider(
            objective, 3, m, 0.0);

    auto sparse_exact =
        quadra::laplace::make_had_quadra_sparse_exact_hdot_provider(
            objective,
            3,
            m,
            quadra::laplace::diagonal_random_hessian_pattern(m),
            0.0);

    std::cout << "\n=== sparse exact Hdot provider: m = "
              << m << " ===\n";

    for (int j = 0; j < 3; ++j) {
        const Eigen::MatrixXd dense = dense_exact(theta, uhat, j);
        const Eigen::SparseMatrix<double> sparse =
            sparse_exact.sparse(theta, uhat, j);
        const Eigen::MatrixXd sparse_dense = Eigen::MatrixXd(sparse);

        const double max_abs = (dense - sparse_dense).cwiseAbs().maxCoeff();

        std::cout << "theta_index " << j
                  << " max |dense exact - sparse exact| = "
                  << std::scientific << std::setprecision(10)
                  << max_abs
                  << " ; sparse nnz = "
                  << sparse.nonZeros() << "\n";

        if (max_abs > 1.0e-12) {
            throw std::runtime_error(
                "sparse exact Hdot does not match dense exact Hdot.");
        }
    }

    const Eigen::VectorXd grad_joint =
        joint_envelope_gradient(theta, uhat);

    auto Hfn = [](const Eigen::VectorXd& th,
                  const Eigen::VectorXd& uh) {
        return hessian_uu(th, uh);
    };

    const Eigen::VectorXd grad_sparse =
        quadra::laplace::full_exact_laplace_gradient_hdot(
            grad_joint, Hfn, sparse_exact, theta, uhat);

    const Eigen::VectorXd grad_dense =
        quadra::laplace::full_exact_laplace_gradient_hdot(
            grad_joint, Hfn, dense_exact, theta, uhat);

    const Eigen::VectorXd grad_fd_obj =
        central_fd_gradient(
            [&](const Eigen::VectorXd& th) {
                return laplace_objective(th, m);
            },
            theta,
            1.0e-5);

    const double diff_dense = (grad_sparse - grad_dense).cwiseAbs().maxCoeff();
    const double diff_fd = (grad_sparse - grad_fd_obj).cwiseAbs().maxCoeff();

    std::cout << "gradient max |sparse exact - dense exact| = "
              << diff_dense << "\n";
    std::cout << "gradient max |sparse exact - FD objective| = "
              << diff_fd << "\n";

    if (diff_dense > 1.0e-12) {
        throw std::runtime_error(
            "sparse exact gradient does not match dense exact gradient.");
    }

    if (diff_fd > 1.0e-5) {
        throw std::runtime_error(
            "sparse exact gradient does not match FD objective.");
    }
}

}  // namespace

int main() {
    Eigen::VectorXd theta(3);

    theta << 1.1, std::log(0.45), std::log(0.8);
    run_case(5, theta);
    run_case(10, theta);
    run_case(25, theta);

    theta << 1.25, std::log(0.7), std::log(0.35);
    run_case(8, theta);

    std::cout << "\nsparse exact Hdot-provider validation passed\n";
    return 0;
}
