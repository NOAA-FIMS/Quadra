#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../core/laplace/exact_laplace_gradient_adapter.hpp"

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
    if (ldlt.info() != Eigen::Success) {
        throw std::runtime_error("LDLT failed in uhat_for.");
    }

    return ldlt.solve(rhs);
}

double joint_objective_double(const Eigen::VectorXd& theta,
                              const Eigen::VectorXd& uhat) {
    const int m = static_cast<int>(uhat.size());
    BandedObjective objective{m};

    std::vector<double> x(static_cast<size_t>(4 + m));
    for (int j = 0; j < 4; ++j) {
        x[static_cast<size_t>(j)] = theta[j];
    }
    for (int i = 0; i < m; ++i) {
        x[static_cast<size_t>(4 + i)] = uhat[i];
    }

    return objective(x);
}

double laplace_objective(const Eigen::VectorXd& theta, int m) {
    const Eigen::VectorXd uhat = uhat_for(theta, m);
    const double joint = joint_objective_double(theta, uhat);
    const Eigen::MatrixXd H = Huu(theta, uhat);

    Eigen::LDLT<Eigen::MatrixXd> ldlt(H);
    if (ldlt.info() != Eigen::Success) {
        throw std::runtime_error("LDLT failed in laplace_objective.");
    }

    const double logdet = ldlt.vectorD().array().log().sum();
    return joint + 0.5 * logdet;
}

Eigen::VectorXd joint_envelope_gradient(const Eigen::VectorXd& theta,
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

Eigen::VectorXd central_fd_gradient(const Eigen::VectorXd& theta,
                                    int m,
                                    double step = 1.0e-5) {
    Eigen::VectorXd g(theta.size());

    for (int j = 0; j < theta.size(); ++j) {
        const double h = step * std::max(1.0, std::abs(theta[j]));
        Eigen::VectorXd plus = theta;
        Eigen::VectorXd minus = theta;

        plus[j] += h;
        minus[j] -= h;

        g[j] = (laplace_objective(plus, m) -
                laplace_objective(minus, m)) / (2.0 * h);
    }

    return g;
}

void run_case(int m, const Eigen::VectorXd& theta) {
    const Eigen::VectorXd uhat = uhat_for(theta, m);
    const Eigen::VectorXd gj = joint_envelope_gradient(theta, uhat);

    BandedObjective objective{m};

    auto Hfn = [](const Eigen::VectorXd& th, const Eigen::VectorXd& uh) {
        return Huu(th, uh);
    };

    auto adapter =
        quadra::laplace::make_exact_laplace_gradient_adapter(
            objective,
            Hfn,
            4,
            m,
            tridiagonal_pattern(m),
            theta,
            uhat);

    quadra::laplace::ExactLaplaceGradientInputs inputs;
    inputs.theta = theta;
    inputs.uhat = uhat;
    inputs.joint_envelope_gradient = gj;

    const Eigen::VectorXd g_exact = adapter.gradient(inputs);
    const Eigen::VectorXd g_fd = central_fd_gradient(theta, m);

    const Eigen::VectorXd abs_diff = (g_exact - g_fd).cwiseAbs();
    const double max_abs = abs_diff.maxCoeff();

    std::cout << "\n=== adapter FD validation: m = " << m << " ===\n";
    std::cout << std::scientific << std::setprecision(10);
    std::cout << "exact = " << g_exact.transpose() << "\n";
    std::cout << "fd    = " << g_fd.transpose() << "\n";
    std::cout << "diff  = " << abs_diff.transpose() << "\n";
    std::cout << "max abs diff = " << max_abs << "\n";

    if (max_abs > 2.0e-7) {
        throw std::runtime_error(
            "adapter exact gradient failed FD Laplace objective validation.");
    }
}

}  // namespace

int main() {
    Eigen::VectorXd theta(4);

    theta << 1.1, std::log(0.45), std::log(0.75), std::log(1.35);
    run_case(10, theta);
    run_case(25, theta);
    run_case(50, theta);

    theta << 1.25, std::log(0.70), std::log(0.55), std::log(0.95);
    run_case(20, theta);

    std::cout << "\nexact Laplace gradient adapter FD validation passed\n";
    return 0;
}
