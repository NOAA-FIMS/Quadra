#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../core/laplace/had_quadra_replay_reuse_implicit_hdot_provider.hpp"

DECLARE_ADGRAPH()

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

// Non-quadratic random-effect objective.
//
// The beta * exp(u_i) term makes H_uu depend on uhat, so the implicit
// du*/dtheta contribution to dH_uu/dtheta is required.
struct NonlinearDiagObjective {
    int m;

    template <class T>
    T operator()(const std::vector<T>& x) const {
        const T mu = x[0];
        const T log_sigma = x[1];
        const T log_lambda0 = x[2];
        const T log_beta = x[3];

        const T inv_sigma2 = exp(T(-2.0) * log_sigma);
        const T lambda0 = exp(log_lambda0);
        const T beta = exp(log_beta);

        T nll = T(0.0);

        for (int i = 0; i < m; ++i) {
            const double xd = static_cast<double>(i + 1);
            const T y = T(0.8 + 0.08 * std::sin(0.4 * xd)
                              + 0.05 * std::cos(0.9 * xd));
            const T u = x[4 + i];
            const T resid = y - mu - u;

            nll = nll
                + T(0.5) * resid * resid * inv_sigma2
                + log_sigma
                + T(0.5 * std::log(2.0 * kPi))
                + T(0.5) * lambda0 * u * u
                + beta * exp(u);
        }

        return nll;
    }
};

quadra::laplace::RandomHessianPattern diagonal_pattern(int m) {
    quadra::laplace::RandomHessianPattern pattern;
    pattern.reserve(static_cast<size_t>(m));
    for (int i = 0; i < m; ++i) {
        pattern.emplace_back(i, i);
    }
    return pattern;
}

Eigen::VectorXd make_y(int m) {
    Eigen::VectorXd y(m);
    for (int i = 0; i < m; ++i) {
        const double x = static_cast<double>(i + 1);
        y[i] = 0.8 + 0.08 * std::sin(0.4 * x)
                   + 0.05 * std::cos(0.9 * x);
    }
    return y;
}

Eigen::VectorXd solve_uhat(const Eigen::VectorXd& theta, int m) {
    const double mu = theta[0];
    const double inv_sigma2 = std::exp(-2.0 * theta[1]);
    const double lambda0 = std::exp(theta[2]);
    const double beta = std::exp(theta[3]);
    const Eigen::VectorXd y = make_y(m);

    Eigen::VectorXd u = Eigen::VectorXd::Zero(m);

    for (int iter = 0; iter < 80; ++iter) {
        double max_step = 0.0;

        for (int i = 0; i < m; ++i) {
            const double resid = y[i] - mu - u[i];
            const double expu = std::exp(u[i]);

            const double g =
                -resid * inv_sigma2
                + lambda0 * u[i]
                + beta * expu;

            const double H =
                inv_sigma2 + lambda0 + beta * expu;

            const double step = g / H;
            u[i] -= step;
            max_step = std::max(max_step, std::abs(step));
        }

        if (max_step < 1.0e-13) {
            break;
        }
    }

    return u;
}

Eigen::MatrixXd Huu(const Eigen::VectorXd& theta,
                    const Eigen::VectorXd& uhat) {
    const int m = static_cast<int>(uhat.size());
    const double inv_sigma2 = std::exp(-2.0 * theta[1]);
    const double lambda0 = std::exp(theta[2]);
    const double beta = std::exp(theta[3]);

    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(m, m);

    for (int i = 0; i < m; ++i) {
        H(i, i) = inv_sigma2 + lambda0 + beta * std::exp(uhat[i]);
    }

    return H;
}

Eigen::VectorXd joint_envelope_gradient(const Eigen::VectorXd& theta,
                                        const Eigen::VectorXd& uhat) {
    const int m = static_cast<int>(uhat.size());
    const double mu = theta[0];
    const double inv_sigma2 = std::exp(-2.0 * theta[1]);
    const double lambda0 = std::exp(theta[2]);
    const double beta = std::exp(theta[3]);
    const Eigen::VectorXd y = make_y(m);

    Eigen::VectorXd g(4);
    g.setZero();

    for (int i = 0; i < m; ++i) {
        const double resid = y[i] - mu - uhat[i];
        const double expu = std::exp(uhat[i]);

        g[0] += -resid * inv_sigma2;
        g[1] += 1.0 - resid * resid * inv_sigma2;
        g[2] += 0.5 * lambda0 * uhat[i] * uhat[i];
        g[3] += beta * expu;
    }

    return g;
}

Eigen::MatrixXd du_dtheta_implicit(const Eigen::VectorXd& theta,
                                   const Eigen::VectorXd& uhat) {
    const int m = static_cast<int>(uhat.size());
    const double mu = theta[0];
    const double inv_sigma2 = std::exp(-2.0 * theta[1]);
    const double lambda0 = std::exp(theta[2]);
    const double beta = std::exp(theta[3]);
    const Eigen::VectorXd y = make_y(m);

    Eigen::MatrixXd du(m, 4);
    du.setZero();

    for (int i = 0; i < m; ++i) {
        const double resid = y[i] - mu - uhat[i];
        const double expu = std::exp(uhat[i]);
        const double H = inv_sigma2 + lambda0 + beta * expu;

        // f_u theta partials at fixed u.
        const double f_u_mu = inv_sigma2;
        const double f_u_log_sigma = 2.0 * resid * inv_sigma2;
        const double f_u_log_lambda0 = lambda0 * uhat[i];
        const double f_u_log_beta = beta * expu;

        du(i, 0) = -f_u_mu / H;
        du(i, 1) = -f_u_log_sigma / H;
        du(i, 2) = -f_u_log_lambda0 / H;
        du(i, 3) = -f_u_log_beta / H;
    }

    return du;
}

double joint_objective_double(const Eigen::VectorXd& theta,
                              const Eigen::VectorXd& uhat) {
    const int m = static_cast<int>(uhat.size());
    NonlinearDiagObjective objective{m};

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
    const Eigen::VectorXd uhat = solve_uhat(theta, m);
    const double joint = joint_objective_double(theta, uhat);
    const Eigen::MatrixXd H = Huu(theta, uhat);

    Eigen::LDLT<Eigen::MatrixXd> ldlt(H);
    if (ldlt.info() != Eigen::Success) {
        throw std::runtime_error("LDLT failed in laplace objective.");
    }

    const double logdet = ldlt.vectorD().array().log().sum();
    return joint + 0.5 * logdet;
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

        g[j] =
            (laplace_objective(plus, m) -
             laplace_objective(minus, m)) /
            (2.0 * h);
    }

    return g;
}

// Analytic total Hdot for extra validation.
Eigen::MatrixXd analytic_total_hdot(const Eigen::VectorXd& theta,
                                    const Eigen::VectorXd& uhat,
                                    const Eigen::MatrixXd& du,
                                    int theta_index) {
    const int m = static_cast<int>(uhat.size());
    const double inv_sigma2 = std::exp(-2.0 * theta[1]);
    const double lambda0 = std::exp(theta[2]);
    const double beta = std::exp(theta[3]);

    Eigen::MatrixXd Hdot = Eigen::MatrixXd::Zero(m, m);

    for (int i = 0; i < m; ++i) {
        const double expu = std::exp(uhat[i]);

        double value = beta * expu * du(i, theta_index);

        if (theta_index == 1) {
            value += -2.0 * inv_sigma2;
        } else if (theta_index == 2) {
            value += lambda0;
        } else if (theta_index == 3) {
            value += beta * expu;
        }

        Hdot(i, i) = value;
    }

    return Hdot;
}

void run_case(int m, const Eigen::VectorXd& theta) {
    const Eigen::VectorXd uhat = solve_uhat(theta, m);
    const Eigen::VectorXd gj = joint_envelope_gradient(theta, uhat);
    const Eigen::MatrixXd du = du_dtheta_implicit(theta, uhat);

    NonlinearDiagObjective objective{m};

    auto provider =
        quadra::laplace::make_had_quadra_replay_reuse_implicit_hdot_provider(
            objective,
            4,
            m,
            diagonal_pattern(m),
            std::vector<int>{0, 1, 2, 3},
            0.0);

    const auto Hdots = provider.compute_all_sparse(theta, uhat, du);

    for (int j = 0; j < 4; ++j) {
        const Eigen::MatrixXd got = Eigen::MatrixXd(Hdots[static_cast<size_t>(j)]);
        const Eigen::MatrixXd expected = analytic_total_hdot(theta, uhat, du, j);
        const double max_abs = (got - expected).cwiseAbs().maxCoeff();

        std::cout << "m=" << m
                  << " theta_index=" << j
                  << " max |AD Hdot_total - analytic| = "
                  << std::scientific << std::setprecision(10)
                  << max_abs << "\n";

        if (max_abs > 1.0e-10) {
            throw std::runtime_error(
                "implicit total Hdot does not match analytic total Hdot.");
        }
    }

    auto Hfn = [](const Eigen::VectorXd& th, const Eigen::VectorXd& uh) {
        return Huu(th, uh);
    };

    const Eigen::VectorXd g_exact =
        quadra::laplace::full_exact_laplace_gradient_implicit_cached_trace(
            gj,
            Hfn,
            provider,
            theta,
            uhat,
            du);

    const Eigen::VectorXd g_fd = central_fd_gradient(theta, m);

    const Eigen::VectorXd diff = (g_exact - g_fd).cwiseAbs();
    const double max_abs = diff.maxCoeff();

    std::cout << "\n=== implicit du Laplace gradient FD validation: m = "
              << m << " ===\n";
    std::cout << std::scientific << std::setprecision(10);
    std::cout << "exact = " << g_exact.transpose() << "\n";
    std::cout << "fd    = " << g_fd.transpose() << "\n";
    std::cout << "diff  = " << diff.transpose() << "\n";
    std::cout << "max abs diff = " << max_abs << "\n\n";

    if (max_abs > 2.0e-7) {
        throw std::runtime_error(
            "implicit exact Laplace gradient failed FD validation.");
    }
}

}  // namespace

int main() {
    Eigen::VectorXd theta(4);

    theta << 0.65, std::log(0.55), std::log(0.85), std::log(0.30);
    run_case(8, theta);
    run_case(25, theta);

    theta << 0.72, std::log(0.75), std::log(0.60), std::log(0.45);
    run_case(15, theta);

    std::cout << "implicit du exact Laplace gradient FD validation passed\n";
    return 0;
}
