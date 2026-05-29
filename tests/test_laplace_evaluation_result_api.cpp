#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../core/laplace/laplace_evaluation_result.hpp"

DECLARE_ADGRAPH()

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

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
                -resid * inv_sigma2 + lambda0 * u[i] + beta * expu;
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

Eigen::VectorXd f_u_theta_column(const Eigen::VectorXd& theta,
                                 const Eigen::VectorXd& uhat,
                                 int theta_index) {
    const int m = static_cast<int>(uhat.size());
    const double inv_sigma2 = std::exp(-2.0 * theta[1]);
    const double lambda0 = std::exp(theta[2]);
    const double beta = std::exp(theta[3]);
    const Eigen::VectorXd y = make_y(m);

    Eigen::VectorXd col(m);
    col.setZero();

    for (int i = 0; i < m; ++i) {
        const double resid = y[i] - theta[0] - uhat[i];
        const double expu = std::exp(uhat[i]);

        if (theta_index == 0) {
            col[i] = inv_sigma2;
        } else if (theta_index == 1) {
            col[i] = 2.0 * resid * inv_sigma2;
        } else if (theta_index == 2) {
            col[i] = lambda0 * uhat[i];
        } else if (theta_index == 3) {
            col[i] = beta * expu;
        }
    }

    return col;
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

        g[j] = (laplace_objective(plus, m) -
                laplace_objective(minus, m)) / (2.0 * h);
    }

    return g;
}

void run_case(int m, const Eigen::VectorXd& theta) {
    const Eigen::VectorXd uhat = solve_uhat(theta, m);
    const Eigen::MatrixXd H = Huu(theta, uhat);
    const double joint = joint_objective_double(theta, uhat);
    const Eigen::VectorXd gj = joint_envelope_gradient(theta, uhat);

    NonlinearDiagObjective objective{m};

    auto cross = [theta, uhat](int theta_index) -> Eigen::VectorXd {
        return f_u_theta_column(theta, uhat, theta_index);
    };

    const auto result =
        quadra::laplace::evaluate_laplace_objective_and_gradient(
            objective,
            cross,
            4,
            m,
            diagonal_pattern(m),
            std::vector<int>{0, 1, 2, 3},
            theta,
            uhat,
            H,
            joint,
            gj);

    const double fd_obj = laplace_objective(theta, m);
    const Eigen::VectorXd fd_grad = central_fd_gradient(theta, m);

    const double obj_diff = std::abs(result.objective - fd_obj);
    const double grad_diff = (result.gradient - fd_grad).cwiseAbs().maxCoeff();

    std::cout << "\n=== LaplaceEvaluationResult validation: m = "
              << m << " ===\n";
    std::cout << std::scientific << std::setprecision(10);
    std::cout << "objective      = " << result.objective << "\n";
    std::cout << "fd objective   = " << fd_obj << "\n";
    std::cout << "obj diff       = " << obj_diff << "\n";
    std::cout << "gradient       = " << result.gradient.transpose() << "\n";
    std::cout << "fd gradient    = " << fd_grad.transpose() << "\n";
    std::cout << "grad max diff  = " << grad_diff << "\n";

    if (obj_diff > 1.0e-12) {
        throw std::runtime_error("Laplace objective mismatch.");
    }

    if (grad_diff > 2.0e-7) {
        throw std::runtime_error("Laplace gradient mismatch.");
    }

    if (result.active_direction_count != 4) {
        throw std::runtime_error("active direction count mismatch.");
    }

    if (result.theta_dim != 4 || result.random_dim != m) {
        throw std::runtime_error("dimension metadata mismatch.");
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

    std::cout << "\nLaplaceEvaluationResult API tests passed\n";
    return 0;
}
