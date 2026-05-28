#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#include "../core/laplace/full_exact_laplace_gradient_fd.hpp"
#include "../core/laplace/full_exact_laplace_gradient_hdot.hpp"

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

void expect_near(double got, double expected, double tol, const char* label) {
    const double err = std::abs(got - expected);
    if (!(err <= tol)) {
        std::cerr << "FAILED: " << label << "\n"
                  << "  got      = " << got << "\n"
                  << "  expected = " << expected << "\n"
                  << "  abs err  = " << err << "\n"
                  << "  tol      = " << tol << "\n";
        throw std::runtime_error(label);
    }
}

class GaussianRandomEffectsModel {
public:
    explicit GaussianRandomEffectsModel(int m) : y_(m) {
        for (int i = 0; i < m; ++i) {
            const double x = static_cast<double>(i + 1);
            y_[i] = 1.35 + 0.15 * std::sin(0.7 * x) + 0.03 * std::cos(1.3 * x);
        }
    }

    Eigen::VectorXd optimize_random_effects(const Eigen::VectorXd& theta) const {
        const double mu = theta[0];
        const double sigma2 = std::exp(2.0 * theta[1]);
        const double tau2 = std::exp(2.0 * theta[2]);
        const double shrink = tau2 / (sigma2 + tau2);

        Eigen::VectorXd uhat(y_.size());
        for (int i = 0; i < y_.size(); ++i) {
            uhat[i] = shrink * (y_[i] - mu);
        }
        return uhat;
    }

    double joint_nll(const Eigen::VectorXd& theta,
                     const Eigen::VectorXd& u) const {
        const double mu = theta[0];
        const double sigma2 = std::exp(2.0 * theta[1]);
        const double tau2 = std::exp(2.0 * theta[2]);

        double nll = 0.0;
        for (int i = 0; i < y_.size(); ++i) {
            const double resid = y_[i] - mu - u[i];

            nll += 0.5 * resid * resid / sigma2
                   + theta[1]
                   + 0.5 * std::log(2.0 * kPi);

            nll += 0.5 * u[i] * u[i] / tau2
                   + theta[2]
                   + 0.5 * std::log(2.0 * kPi);
        }
        return nll;
    }

    Eigen::VectorXd joint_envelope_gradient(
        const Eigen::VectorXd& theta,
        const Eigen::VectorXd& uhat) const {
        const double mu = theta[0];
        const double sigma2 = std::exp(2.0 * theta[1]);
        const double tau2 = std::exp(2.0 * theta[2]);

        Eigen::VectorXd g(3);
        g.setZero();

        for (int i = 0; i < y_.size(); ++i) {
            const double resid = y_[i] - mu - uhat[i];

            g[0] += -resid / sigma2;
            g[1] += 1.0 - resid * resid / sigma2;
            g[2] += 1.0 - uhat[i] * uhat[i] / tau2;
        }

        return g;
    }

    Eigen::MatrixXd hessian_uu_at_fixed_uhat(
        const Eigen::VectorXd& theta,
        const Eigen::VectorXd& uhat) const {
        (void)uhat;

        const double sigma2 = std::exp(2.0 * theta[1]);
        const double tau2 = std::exp(2.0 * theta[2]);
        const double h = 1.0 / sigma2 + 1.0 / tau2;

        return h * Eigen::MatrixXd::Identity(y_.size(), y_.size());
    }

    // Analytic Hdot_j = dH_uu/dtheta_j.
    //
    // H_uu = (exp(-2 log_sigma) + exp(-2 log_tau)) I
    //
    // dH/dmu        = 0
    // dH/dlog_sigma = -2 exp(-2 log_sigma) I
    // dH/dlog_tau   = -2 exp(-2 log_tau) I
    Eigen::MatrixXd analytic_hdot(
        const Eigen::VectorXd& theta,
        const Eigen::VectorXd& uhat,
        int theta_index) const {
        (void)uhat;

        Eigen::MatrixXd Hdot =
            Eigen::MatrixXd::Zero(y_.size(), y_.size());

        if (theta_index == 1) {
            Hdot.diagonal().array() =
                -2.0 * std::exp(-2.0 * theta[1]);
        } else if (theta_index == 2) {
            Hdot.diagonal().array() =
                -2.0 * std::exp(-2.0 * theta[2]);
        }

        return Hdot;
    }

    double laplace_objective(const Eigen::VectorXd& theta) const {
        const Eigen::VectorXd uhat = optimize_random_effects(theta);
        const double joint = joint_nll(theta, uhat);
        const Eigen::MatrixXd H = hessian_uu_at_fixed_uhat(theta, uhat);
        return joint + 0.5 * std::log(H.determinant());
    }

private:
    Eigen::VectorXd y_;
};

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
    GaussianRandomEffectsModel model(m);
    const Eigen::VectorXd uhat = model.optimize_random_effects(theta);
    const Eigen::VectorXd joint_grad =
        model.joint_envelope_gradient(theta, uhat);

    auto hessian = [&model](const Eigen::VectorXd& th,
                            const Eigen::VectorXd& uh) {
        return model.hessian_uu_at_fixed_uhat(th, uh);
    };

    auto hessian_theta_only = [&](const Eigen::VectorXd& th) {
        return model.hessian_uu_at_fixed_uhat(th, uhat);
    };

    auto hdot = [&model](const Eigen::VectorXd& th,
                         const Eigen::VectorXd& uh,
                         int j) {
        return model.analytic_hdot(th, uh, j);
    };

    const Eigen::VectorXd grad_hdot =
        quadra::laplace::full_exact_laplace_gradient_hdot(
            joint_grad,
            hessian,
            hdot,
            theta,
            uhat);

    quadra::laplace::FullExactLaplaceGradientFDOptions fd_options;
    fd_options.step = 1.0e-6;
    fd_options.relative_step = true;

    const Eigen::VectorXd grad_fd_trace =
        quadra::laplace::full_exact_laplace_gradient_fd(
            joint_grad,
            hessian_theta_only,
            theta,
            fd_options);

    const Eigen::VectorXd grad_fd_objective =
        central_fd_gradient(
            [&](const Eigen::VectorXd& th) {
                return model.laplace_objective(th);
            },
            theta);

    const Eigen::VectorXd diff_trace =
        (grad_hdot - grad_fd_trace).cwiseAbs();
    const Eigen::VectorXd diff_objective =
        (grad_hdot - grad_fd_objective).cwiseAbs();

    std::cout << "\n=== Hdot provider validation: m = " << m << " ===\n";
    std::cout << "theta = [" << theta.transpose() << "]\n\n";
    std::cout << std::setw(12) << "parameter"
              << std::setw(20) << "hdot exact"
              << std::setw(20) << "fd trace"
              << std::setw(20) << "fd objective"
              << std::setw(20) << "|hdot-fdobj|"
              << "\n";

    const char* names[3] = {"mu", "log_sigma", "log_tau"};

    std::cout << std::scientific << std::setprecision(10);
    for (int j = 0; j < 3; ++j) {
        std::cout << std::setw(12) << names[j]
                  << std::setw(20) << grad_hdot[j]
                  << std::setw(20) << grad_fd_trace[j]
                  << std::setw(20) << grad_fd_objective[j]
                  << std::setw(20) << diff_objective[j]
                  << "\n";
    }

    std::cout << "max |hdot - fd trace|     = "
              << diff_trace.maxCoeff() << "\n";
    std::cout << "max |hdot - fd objective| = "
              << diff_objective.maxCoeff() << "\n";

    if (diff_trace.maxCoeff() > 1.0e-7) {
        throw std::runtime_error("Hdot provider does not match FD trace path.");
    }

    if (diff_objective.maxCoeff() > 1.0e-5) {
        throw std::runtime_error("Hdot provider does not match FD objective.");
    }
}

}  // namespace

int main() {
    Eigen::VectorXd theta(3);

    theta << 1.1, std::log(0.45), std::log(0.8);
    run_case(10, theta);
    run_case(50, theta);
    run_case(100, theta);

    theta << 1.25, std::log(0.7), std::log(0.35);
    run_case(25, theta);

    std::cout << "\nexact Laplace gradient Hdot-provider validation passed\n";
    return 0;
}
