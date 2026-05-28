#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "../core/laplace/laplace_evaluator_exact_gradient_fd_adapter.hpp"

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

struct GradientComparison {
    Eigen::VectorXd exact;
    Eigen::VectorXd fd;
    Eigen::VectorXd abs_diff;
    Eigen::VectorXd rel_diff;
    double max_abs_diff = 0.0;
    double max_rel_diff = 0.0;
};

class GaussianRandomEffectsModel {
public:
    explicit GaussianRandomEffectsModel(int m) : y_(m) {
        if (m <= 0) {
            throw std::invalid_argument("m must be positive.");
        }

        // Deterministic, mildly structured data. No RNG required.
        for (int i = 0; i < m; ++i) {
            const double x = static_cast<double>(i + 1);
            y_[i] = 1.35 + 0.15 * std::sin(0.7 * x) + 0.03 * std::cos(1.3 * x);
        }
    }

    int random_dim() const {
        return static_cast<int>(y_.size());
    }

    // theta = [mu, log_sigma, log_tau]
    double joint_nll(const Eigen::VectorXd& theta,
                     const Eigen::VectorXd& u) const {
        check_theta(theta);
        if (u.size() != y_.size()) {
            throw std::invalid_argument("u length must equal y length.");
        }

        const double mu = theta[0];
        const double sigma = std::exp(theta[1]);
        const double tau = std::exp(theta[2]);
        const double sigma2 = sigma * sigma;
        const double tau2 = tau * tau;

        double nll = 0.0;

        for (int i = 0; i < y_.size(); ++i) {
            const double resid = y_[i] - mu - u[i];
            nll += 0.5 * resid * resid / sigma2 + std::log(sigma)
                   + 0.5 * std::log(2.0 * kPi);

            nll += 0.5 * u[i] * u[i] / tau2 + std::log(tau)
                   + 0.5 * std::log(2.0 * kPi);
        }

        return nll;
    }

    // Closed-form inner optimum for u_i:
    //
    // minimize
    //   0.5 (y_i - mu - u_i)^2 / sigma^2 + 0.5 u_i^2 / tau^2
    //
    // gives
    //   uhat_i = tau^2 / (sigma^2 + tau^2) * (y_i - mu)
    Eigen::VectorXd optimize_random_effects(const Eigen::VectorXd& theta) const {
        check_theta(theta);

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

    // H_uu(theta, uhat). For this Gaussian model it does not depend on uhat,
    // but the argument is retained to match the real evaluator contract.
    Eigen::MatrixXd hessian_uu_at_fixed_uhat(
        const Eigen::VectorXd& theta,
        const Eigen::VectorXd& uhat) const {
        check_theta(theta);
        if (uhat.size() != y_.size()) {
            throw std::invalid_argument("uhat length must equal y length.");
        }

        const double sigma2 = std::exp(2.0 * theta[1]);
        const double tau2 = std::exp(2.0 * theta[2]);

        const double h = 1.0 / sigma2 + 1.0 / tau2;
        return h * Eigen::MatrixXd::Identity(y_.size(), y_.size());
    }

    // Envelope joint gradient wrt theta at fixed optimized uhat.
    //
    // Do not differentiate through uhat here.
    Eigen::VectorXd joint_envelope_gradient(
        const Eigen::VectorXd& theta,
        const Eigen::VectorXd& uhat) const {
        check_theta(theta);
        if (uhat.size() != y_.size()) {
            throw std::invalid_argument("uhat length must equal y length.");
        }

        const double mu = theta[0];
        const double sigma2 = std::exp(2.0 * theta[1]);
        const double tau2 = std::exp(2.0 * theta[2]);

        Eigen::VectorXd g(3);
        g.setZero();

        for (int i = 0; i < y_.size(); ++i) {
            const double resid = y_[i] - mu - uhat[i];

            // d/dmu: 0.5 resid^2 / sigma^2
            // resid = y - mu - u
            g[0] += -resid / sigma2;

            // d/dlog_sigma:
            // 0.5 resid^2 exp(-2 log_sigma) + log_sigma
            // derivative = -resid^2 / sigma^2 + 1
            g[1] += 1.0 - resid * resid / sigma2;

            // d/dlog_tau:
            // 0.5 u^2 exp(-2 log_tau) + log_tau
            // derivative = -u^2 / tau^2 + 1
            g[2] += 1.0 - uhat[i] * uhat[i] / tau2;
        }

        return g;
    }

    double laplace_objective(const Eigen::VectorXd& theta) const {
        const Eigen::VectorXd uhat = optimize_random_effects(theta);
        const double joint = joint_nll(theta, uhat);
        const Eigen::MatrixXd H = hessian_uu_at_fixed_uhat(theta, uhat);

        Eigen::LDLT<Eigen::MatrixXd> ldlt(H);
        if (ldlt.info() != Eigen::Success) {
            throw std::runtime_error("LDLT factorization failed.");
        }

        // H is diagonal/SPD in this test model. Use determinant for simple,
        // transparent validation output.
        const double det = H.determinant();
        if (!(det > 0.0) || !std::isfinite(det)) {
            throw std::runtime_error("H determinant must be positive and finite.");
        }

        // Constants do not matter for gradients. We include the 0.5 log |H|
        // term only because that is the piece whose derivative is being tested.
        return joint + 0.5 * std::log(det);
    }

    Eigen::VectorXd exact_gradient_fd_trace(const Eigen::VectorXd& theta) const {
        const Eigen::VectorXd uhat = optimize_random_effects(theta);

        quadra::laplace::FullExactLaplaceGradientFDOptions options;
        options.step = 1.0e-6;
        options.relative_step = true;

        auto adapter =
            quadra::laplace::make_laplace_evaluator_exact_gradient_fd_adapter(
                [this](const Eigen::VectorXd& th,
                       const Eigen::VectorXd& uh) {
                    return this->joint_envelope_gradient(th, uh);
                },
                [this](const Eigen::VectorXd& th,
                       const Eigen::VectorXd& uh) {
                    return this->hessian_uu_at_fixed_uhat(th, uh);
                },
                options);

        return adapter(theta, uhat);
    }

private:
    Eigen::VectorXd y_;

    static void check_theta(const Eigen::VectorXd& theta) {
        if (theta.size() != 3) {
            throw std::invalid_argument("theta must have length 3.");
        }
        for (int j = 0; j < theta.size(); ++j) {
            if (!std::isfinite(theta[j])) {
                throw std::invalid_argument("theta contains non-finite values.");
            }
        }
    }
};

Eigen::VectorXd central_fd_gradient(
    const std::function<double(const Eigen::VectorXd&)>& f,
    const Eigen::VectorXd& theta,
    double step = 1.0e-5) {
    Eigen::VectorXd g(theta.size());
    g.setZero();

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

GradientComparison compare_gradient(const GaussianRandomEffectsModel& model,
                                    const Eigen::VectorXd& theta) {
    GradientComparison out;

    out.exact = model.exact_gradient_fd_trace(theta);
    out.fd = central_fd_gradient(
        [&](const Eigen::VectorXd& th) {
            return model.laplace_objective(th);
        },
        theta,
        1.0e-5);

    out.abs_diff = (out.exact - out.fd).cwiseAbs();
    out.rel_diff.resize(theta.size());

    for (int j = 0; j < theta.size(); ++j) {
        const double denom =
            std::max(1.0, std::max(std::abs(out.exact[j]), std::abs(out.fd[j])));
        out.rel_diff[j] = out.abs_diff[j] / denom;
    }

    out.max_abs_diff = out.abs_diff.maxCoeff();
    out.max_rel_diff = out.rel_diff.maxCoeff();

    return out;
}

void print_comparison(int m,
                      const Eigen::VectorXd& theta,
                      const GradientComparison& cmp) {
    std::cout << "\n=== Gaussian random-effects validation: m = " << m
              << " ===\n";
    std::cout << "theta = [" << theta.transpose() << "]\n\n";

    std::cout << std::setw(12) << "parameter"
              << std::setw(20) << "exact"
              << std::setw(20) << "fd objective"
              << std::setw(20) << "abs diff"
              << std::setw(20) << "rel diff"
              << "\n";

    const char* names[3] = {"mu", "log_sigma", "log_tau"};

    std::cout << std::scientific << std::setprecision(10);
    for (int j = 0; j < 3; ++j) {
        std::cout << std::setw(12) << names[j]
                  << std::setw(20) << cmp.exact[j]
                  << std::setw(20) << cmp.fd[j]
                  << std::setw(20) << cmp.abs_diff[j]
                  << std::setw(20) << cmp.rel_diff[j]
                  << "\n";
    }

    std::cout << "\nmax abs diff = " << cmp.max_abs_diff << "\n";
    std::cout << "max rel diff = " << cmp.max_rel_diff << "\n";
}

void run_case(int m, const Eigen::VectorXd& theta) {
    GaussianRandomEffectsModel model(m);
    const GradientComparison cmp = compare_gradient(model, theta);
    print_comparison(m, theta, cmp);

    // This is a full objective-vs-gradient FD comparison, so do not set this
    // too aggressively. 1e-5 is already quite good for nested finite
    // differences involving a log determinant.
    const double max_abs_tol = 1.0e-5;
    const double max_rel_tol = 1.0e-5;

    if (cmp.max_abs_diff > max_abs_tol && cmp.max_rel_diff > max_rel_tol) {
        std::cerr << "\nFAILED: exact Laplace gradient did not match FD "
                  << "of full Laplace objective.\n"
                  << "  max_abs_diff = " << cmp.max_abs_diff << "\n"
                  << "  max_rel_diff = " << cmp.max_rel_diff << "\n";
        throw std::runtime_error("gradient validation failed");
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

    std::cout << "\nexact Laplace gradient vs FD objective validation passed\n";
    return 0;
}
