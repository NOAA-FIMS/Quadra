#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../core/had_quadra.hpp"
#include "../core/laplace/exact_gradient_workspace.hpp"
#include "../core/had_graph_workspace.hpp"
#include "../core/laplace/sparse_huu_factorization.hpp"

DECLARE_ADGRAPH()

namespace {

using Clock = std::chrono::steady_clock;
using had::AReal;
using had::Real;

constexpr double kPi = 3.141592653589793238462643383279502884;

double ms_between(const Clock::time_point& a, const Clock::time_point& b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
}

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
        if (step.lpNorm<Eigen::Infinity>() < 1.0e-12) break;
    }

    return u;
}

double joint_objective_double(const Eigen::VectorXd& theta,
                              const Eigen::VectorXd& uhat) {
    const int m = static_cast<int>(uhat.size());
    SparseRw1Objective objective{m};

    std::vector<double> x(static_cast<size_t>(5 + m));
    for (int j = 0; j < 5; ++j) x[static_cast<size_t>(j)] = theta[j];
    for (int i = 0; i < m; ++i) x[static_cast<size_t>(5 + i)] = uhat[i];

    return objective(x);
}

Eigen::VectorXd joint_envelope_gradient(const Eigen::VectorXd& theta,
                                        const Eigen::VectorXd& uhat) {
    const int m = static_cast<int>(uhat.size());
    const double mu = theta[0];
    const double inv_sigma2 = std::exp(-2.0 * theta[1]);
    const double lambda0 = std::exp(theta[2]);
    const double lambda_rw = std::exp(theta[3]);
    const double beta = std::exp(theta[4]);
    const Eigen::VectorXd y = make_y(m);

    Eigen::VectorXd g = Eigen::VectorXd::Zero(5);

    for (int i = 0; i < m; ++i) {
        const double resid = y[i] - mu - uhat[i];
        const double expu = std::exp(uhat[i]);

        g[0] += -resid * inv_sigma2;
        g[1] += 1.0 - resid * resid * inv_sigma2;
        g[2] += 0.5 * lambda0 * uhat[i] * uhat[i];
        g[4] += beta * expu;
    }

    for (int i = 1; i < m; ++i) {
        const double diff = uhat[i] - uhat[i - 1];
        g[3] += 0.5 * lambda_rw * diff * diff;
    }

    return g;
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
        for (int i = 0; i < m; ++i) col[i] = lambda0 * uhat[i];
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
        for (int i = 0; i < m; ++i) col[i] = beta * std::exp(uhat[i]);
        return col;
    }

    return col;
}

struct ExactGradientResult {
    double objective = 0.0;
    Eigen::VectorXd gradient;
};

ExactGradientResult compute_exact_gradient_from_current_graph(
    const std::vector<AReal>& x,
    int m,
    int K,
    const Eigen::MatrixXd& Hinv,
    const Eigen::VectorXd& joint_grad,
    double joint_objective,
    double logdet) {

    Eigen::VectorXd grad = joint_grad;

    for (int k = 0; k < K; ++k) {
        double trace = 0.0;

        for (int i = 0; i < m; ++i) {
            const double hdot_ii =
                had::GetAdjointDotBatch(
                    x[static_cast<size_t>(5 + i)],
                    x[static_cast<size_t>(5 + i)],
                    k);

            trace += Hinv(i, i) * hdot_ii;

            if (i > 0) {
                const double hdot_i_im1 =
                    had::GetAdjointDotBatch(
                        x[static_cast<size_t>(5 + i)],
                        x[static_cast<size_t>(5 + i - 1)],
                        k);

                trace += 2.0 * Hinv(i, i - 1) * hdot_i_im1;
            }
        }

        grad[k] += 0.5 * trace;
    }

    ExactGradientResult out;
    out.objective = joint_objective + 0.5 * logdet;
    out.gradient = grad;
    return out;
}


struct ReusableRw1HadWorkspace {
    int m = 0;
    int K = 0;
    quadra::HadGraphWorkspace had_workspace;
    std::vector<AReal> x;
    AReal f;

    ReusableRw1HadWorkspace(int m_, int K_,
                            const Eigen::VectorXd& theta,
                            const Eigen::VectorXd& uhat)
        : m(m_), K(K_), x(static_cast<size_t>(5 + m_)) {
        build(theta, uhat);
    }

    void build(const Eigen::VectorXd& theta,
               const Eigen::VectorXd& uhat) {
        f = had_workspace.Build([&]() {
            for (int j = 0; j < 5; ++j) {
                x[static_cast<size_t>(j)] = AReal(theta[j]);
            }
            for (int i = 0; i < m; ++i) {
                x[static_cast<size_t>(5 + i)] = AReal(uhat[i]);
            }

            SparseRw1Objective objective{m};
            return objective(x);
        });

        had_workspace.ResizeDirectionalBatch(static_cast<std::size_t>(K));
    }

    void propagate_base_adjoint() {
        had_workspace.PropagateAdjoint(f.varId);
    }

    void seed_directions(const Eigen::VectorXd& theta,
                         const Eigen::VectorXd& uhat,
                         quadra::laplace::SparseHuuFactorization& factor) {
        had_workspace.Activate();
        had_workspace.ResizeDirectionalBatch(static_cast<std::size_t>(K));

        auto direction_provider = [&](int theta_index) -> Eigen::VectorXd {
            return -factor.solve(f_u_theta_column(theta, uhat, theta_index));
        };

        for (int k = 0; k < K; ++k) {
            const int theta_index = k % 5;
            const Eigen::VectorXd udir = direction_provider(theta_index);

            for (int j = 0; j < 5; ++j) {
                had::SetARealDotBatch(x[static_cast<size_t>(j)], k,
                                      j == theta_index ? 1.0 : 0.0);
            }
            for (int i = 0; i < m; ++i) {
                had::SetARealDotBatch(x[static_cast<size_t>(5 + i)], k, udir[i]);
            }
        }
    }

    void propagate_directional_batch() {
        had_workspace.PropagateAdjointDirectionalBatch();
    }

    int vertex_count() const {
        return static_cast<int>(had_workspace.Graph().vertices.size());
    }
};


ExactGradientResult exact_gradient_rebuild(
    int m,
    int K,
    const Eigen::VectorXd& theta,
    const Eigen::VectorXd& uhat,
    const Eigen::MatrixXd& Hinv,
    const Eigen::VectorXd& joint_grad,
    double joint_objective,
    double logdet,
    quadra::laplace::SparseHuuFactorization& factor) {

    had::ADGraph graph;
    had::g_ADGraph = &graph;

    std::vector<AReal> x(static_cast<size_t>(5 + m));
    for (int j = 0; j < 5; ++j) x[static_cast<size_t>(j)] = AReal(theta[j]);
    for (int i = 0; i < m; ++i) x[static_cast<size_t>(5 + i)] = AReal(uhat[i]);

    SparseRw1Objective objective{m};
    AReal f = objective(x);

    graph.vertices[f.varId].w = Real(1.0);
    had::PropagateAdjoint();

    had::ResizeDirectionalBatch(K);

    auto direction_provider = [&](int theta_index) -> Eigen::VectorXd {
        return -factor.solve(f_u_theta_column(theta, uhat, theta_index));
    };

    for (int k = 0; k < K; ++k) {
        const int theta_index = k % 5;
        const Eigen::VectorXd udir = direction_provider(theta_index);

        for (int j = 0; j < 5; ++j) {
            had::SetARealDotBatch(x[static_cast<size_t>(j)], k,
                                  j == theta_index ? 1.0 : 0.0);
        }
        for (int i = 0; i < m; ++i) {
            had::SetARealDotBatch(x[static_cast<size_t>(5 + i)], k, udir[i]);
        }
    }

    // PropagateAdjoint() consumes/clears reverse adjoints while building the
    // base Hessian. Match HadGraphWorkspace behavior by reseeding the output
    // adjoint before the directional reverse sweep.
    graph.vertices[f.varId].w = Real(1.0);
    had::PropagateAdjointDirectionalBatch();

    return compute_exact_gradient_from_current_graph(
        x, m, K, Hinv, joint_grad, joint_objective, logdet);
}


struct ExactGradientWorkspaceRw1Adapter {
    int m = 0;
    int K = 0;
    quadra::laplace::ExactGradientWorkspace workspace;
    std::vector<AReal> theta_vars;
    std::vector<AReal> random_vars;
    AReal f;

    ExactGradientWorkspaceRw1Adapter(int m_, int K_,
                                     const Eigen::VectorXd& theta,
                                     const Eigen::VectorXd& uhat)
        : m(m_),
          K(K_),
          theta_vars(5),
          random_vars(static_cast<std::size_t>(m_)) {
        build(theta, uhat);
    }

    void build(const Eigen::VectorXd& theta,
               const Eigen::VectorXd& uhat) {
        f = workspace.Build(
            [&]() {
                for (int j = 0; j < 5; ++j) {
                    theta_vars[static_cast<std::size_t>(j)] = AReal(theta[j]);
                }
                for (int i = 0; i < m; ++i) {
                    random_vars[static_cast<std::size_t>(i)] = AReal(uhat[i]);
                }

                std::vector<AReal> x;
                x.reserve(static_cast<std::size_t>(5 + m));

                for (auto& v : theta_vars) {
                    x.push_back(v);
                }
                for (auto& v : random_vars) {
                    x.push_back(v);
                }

                SparseRw1Objective objective{m};
                return objective(x);
            },
            &theta_vars,
            &random_vars);

        workspace.ResizeDirectionalBatch(static_cast<std::size_t>(K));
    }

    void propagate_base_adjoint() {
        workspace.PropagateBaseAdjoint();
    }

    void seed_directions(const Eigen::VectorXd& theta,
                         const Eigen::VectorXd& uhat,
                         quadra::laplace::SparseHuuFactorization& factor) {
        workspace.SeedTotalDirections(
            static_cast<std::size_t>(K),
            [&](std::size_t k,
                Eigen::VectorXd& theta_direction,
                Eigen::VectorXd& random_direction) {
                const int theta_index = static_cast<int>(k % 5);

                theta_direction = Eigen::VectorXd::Zero(5);
                theta_direction[theta_index] = 1.0;

                random_direction =
                    -factor.solve(
                        f_u_theta_column(theta, uhat, theta_index));
            });
    }

    void propagate_directional_batch() {
        workspace.PropagateDirectionalBatch();
    }

    ExactGradientResult exact_gradient(const Eigen::MatrixXd& Hinv,
                                       const Eigen::VectorXd& joint_grad,
                                       double joint_objective,
                                       double logdet) {
        const auto pattern =
            quadra::laplace::MakeTridiagonalHdotPattern(m);

        const auto assembled =
            workspace.AssembleExactGradient(
                joint_objective,
                logdet,
                joint_grad,
                [&](int row, int col) {
                    return Hinv(row, col);
                },
                pattern);

        ExactGradientResult out;
        out.objective = assembled.objective;
        out.gradient = assembled.gradient;
        return out;
    }

    int vertex_count() const {
        return static_cast<int>(
            workspace.HadWorkspace().Graph().vertices.size());
    }
};

struct Row {
    double rebuild_ms = 0.0;
    double reuse_ms = 0.0;
    double speedup = 0.0;
    double grad_diff = 0.0;
    double obj_diff = 0.0;
    int vertices = 0;
};

Row run_case(int m, int K, int reps) {
    Eigen::VectorXd theta(5);
    theta << 0.55, std::log(0.65), std::log(0.55), std::log(0.90), std::log(0.25);

    const Eigen::VectorXd uhat = solve_uhat(theta, m);
    const Eigen::SparseMatrix<double> H = Huu_sparse_direct(theta, uhat);

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(H);
    const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(m, m);
    const Eigen::MatrixXd Hinv = ldlt.solve(I);

    double logdet = 0.0;
    Eigen::MatrixXd Hdense(H);
    Eigen::LLT<Eigen::MatrixXd> llt(Hdense);
    const auto& L = llt.matrixL();
    for (int i = 0; i < m; ++i) logdet += 2.0 * std::log(L(i, i));

    const double joint_obj = joint_objective_double(theta, uhat);
    const Eigen::VectorXd joint_grad = joint_envelope_gradient(theta, uhat);

    quadra::laplace::SparseHuuFactorization factor(H);

    Row out;
    ExactGradientResult last_rebuild;
    ExactGradientResult last_reuse;

    const auto t0 = Clock::now();
    for (int r = 0; r < reps; ++r) {
        last_rebuild = exact_gradient_rebuild(
            m, K, theta, uhat, Hinv, joint_grad, joint_obj, logdet, factor);
    }
    const auto t1 = Clock::now();

    ExactGradientWorkspaceRw1Adapter workspace(m, K, theta, uhat);
    out.vertices = workspace.vertex_count();

    const auto t2 = Clock::now();
    for (int r = 0; r < reps; ++r) {
        workspace.propagate_base_adjoint();
        workspace.seed_directions(theta, uhat, factor);
        workspace.propagate_directional_batch();

        last_reuse = workspace.exact_gradient(
            Hinv, joint_grad, joint_obj, logdet);
    }
    const auto t3 = Clock::now();

    out.rebuild_ms = ms_between(t0, t1) / static_cast<double>(reps);
    out.reuse_ms = ms_between(t2, t3) / static_cast<double>(reps);
    out.speedup = out.rebuild_ms / out.reuse_ms;
    out.grad_diff = (last_rebuild.gradient - last_reuse.gradient).cwiseAbs().maxCoeff();
    out.obj_diff = std::abs(last_rebuild.objective - last_reuse.objective);

    return out;
}

}  // namespace

int main(int argc, char** argv) {
    int reps = 10;
    if (argc > 1) reps = std::stoi(argv[1]);

    const std::vector<int> m_values = {100, 250, 500};
    const int K = 5;

    std::cout << "Sparse RW1 exact-gradient graph reuse benchmark\n";
    std::cout << "reps per case = " << reps << "\n";
    std::cout << "K = " << K << "\n\n";

    std::cout << std::setw(8) << "m"
              << std::setw(12) << "vertices"
              << std::setw(16) << "rebuild ms"
              << std::setw(16) << "reuse ms"
              << std::setw(14) << "speedup"
              << std::setw(16) << "grad diff"
              << std::setw(16) << "obj diff"
              << "\n";

    std::cout << std::scientific << std::setprecision(6);

    for (int m : m_values) {
        const Row r = run_case(m, K, reps);

        std::cout << std::setw(8) << m
                  << std::setw(12) << r.vertices
                  << std::setw(16) << r.rebuild_ms
                  << std::setw(16) << r.reuse_ms
                  << std::setw(14) << r.speedup
                  << std::setw(16) << r.grad_diff
                  << std::setw(16) << r.obj_diff
                  << "\n";
    }

    std::cout << "\nBenchmark complete.\n";
    return 0;
}
