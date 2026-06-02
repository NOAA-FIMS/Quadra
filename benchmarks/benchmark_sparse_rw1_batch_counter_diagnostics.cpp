#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../core/had_quadra.hpp"
#include "../core/laplace/sparse_huu_factorization.hpp"

DECLARE_ADGRAPH()

namespace {

using Clock = std::chrono::steady_clock;
using had::AReal;
using had::Real;

constexpr double kPi = 3.141592653589793238462643383279502884;

struct SparseRw1Objective {
  int m;

  template <class T> T operator()(const std::vector<T> &x) const {
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
      const T y =
          T(0.6 + 0.10 * std::sin(0.21 * xd) + 0.07 * std::cos(0.47 * xd));
      const T u = x[5 + i];
      const T resid = y - mu - u;

      nll = nll + T(0.5) * resid * resid * inv_sigma2 + log_sigma +
            T(0.5 * std::log(2.0 * kPi)) + T(0.5) * lambda0 * u * u +
            beta * exp(u);
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
    y[i] = 0.6 + 0.10 * std::sin(0.21 * x) + 0.07 * std::cos(0.47 * x);
  }
  return y;
}

Eigen::SparseMatrix<double> Huu_sparse_direct(const Eigen::VectorXd &theta,
                                              const Eigen::VectorXd &u) {
  const int m = static_cast<int>(u.size());
  const double inv_sigma2 = std::exp(-2.0 * theta[1]);
  const double lambda0 = std::exp(theta[2]);
  const double lambda_rw = std::exp(theta[3]);
  const double beta = std::exp(theta[4]);

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(static_cast<size_t>(3 * m - 2));

  for (int i = 0; i < m; ++i) {
    double diag = inv_sigma2 + lambda0 + beta * std::exp(u[i]);
    if (i > 0) {
      diag += lambda_rw;
    }
    if (i + 1 < m) {
      diag += lambda_rw;
    }
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

Eigen::VectorXd random_gradient(const Eigen::VectorXd &theta,
                                const Eigen::VectorXd &u) {
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
    g[i] += -resid * inv_sigma2 + lambda0 * u[i] + beta * std::exp(u[i]);
  }

  for (int i = 1; i < m; ++i) {
    const double diff = u[i] - u[i - 1];
    g[i] += lambda_rw * diff;
    g[i - 1] -= lambda_rw * diff;
  }

  return g;
}

Eigen::VectorXd solve_uhat(const Eigen::VectorXd &theta, int m) {
  Eigen::VectorXd u = Eigen::VectorXd::Zero(m);

  for (int iter = 0; iter < 80; ++iter) {
    const Eigen::VectorXd g = random_gradient(theta, u);
    Eigen::LDLT<Eigen::MatrixXd> ldlt(
        Eigen::MatrixXd(Huu_sparse_direct(theta, u)));
    const Eigen::VectorXd step = ldlt.solve(g);
    u -= step;
    if (step.lpNorm<Eigen::Infinity>() < 1.0e-12) {
      break;
    }
  }

  return u;
}

Eigen::VectorXd f_u_theta_column(const Eigen::VectorXd &theta,
                                 const Eigen::VectorXd &uhat, int theta_index) {
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

template <class Fn> double elapsed_ms(Fn &&fn) {
  const auto start = Clock::now();
  volatile double sink = fn();
  (void)sink;
  const auto end = Clock::now();
  return std::chrono::duration<double, std::milli>(end - start).count();
}

struct Row {
  double total_ms = 0.0;
  std::uint64_t queries = 0;
  std::uint64_t pushdots = 0;
  std::uint64_t inserts = 0;
  int vertices = 0;
  int so_edges_nnz = 0;
};

Row run_case(int m, int K) {
  Eigen::VectorXd theta(5);
  theta << 0.55, std::log(0.65), std::log(0.55), std::log(0.90), std::log(0.25);

  const Eigen::VectorXd uhat = solve_uhat(theta, m);
  const Eigen::SparseMatrix<double> H = Huu_sparse_direct(theta, uhat);
  quadra::laplace::SparseHuuFactorization factor(H);

  auto direction_provider = [&](int theta_index) -> Eigen::VectorXd {
    return -factor.solve(f_u_theta_column(theta, uhat, theta_index));
  };

  had::ADGraph graph;
  had::g_ADGraph = &graph;

  std::vector<AReal> x(static_cast<size_t>(5 + m));

  for (int j = 0; j < 5; ++j) {
    x[static_cast<size_t>(j)] = AReal(theta[j]);
  }
  for (int i = 0; i < m; ++i) {
    x[static_cast<size_t>(5 + i)] = AReal(uhat[i]);
  }

  SparseRw1Objective objective{m};
  AReal f = objective(x);

  had::g_ADGraph->vertices[f.varId].w = 1.0;
  had::PropagateAdjoint();

  had::ResizeDirectionalBatch(K);

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

  Row out;
  out.vertices = static_cast<int>(had::g_ADGraph->vertices.size());

  out.total_ms = elapsed_ms([&]() {
    had::PropagateAdjointDirectionalBatch();
    return had::GetAdjointDotBatch(x.back(), x.front(), 0);
  });

  out.queries = had::g_batch_query_count;
  out.pushdots = had::g_batch_pushdot_count;
  out.inserts = had::g_batch_insert_count;

  int nnz = 0;
  for (const auto &tree : had::g_ADGraph->soEdges) {
    nnz += static_cast<int>(tree.nodes.size());
  }
  out.so_edges_nnz = nnz;

  return out;
}

} // namespace

int main() {
  const std::vector<int> m_values = {10, 25, 50, 100, 250};
  const std::vector<int> K_values = {1, 2, 4, 5};

  std::cout << "Sparse RW1 HAD batch counter diagnostics\n\n";

  std::cout << std::setw(8) << "m" << std::setw(8) << "K" << std::setw(14)
            << "total ms" << std::setw(12) << "vertices" << std::setw(12)
            << "so nnz" << std::setw(14) << "queries" << std::setw(14)
            << "pushdots" << std::setw(14) << "inserts" << std::setw(14)
            << "q/K" << std::setw(14) << "pd/K" << "\n";

  std::cout << std::scientific << std::setprecision(6);

  for (int m : m_values) {
    for (int K : K_values) {
      const Row r = run_case(m, K);

      std::cout << std::setw(8) << m << std::setw(8) << K << std::setw(14)
                << r.total_ms << std::setw(12) << r.vertices << std::setw(12)
                << r.so_edges_nnz << std::setw(14) << r.queries << std::setw(14)
                << r.pushdots << std::setw(14) << r.inserts << std::setw(14)
                << static_cast<double>(r.queries) / K << std::setw(14)
                << static_cast<double>(r.pushdots) / K << "\n";
    }
  }

  std::cout << "\nBenchmark complete.\n";
  return 0;
}
