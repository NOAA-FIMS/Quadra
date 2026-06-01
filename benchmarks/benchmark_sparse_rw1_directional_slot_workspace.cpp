#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../core/laplace/exact_gradient_workspace.hpp"
#include "../core/laplace/sparse_huu_factorization.hpp"

DECLARE_ADGRAPH()

using Clock = std::chrono::steady_clock;
using had::AReal;
using had::Real;

static double ms_between(const Clock::time_point& a, const Clock::time_point& b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

struct SparseRw1Objective {
  int m;
  template <class T>
  T operator()(const std::vector<T>& x) const {
    const T mu = x[0], log_sigma = x[1], log_lambda0 = x[2],
            log_lambda_rw = x[3], log_beta = x[4];
    const T inv_sigma2 = exp(T(-2.0) * log_sigma);
    const T lambda0 = exp(log_lambda0);
    const T lambda_rw = exp(log_lambda_rw);
    const T beta = exp(log_beta);
    T nll = T(0.0);
    for (int i = 0; i < m; ++i) {
      const double xd = static_cast<double>(i + 1);
      const T y = T(0.6 + 0.10 * std::sin(0.21 * xd) + 0.07 * std::cos(0.47 * xd));
      const T u = x[5 + i];
      const T resid = y - mu - u;
      nll = nll + T(0.5) * resid * resid * inv_sigma2 + log_sigma
          + T(0.5 * std::log(2.0 * 3.14159265358979323846))
          + T(0.5) * lambda0 * u * u + beta * exp(u);
    }
    for (int i = 1; i < m; ++i) {
      const T diff = x[5 + i] - x[5 + i - 1];
      nll = nll + T(0.5) * lambda_rw * diff * diff;
    }
    return nll;
  }
};

static Eigen::VectorXd make_y(int m) {
  Eigen::VectorXd y(m);
  for (int i = 0; i < m; ++i) {
    double x = static_cast<double>(i + 1);
    y[i] = 0.6 + 0.10 * std::sin(0.21 * x) + 0.07 * std::cos(0.47 * x);
  }
  return y;
}

static Eigen::SparseMatrix<double> Huu_sparse_direct(const Eigen::VectorXd& theta,
                                                     const Eigen::VectorXd& u) {
  int m = static_cast<int>(u.size());
  double inv_sigma2 = std::exp(-2.0 * theta[1]);
  double lambda0 = std::exp(theta[2]);
  double lambda_rw = std::exp(theta[3]);
  double beta = std::exp(theta[4]);
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

static Eigen::VectorXd random_gradient(const Eigen::VectorXd& theta,
                                       const Eigen::VectorXd& u) {
  int m = static_cast<int>(u.size());
  double mu = theta[0], inv_sigma2 = std::exp(-2.0 * theta[1]);
  double lambda0 = std::exp(theta[2]), lambda_rw = std::exp(theta[3]);
  double beta = std::exp(theta[4]);
  Eigen::VectorXd y = make_y(m);
  Eigen::VectorXd g = Eigen::VectorXd::Zero(m);
  for (int i = 0; i < m; ++i) {
    double resid = y[i] - mu - u[i];
    g[i] += -resid * inv_sigma2 + lambda0 * u[i] + beta * std::exp(u[i]);
  }
  for (int i = 1; i < m; ++i) {
    double diff = u[i] - u[i - 1];
    g[i] += lambda_rw * diff;
    g[i - 1] -= lambda_rw * diff;
  }
  return g;
}

static Eigen::VectorXd solve_uhat(const Eigen::VectorXd& theta, int m) {
  Eigen::VectorXd u = Eigen::VectorXd::Zero(m);
  for (int iter = 0; iter < 80; ++iter) {
    Eigen::VectorXd g = random_gradient(theta, u);
    Eigen::LDLT<Eigen::MatrixXd> ldlt(Eigen::MatrixXd(Huu_sparse_direct(theta, u)));
    Eigen::VectorXd step = ldlt.solve(g);
    u -= step;
    if (step.lpNorm<Eigen::Infinity>() < 1.0e-12) break;
  }
  return u;
}

static Eigen::VectorXd f_u_theta_column(const Eigen::VectorXd& theta,
                                        const Eigen::VectorXd& uhat,
                                        int j) {
  int m = static_cast<int>(uhat.size());
  double inv_sigma2 = std::exp(-2.0 * theta[1]);
  double lambda0 = std::exp(theta[2]), lambda_rw = std::exp(theta[3]);
  double beta = std::exp(theta[4]);
  Eigen::VectorXd y = make_y(m);
  Eigen::VectorXd col = Eigen::VectorXd::Zero(m);
  if (j == 0) { col.array() = inv_sigma2; return col; }
  if (j == 1) {
    for (int i = 0; i < m; ++i) col[i] = 2.0 * (y[i] - theta[0] - uhat[i]) * inv_sigma2;
    return col;
  }
  if (j == 2) { for (int i = 0; i < m; ++i) col[i] = lambda0 * uhat[i]; return col; }
  if (j == 3) {
    for (int i = 1; i < m; ++i) {
      double diff = uhat[i] - uhat[i - 1];
      col[i] += lambda_rw * diff; col[i - 1] -= lambda_rw * diff;
    }
    return col;
  }
  if (j == 4) { for (int i = 0; i < m; ++i) col[i] = beta * std::exp(uhat[i]); }
  return col;
}

struct SelectedTriInv {
  Eigen::VectorXd diag, sub;
  double operator()(int row, int col) const {
    if (row == col) return diag[row];
    int hi = std::max(row, col), lo = std::min(row, col);
    if (hi == lo + 1) return sub[hi - 1];
    return 0.0;
  }
};

static SelectedTriInv selected_inverse(Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>& ldlt, int m) {
  SelectedTriInv out;
  out.diag = Eigen::VectorXd::Zero(m);
  out.sub = Eigen::VectorXd::Zero(std::max(0, m - 1));
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(m);
  for (int j = 0; j < m; ++j) {
    rhs.setZero(); rhs[j] = 1.0;
    Eigen::VectorXd col = ldlt.solve(rhs);
    out.diag[j] = col[j];
    if (j > 0) out.sub[j - 1] = col[j - 1];
  }
  return out;
}

struct Row { double reverse=0, btree_trace=0, slot_read=0, slot_trace=0, diff=0; double queries=0, pushdots=0, inserts=0; };

static Row run_case(int m, int K, int reps) {
  Eigen::VectorXd theta(5);
  theta << 0.55, std::log(0.65), std::log(0.55), std::log(0.90), std::log(0.25);
  Eigen::VectorXd uhat = solve_uhat(theta, m);
  Eigen::SparseMatrix<double> H = Huu_sparse_direct(theta, uhat);
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(H);
  SelectedTriInv hinv = selected_inverse(ldlt, m);
  quadra::laplace::SparseHuuFactorization factor(H);
  Row avg;

  for (int r = 0; r < reps; ++r) {
    quadra::laplace::ExactGradientWorkspace ws;
    std::vector<AReal> th(5), u(static_cast<size_t>(m));
    ws.Build([&]() {
      for (int j = 0; j < 5; ++j) th[j] = AReal(theta[j]);
      for (int i = 0; i < m; ++i) u[static_cast<size_t>(i)] = AReal(uhat[i]);
      std::vector<AReal> x; x.reserve(static_cast<size_t>(5 + m));
      for (auto& v : th) x.push_back(v);
      for (auto& v : u) x.push_back(v);
      return SparseRw1Objective{m}(x);
    }, &th, &u);

    ws.PropagateBaseAdjoint();
    ws.SeedTotalDirections(static_cast<size_t>(K),
      [&](size_t k, Eigen::VectorXd& td, Eigen::VectorXd& ud) {
        int j = static_cast<int>(k % 5);
        td = Eigen::VectorXd::Zero(5); td[j] = 1.0;
        ud = -factor.solve(f_u_theta_column(theta, uhat, j));
      });

    had::ResetBatchDirectionalCounters();
    auto t0 = Clock::now();
    ws.PropagateDirectionalBatch();
    auto t1 = Clock::now();

    avg.queries += static_cast<double>(had::g_batch_query_count);
    avg.pushdots += static_cast<double>(had::g_batch_pushdot_count);
    avg.inserts += static_cast<double>(had::g_batch_insert_count);

    auto b0 = Clock::now();
    Eigen::VectorXd tb = Eigen::VectorXd::Zero(K);
    for (int k = 0; k < K; ++k) {
      for (int i = 0; i < m; ++i) {
        tb[k] += hinv.diag[i] * had::GetAdjointDotBatch(u[i], u[i], k);
        if (i > 0) tb[k] += 2.0 * hinv.sub[i-1] * had::GetAdjointDotBatch(u[i], u[i-1], k);
      }
    }
    auto b1 = Clock::now();

    auto s0 = Clock::now();
    std::vector<std::vector<double>> diag(static_cast<size_t>(K), std::vector<double>(static_cast<size_t>(m)));
    std::vector<std::vector<double>> sub(static_cast<size_t>(K), std::vector<double>(static_cast<size_t>(std::max(0,m-1))));
    for (int k = 0; k < K; ++k) {
      for (int i = 0; i < m; ++i) {
        diag[k][i] = had::GetAdjointDotBatch(u[i], u[i], k);
        if (i > 0) sub[k][i-1] = had::GetAdjointDotBatch(u[i], u[i-1], k);
      }
    }
    auto s1 = Clock::now();

    auto s2 = Clock::now();
    Eigen::VectorXd ts = Eigen::VectorXd::Zero(K);
    for (int k = 0; k < K; ++k) {
      for (int i = 0; i < m; ++i) {
        ts[k] += hinv.diag[i] * diag[k][i];
        if (i > 0) ts[k] += 2.0 * hinv.sub[i-1] * sub[k][i-1];
      }
    }
    auto s3 = Clock::now();

    avg.reverse += ms_between(t0,t1);
    avg.btree_trace += ms_between(b0,b1);
    avg.slot_read += ms_between(s0,s1);
    avg.slot_trace += ms_between(s2,s3);
    avg.diff = std::max(avg.diff, (tb-ts).cwiseAbs().maxCoeff());
  }

  avg.reverse /= reps; avg.btree_trace /= reps; avg.slot_read /= reps; avg.slot_trace /= reps;
  avg.queries /= reps; avg.pushdots /= reps; avg.inserts /= reps;
  return avg;
}

int main(int argc, char** argv) {
  int reps = argc > 1 ? std::stoi(argv[1]) : 10;
  int K = 5;
  std::cout << "Sparse RW1 directional slot workspace benchmark\\n";
  std::cout << "reps per case = " << reps << "\\nK = " << K << "\\n\\n";
  std::cout << std::setw(8) << "m" << std::setw(14) << "reverse"
            << std::setw(14) << "btree trace" << std::setw(14) << "slot read"
            << std::setw(14) << "slot trace" << std::setw(14) << "queries"
            << std::setw(14) << "pushdots" << std::setw(14) << "inserts"
            << std::setw(14) << "max diff" << "\\n";
  std::cout << std::scientific << std::setprecision(6);
  for (int m : {100,250,500}) {
    Row r = run_case(m, K, reps);
    std::cout << std::setw(8) << m << std::setw(14) << r.reverse
              << std::setw(14) << r.btree_trace << std::setw(14) << r.slot_read
              << std::setw(14) << r.slot_trace << std::setw(14) << r.queries
              << std::setw(14) << r.pushdots << std::setw(14) << r.inserts
              << std::setw(14) << r.diff << "\\n";
  }
  std::cout << "\\nBenchmark complete.\\n";
}
