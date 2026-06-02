#include "../core/laplace/had_quadra_replay_reuse_lazy_implicit_hdot_provider.hpp"
#include "../core/laplace/sparse_huu_factorization.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
DECLARE_ADGRAPH()
using Clock = std::chrono::steady_clock;
constexpr double kPi = 3.141592653589793238462643383279502884;
struct SparseRw1Objective {
  int m;
  template <class T> T operator()(const std::vector<T> &x) const {
    T mu = x[0], ls = x[1], ll0 = x[2], llrw = x[3], lb = x[4];
    T is2 = exp(T(-2.0) * ls), l0 = exp(ll0), lrw = exp(llrw), b = exp(lb),
      nll = T(0);
    for (int i = 0; i < m; ++i) {
      double xd = i + 1;
      T y = T(0.6 + 0.10 * std::sin(0.21 * xd) + 0.07 * std::cos(0.47 * xd));
      T u = x[5 + i];
      T r = y - mu - u;
      nll += T(.5) * r * r * is2 + ls + T(.5 * std::log(2.0 * kPi)) +
             T(.5) * l0 * u * u + b * exp(u);
    }
    for (int i = 1; i < m; ++i) {
      T d = x[5 + i] - x[5 + i - 1];
      nll += T(.5) * lrw * d * d;
    }
    return nll;
  }
};
quadra::laplace::RandomHessianPattern pat(int m) {
  quadra::laplace::RandomHessianPattern p;
  for (int i = 0; i < m; ++i) {
    p.emplace_back(i, i);
    if (i > 0)
      p.emplace_back(i, i - 1);
  }
  return p;
}
Eigen::VectorXd yvec(int m) {
  Eigen::VectorXd y(m);
  for (int i = 0; i < m; ++i) {
    double x = i + 1;
    y[i] = 0.6 + 0.10 * std::sin(0.21 * x) + 0.07 * std::cos(0.47 * x);
  }
  return y;
}
Eigen::SparseMatrix<double> Huu(const Eigen::VectorXd &t,
                                const Eigen::VectorXd &u) {
  int m = u.size();
  double is2 = std::exp(-2 * t[1]), l0 = std::exp(t[2]), lrw = std::exp(t[3]),
         b = std::exp(t[4]);
  std::vector<Eigen::Triplet<double>> tr;
  tr.reserve(3 * m - 2);
  for (int i = 0; i < m; ++i) {
    double d = is2 + l0 + b * std::exp(u[i]);
    if (i > 0)
      d += lrw;
    if (i + 1 < m)
      d += lrw;
    tr.emplace_back(i, i, d);
  }
  for (int i = 1; i < m; ++i) {
    tr.emplace_back(i, i - 1, -lrw);
    tr.emplace_back(i - 1, i, -lrw);
  }
  Eigen::SparseMatrix<double> H(m, m);
  H.setFromTriplets(tr.begin(), tr.end());
  H.makeCompressed();
  return H;
}
Eigen::VectorXd gu(const Eigen::VectorXd &t, const Eigen::VectorXd &u) {
  int m = u.size();
  double mu = t[0], is2 = std::exp(-2 * t[1]), l0 = std::exp(t[2]),
         lrw = std::exp(t[3]), b = std::exp(t[4]);
  auto y = yvec(m);
  Eigen::VectorXd g = Eigen::VectorXd::Zero(m);
  for (int i = 0; i < m; ++i) {
    double r = y[i] - mu - u[i];
    g[i] += -r * is2 + l0 * u[i] + b * std::exp(u[i]);
  }
  for (int i = 1; i < m; ++i) {
    double d = u[i] - u[i - 1];
    g[i] += lrw * d;
    g[i - 1] -= lrw * d;
  }
  return g;
}
Eigen::VectorXd solve_u(const Eigen::VectorXd &t, int m) {
  Eigen::VectorXd u = Eigen::VectorXd::Zero(m);
  for (int it = 0; it < 80; ++it) {
    Eigen::LDLT<Eigen::MatrixXd> ldlt(Eigen::MatrixXd(Huu(t, u)));
    Eigen::VectorXd s = ldlt.solve(gu(t, u));
    u -= s;
    if (s.lpNorm<Eigen::Infinity>() < 1e-12)
      break;
  }
  return u;
}
Eigen::VectorXd futheta(const Eigen::VectorXd &t, const Eigen::VectorXd &u,
                        int j) {
  int m = u.size();
  double is2 = std::exp(-2 * t[1]), l0 = std::exp(t[2]), lrw = std::exp(t[3]),
         b = std::exp(t[4]);
  auto y = yvec(m);
  Eigen::VectorXd c = Eigen::VectorXd::Zero(m);
  if (j == 0) {
    c.array() = is2;
    return c;
  }
  if (j == 1) {
    for (int i = 0; i < m; ++i)
      c[i] = 2 * (y[i] - t[0] - u[i]) * is2;
    return c;
  }
  if (j == 2) {
    for (int i = 0; i < m; ++i)
      c[i] = l0 * u[i];
    return c;
  }
  if (j == 3) {
    for (int i = 1; i < m; ++i) {
      double d = u[i] - u[i - 1];
      c[i] += lrw * d;
      c[i - 1] -= lrw * d;
    }
    return c;
  }
  if (j == 4) {
    for (int i = 0; i < m; ++i)
      c[i] = b * std::exp(u[i]);
    return c;
  }
  return c;
}
template <class F> double mean_ms(F &&f, int reps) {
  auto a = Clock::now();
  for (int r = 0; r < reps; ++r) {
    volatile double z = f();
    (void)z;
  }
  auto b = Clock::now();
  return std::chrono::duration<double, std::milli>(b - a).count() / reps;
}
const char *nm(int j) {
  static const char *n[] = {"mu", "log_sigma", "log_lambda0", "log_lambda_rw",
                            "log_beta"};
  return n[j];
}
int main(int argc, char **argv) {
  int reps = argc > 1 ? std::stoi(argv[1]) : 20;
  std::vector<int> dims = {10, 25, 50, 100, 250, 500};
  std::cout << "Sparse RW1 per-direction Hdot benchmark\nreps per direction = "
            << reps << "\n\n"
            << std::setw(8) << "m" << std::setw(18) << "direction"
            << std::setw(16) << "Hdot ms" << std::setw(12) << "nnz"
            << std::setw(16) << "norm" << "\n"
            << std::scientific << std::setprecision(6);
  for (int m : dims) {
    Eigen::VectorXd th(5);
    th << 0.55, std::log(0.65), std::log(0.55), std::log(0.90), std::log(0.25);
    Eigen::VectorXd u = solve_u(th, m);
    auto H = Huu(th, u);
    quadra::laplace::SparseHuuFactorization fac(H);
    auto dir = [&](int j) { return -fac.solve(futheta(th, u, j)); };
    SparseRw1Objective obj{m};
    auto pattern = pat(m);
    for (int j = 0; j < 5; ++j) {
      double ms = mean_ms(
          [&]() {
            auto p = quadra::laplace::
                make_had_quadra_replay_reuse_lazy_implicit_hdot_provider(
                    obj, dir, 5, m, pattern, std::vector<int>{j}, 0.0);
            auto hs = p.compute_all_sparse(th, u);
            auto &Hdot = hs[(size_t)j];
            return double(Hdot.nonZeros()) + Hdot.norm();
          },
          reps);
      auto p = quadra::laplace::
          make_had_quadra_replay_reuse_lazy_implicit_hdot_provider(
              obj, dir, 5, m, pattern, std::vector<int>{j}, 0.0);
      auto hs = p.compute_all_sparse(th, u);
      auto &Hdot = hs[(size_t)j];
      std::cout << std::setw(8) << m << std::setw(18) << nm(j) << std::setw(16)
                << ms << std::setw(12) << Hdot.nonZeros() << std::setw(16)
                << Hdot.norm() << "\n";
    }
  }
  std::cout << "\nBenchmark complete.\n";
}
