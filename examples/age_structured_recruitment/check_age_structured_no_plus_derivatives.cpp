#define main quadra_age_structured_no_plus_analytic_benchmark_main
#include "benchmark_age_structured_no_plus_analytic_banded.cpp"
#undef main

#include <iomanip>
#include <iostream>

namespace {

double objective_x(const Data &data, const Parameters &par,
                   const Eigen::VectorXd &x) {
  return eval_all(data, par, x).objective;
}

Eigen::VectorXd fd_gradient_check(const Data &data, const Parameters &par,
                                  const Eigen::VectorXd &x) {
  Eigen::VectorXd g(x.size());

  for (int i = 0; i < x.size(); ++i) {
    const double h = 1e-6 * (1.0 + std::abs(x[i]));
    Eigen::VectorXd xp = x;
    Eigen::VectorXd xm = x;
    xp[i] += h;
    xm[i] -= h;

    g[i] =
        (objective_x(data, par, xp) - objective_x(data, par, xm)) / (2.0 * h);
  }

  return g;
}

Eigen::MatrixXd fd_hessian_check(const Data &data, const Parameters &par,
                                 const Eigen::VectorXd &x) {
  const int n = static_cast<int>(x.size());
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(n, n);
  const double f0 = objective_x(data, par, x);

  for (int i = 0; i < n; ++i) {
    const double hi = 1e-5 * (1.0 + std::abs(x[i]));

    Eigen::VectorXd xp = x;
    Eigen::VectorXd xm = x;
    xp[i] += hi;
    xm[i] -= hi;

    H(i, i) =
        (objective_x(data, par, xp) - 2.0 * f0 + objective_x(data, par, xm)) /
        (hi * hi);

    for (int j = i + 1; j < n; ++j) {
      const double hj = 1e-5 * (1.0 + std::abs(x[j]));

      Eigen::VectorXd xpp = x;
      Eigen::VectorXd xpm = x;
      Eigen::VectorXd xmp = x;
      Eigen::VectorXd xmm = x;

      xpp[i] += hi;
      xpp[j] += hj;
      xpm[i] += hi;
      xpm[j] -= hj;
      xmp[i] -= hi;
      xmp[j] += hj;
      xmm[i] -= hi;
      xmm[j] -= hj;

      const double hij =
          (objective_x(data, par, xpp) - objective_x(data, par, xpm) -
           objective_x(data, par, xmp) + objective_x(data, par, xmm)) /
          (4.0 * hi * hj);

      H(i, j) = hij;
      H(j, i) = hij;
    }
  }

  return H;
}

} // namespace

int main(int argc, char **argv) {
  int n = 25;
  int ages = 10;

  if (argc > 1)
    n = std::stoi(argv[1]);
  if (argc > 2)
    ages = std::stoi(argv[2]);

  const Parameters par;
  const Data data = make_data(n, ages, par);

  Eigen::VectorXd x = Eigen::VectorXd::Zero(n);

  // Put the check away from exactly zero too.
  for (int i = 0; i < n; ++i) {
    x[i] = 0.05 * std::sin(2.0 * M_PI * static_cast<double>(i) / 17.0);
  }

  const EvalAll a = eval_all(data, par, x);
  const Eigen::VectorXd g_fd = fd_gradient_check(data, par, x);

  const Eigen::VectorXd g_diff = a.gradient - g_fd;

  Eigen::MatrixXd H_analytic = Eigen::MatrixXd(a.hessian);
  Eigen::MatrixXd H_fd = fd_hessian_check(data, par, x);
  Eigen::MatrixXd H_diff = H_analytic - H_fd;

  double max_abs_H = 0.0;
  int max_i = 0;
  int max_j = 0;

  for (int i = 0; i < H_diff.rows(); ++i) {
    for (int j = 0; j < H_diff.cols(); ++j) {
      const double v = std::abs(H_diff(i, j));
      if (v > max_abs_H) {
        max_abs_H = v;
        max_i = i;
        max_j = j;
      }
    }
  }

  std::cout << std::fixed << std::setprecision(10);
  std::cout << "Age-structured no-plus derivative check\n";
  std::cout << "n = " << n << ", ages = " << ages << "\n\n";

  std::cout << "objective = " << a.objective << "\n";
  std::cout << "grad analytic norm = " << a.gradient.norm() << "\n";
  std::cout << "grad fd norm       = " << g_fd.norm() << "\n";
  std::cout << "grad max abs diff  = " << g_diff.cwiseAbs().maxCoeff() << "\n";
  std::cout << "grad rel diff      = "
            << g_diff.norm() / std::max(1.0, g_fd.norm()) << "\n\n";

  std::cout << "H analytic nnz     = " << a.hessian.nonZeros() << "\n";
  std::cout << "H fd norm          = " << H_fd.norm() << "\n";
  std::cout << "H diff norm        = " << H_diff.norm() << "\n";
  std::cout << "H rel diff         = "
            << H_diff.norm() / std::max(1.0, H_fd.norm()) << "\n";
  std::cout << "H max abs diff     = " << max_abs_H << " at (" << max_i << ","
            << max_j << ")\n\n";

  std::cout << "First 10 gradient entries:\n";
  std::cout << std::setw(6) << "i" << std::setw(18) << "analytic"
            << std::setw(18) << "fd" << std::setw(18) << "diff" << "\n";

  for (int i = 0; i < std::min(n, 10); ++i) {
    std::cout << std::setw(6) << i << std::setw(18) << a.gradient[i]
              << std::setw(18) << g_fd[i] << std::setw(18) << g_diff[i] << "\n";
  }

  std::cout << "\nLargest Hessian diff local row around max:\n";
  const int i0 = std::max(0, max_i - 3);
  const int i1 = std::min(n - 1, max_i + 3);
  const int j0 = std::max(0, max_j - 3);
  const int j1 = std::min(n - 1, max_j + 3);

  for (int i = i0; i <= i1; ++i) {
    for (int j = j0; j <= j1; ++j) {
      std::cout << "(" << i << "," << j << ")" << " A=" << H_analytic(i, j)
                << " FD=" << H_fd(i, j) << " D=" << H_diff(i, j) << "\n";
    }
  }

  return 0;
}
