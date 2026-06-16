#include "../data/pollock_data.hpp"
#include "reports/pollock_reports.hpp"
#include "diagnostics/pollock_utilities.hpp"
#include "model/pollock_parameters.hpp"
#include "model/pollock_constants.hpp"
#include "model/pollock_model.hpp"
#include "../../../../core/optimizer.hpp"
#include "../../../../core/laplace/laplace_structure_report.hpp"
#include "../../../../core/laplace/functional_analysis_report.hpp"
#include <Eigen/Dense>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{



std::vector<std::string> split(const std::string &line)
  {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string x;
    while (std::getline(ss, x, ','))
      out.push_back(x);
    return out;
  }

  std::vector<Obs> read_obs(const std::string &path)
  {
    std::ifstream in(path);
    if (!in)
      throw std::runtime_error("could not open " + path);
    std::string line;
    std::getline(in, line);
    std::vector<Obs> rows;
    while (std::getline(in, line))
    {
      if (line.empty())
        continue;
      auto f = split(line);
      Obs o{std::stoi(f[0]), std::stod(f[1]), std::stod(f[2]), {}};
      for (std::size_t i = 3; i < f.size(); ++i)
        o.age.push_back(std::stod(f[i]));
      rows.push_back(o);
    }
    return rows;
  }

  template <class AD>
  AD logistic(const AD &x)
  {
    return AD(1.0) / (AD(1.0) + exp(-x));
  }







void write_summary(const std::string &path, const quadra::OptResult &fit)
  {
    std::ofstream out(path);
    out << std::setprecision(15);
    out << "field,value\n";
    out << "objective," << fit.value << "\n";
    out << "grad_norm," << fit.grad_norm << "\n";
    out << "iterations," << fit.iterations << "\n";
    out << "converged," << (fit.converged ? "yes" : "no") << "\n";
    out << "message," << fit.message << "\n";
    out << "random_effects," << fit.u_hat.size() << "\n";
  }

#ifdef WALLEYE_POLLOCK_HUU_DIAGNOSTICS
  void pollock_write_huu_diagnostics(const std::string &path,
                                     PollockModel &model,
                                     quadra::ParameterVector &params,
                                     const quadra::OptResult &fit)
  {
    std::ofstream out(path);
    out << std::setprecision(15);
    out << "field,value\n";
    out << "random_effects," << fit.u_hat.size() << "\n";

    if (fit.u_hat.empty())
    {
      out << "available,no\n";
      out << "reason,no random effects\n";
      return;
    }

    try
    {
      const auto fixed_idx = quadra::build_fixed_index(params);
      const auto random_idx = quadra::build_random_index(params);

      for (std::size_t k = 0; k < fixed_idx.size() && k < fit.par.size(); ++k)
      {
        params.params[static_cast<std::size_t>(fixed_idx[k])].value = fit.par[k];
      }

      for (std::size_t k = 0; k < random_idx.size() && k < fit.u_hat.size(); ++k)
      {
        params.params[static_cast<std::size_t>(random_idx[k])].value = fit.u_hat[k];
      }

      had::ADGraph graph;
      quadra::ADScope scope(graph);

      std::vector<quadra::AD> p_full;
      p_full.reserve(static_cast<std::size_t>(params.size()));
      for (int i = 0; i < params.size(); ++i)
      {
        p_full.emplace_back(
            quadra::AD(params.params[static_cast<std::size_t>(i)].value));
      }

      quadra::AD nll = model(p_full);
      scope.backward(nll);

      const auto &pattern = quadra::get_pattern(scope, p_full, random_idx);
      Eigen::SparseMatrix<double> H =
          quadra::extract_sparse_hessian(scope, p_full, random_idx, pattern);

      Eigen::MatrixXd dense = Eigen::MatrixXd(H);
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(dense);

      out << "available,yes\n";
      out << "pattern_entries," << pattern.size() << "\n";
      out << "hessian_nonzeros," << H.nonZeros() << "\n";
      out << "min_diagonal," << dense.diagonal().minCoeff() << "\n";
      out << "max_diagonal," << dense.diagonal().maxCoeff() << "\n";

      if (es.info() == Eigen::Success)
      {
        const auto evals = es.eigenvalues();
        out << "eigen_success,yes\n";
        out << "min_eigenvalue," << evals.minCoeff() << "\n";
        out << "max_eigenvalue," << evals.maxCoeff() << "\n";
        out << "positive_definite," << (evals.minCoeff() > 0.0 ? "yes" : "no") << "\n";
        if (std::abs(evals.minCoeff()) > 0.0)
        {
          out << "condition_number_abs,"
              << std::abs(evals.maxCoeff()) / std::abs(evals.minCoeff()) << "\n";
        }
        else
        {
          out << "condition_number_abs,inf\n";
        }
      }
      else
      {
        out << "eigen_success,no\n";
      }
    }
    catch (const std::exception &e)
    {
      out << "available,no\n";
      out << "reason," << e.what() << "\n";
    }
  }
#endif


void write_fixed_gradient_diagnostics(const std::string &path,
                                      const quadra::OptResult &fit) {
  std::ofstream out(path);
  out << std::setprecision(15);
  out << "parameter,gradient,abs_gradient\n";

  for (std::size_t i = 0; i < fit.fixed_gradient.size(); ++i) {
    const std::string name =
        (i < fit.fixed_gradient_names.size()) ? fit.fixed_gradient_names[i]
                                              : ("fixed_" + std::to_string(i));
    const double g = fit.fixed_gradient[i];
    out << name << "," << g << "," << std::abs(g) << "\n";
  }
}

std::size_t max_fixed_gradient_index(const quadra::OptResult &fit) {
  std::size_t best = 0;
  double best_abs = -1.0;

  for (std::size_t i = 0; i < fit.fixed_gradient.size(); ++i) {
    const double a = std::abs(fit.fixed_gradient[i]);
    if (a > best_abs) {
      best = i;
      best_abs = a;
    }
  }

  return best;
}





void write_fixed_parameter_estimates(const std::string &path,
                                     const quadra::OptResult &fit) {
  std::ofstream out(path);
  out << std::setprecision(15);
  out << "parameter,estimate,exp_estimate\n";

  for (std::size_t i = 0; i < fit.par.size(); ++i) {
    const std::string name =
        (i < fit.fixed_gradient_names.size()) ? fit.fixed_gradient_names[i]
                                              : ("fixed_" + std::to_string(i));
    out << name << "," << fit.par[i] << "," << std::exp(fit.par[i]) << "\n";
  }
}



#ifdef WALLEYE_POLLOCK_FIXED_HESSIAN_DIAGNOSTICS
double pollock_profile_objective_at_fixed(
    PollockModel &model,
    quadra::ParameterVector params,
    const std::vector<int> &fixed_idx,
    const std::vector<int> &random_idx,
    const Eigen::VectorXd &x,
    const quadra::LaplaceOptions &opts) {
  for (std::size_t k = 0; k < fixed_idx.size(); ++k) {
    params.params[static_cast<std::size_t>(fixed_idx[k])].value =
        x[static_cast<Eigen::Index>(k)];
  }

  had::ADGraph graph;

  if (random_idx.empty()) {
    std::vector<double> p_double;
    p_double.reserve(static_cast<std::size_t>(params.size()));
    for (int i = 0; i < params.size(); ++i) {
      p_double.emplace_back(
          params.params[static_cast<std::size_t>(i)].value);
    }
    return model(p_double);
  }

  const auto u_hat = quadra::solve_random_effects_laplace(
      model, params, x, fixed_idx, random_idx, graph);
  const auto res = quadra::laplace_eval_at_u_star(
      model, params, fixed_idx, random_idx, x, u_hat, graph, opts);
  return res.value;
}

void write_fixed_hessian_diagnostics(
    const std::string &summary_path,
    const std::string &matrix_path,
    PollockModel &model,
    const quadra::ParameterVector &params_in,
    const quadra::OptResult &fit,
    const quadra::LaplaceOptions &opts) {
  std::ofstream summary(summary_path);
  summary << std::setprecision(15);
  summary << "field,value\n";

  quadra::ParameterVector params = params_in;
  const auto fixed_idx = quadra::build_fixed_index(params);
  const auto random_idx = quadra::build_random_index(params);

  const Eigen::Index n = static_cast<Eigen::Index>(fixed_idx.size());
  summary << "fixed_effects," << n << "\n";

  if (n == 0 || fit.par.size() != static_cast<std::size_t>(n)) {
    summary << "available,no\n";
    summary << "reason,missing fixed-effect vector\n";
    return;
  }

  try {
    Eigen::VectorXd x(n);
    for (Eigen::Index i = 0; i < n; ++i) {
      x[i] = fit.par[static_cast<std::size_t>(i)];
    }

    const double eps = 1.0e-4;
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(n, n);

    const double f0 = pollock_profile_objective_at_fixed(
        model, params, fixed_idx, random_idx, x, opts);

    for (Eigen::Index i = 0; i < n; ++i) {
      Eigen::VectorXd xp = x;
      Eigen::VectorXd xm = x;
      xp[i] += eps;
      xm[i] -= eps;

      const double fp = pollock_profile_objective_at_fixed(
          model, params, fixed_idx, random_idx, xp, opts);
      const double fm = pollock_profile_objective_at_fixed(
          model, params, fixed_idx, random_idx, xm, opts);

      H(i, i) = (fp - 2.0 * f0 + fm) / (eps * eps);

      for (Eigen::Index j = i + 1; j < n; ++j) {
        Eigen::VectorXd xpp = x;
        Eigen::VectorXd xpm = x;
        Eigen::VectorXd xmp = x;
        Eigen::VectorXd xmm = x;

        xpp[i] += eps;
        xpp[j] += eps;
        xpm[i] += eps;
        xpm[j] -= eps;
        xmp[i] -= eps;
        xmp[j] += eps;
        xmm[i] -= eps;
        xmm[j] -= eps;

        const double fpp = pollock_profile_objective_at_fixed(
            model, params, fixed_idx, random_idx, xpp, opts);
        const double fpm = pollock_profile_objective_at_fixed(
            model, params, fixed_idx, random_idx, xpm, opts);
        const double fmp = pollock_profile_objective_at_fixed(
            model, params, fixed_idx, random_idx, xmp, opts);
        const double fmm = pollock_profile_objective_at_fixed(
            model, params, fixed_idx, random_idx, xmm, opts);

        const double hij = (fpp - fpm - fmp + fmm) / (4.0 * eps * eps);
        H(i, j) = hij;
        H(j, i) = hij;
      }
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);

    summary << "available,yes\n";
    summary << "fd_step," << eps << "\n";
    summary << "profile_objective," << f0 << "\n";
    summary << "min_diagonal," << H.diagonal().minCoeff() << "\n";
    summary << "max_diagonal," << H.diagonal().maxCoeff() << "\n";

    if (es.info() == Eigen::Success) {
      const auto evals = es.eigenvalues();
      summary << "eigen_success,yes\n";
      summary << "min_eigenvalue," << evals.minCoeff() << "\n";
      summary << "max_eigenvalue," << evals.maxCoeff() << "\n";
      summary << "positive_definite,"
              << (evals.minCoeff() > 0.0 ? "yes" : "no") << "\n";

      if (std::abs(evals.minCoeff()) > 0.0) {
        summary << "condition_number_abs,"
                << std::abs(evals.maxCoeff()) / std::abs(evals.minCoeff())
                << "\n";
      } else {
        summary << "condition_number_abs,inf\n";
      }

      summary << "eigenvalues";
      for (Eigen::Index i = 0; i < evals.size(); ++i) {
        summary << (i == 0 ? "," : ";") << evals[i];
      }
      summary << "\n";
    } else {
      summary << "eigen_success,no\n";
    }

    std::ofstream mat(matrix_path);
    mat << "parameter";
    for (std::size_t j = 0; j < fit.fixed_gradient_names.size(); ++j) {
      mat << "," << fit.fixed_gradient_names[j];
    }
    mat << "\n";

    for (Eigen::Index i = 0; i < n; ++i) {
      const std::string row_name =
          (static_cast<std::size_t>(i) < fit.fixed_gradient_names.size())
              ? fit.fixed_gradient_names[static_cast<std::size_t>(i)]
              : ("fixed_" + std::to_string(i));
      mat << row_name;
      for (Eigen::Index j = 0; j < n; ++j) {
        mat << "," << std::setprecision(15) << H(i, j);
      }
      mat << "\n";
    }
  } catch (const std::exception &e) {
    summary << "available,no\n";
    summary << "reason," << e.what() << "\n";
  }
}
#endif




#ifdef WALLEYE_POLLOCK_HUU_MATRIX_DUMP
double pollock_joint_objective_at_x_u(
    PollockModel &model,
    quadra::ParameterVector params,
    const std::vector<int> &fixed_idx,
    const std::vector<int> &random_idx,
    const Eigen::VectorXd &x,
    const std::vector<double> &u) {
  for (std::size_t k = 0; k < fixed_idx.size(); ++k) {
    params.params[static_cast<std::size_t>(fixed_idx[k])].value =
        x[static_cast<Eigen::Index>(k)];
  }
  for (std::size_t k = 0; k < random_idx.size(); ++k) {
    params.params[static_cast<std::size_t>(random_idx[k])].value = u[k];
  }

  std::vector<double> p_double;
  p_double.reserve(static_cast<std::size_t>(params.size()));
  for (int i = 0; i < params.size(); ++i) {
    p_double.emplace_back(params.params[static_cast<std::size_t>(i)].value);
  }

  return model(p_double);
}

Eigen::MatrixXd pollock_fd_huu(
    PollockModel &model,
    quadra::ParameterVector params,
    const quadra::OptResult &fit,
    double eps = 1.0e-4) {
  const auto fixed_idx = quadra::build_fixed_index(params);
  const auto random_idx = quadra::build_random_index(params);

  const Eigen::Index n = static_cast<Eigen::Index>(random_idx.size());
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(n, n);

  if (n == 0 || fit.par.size() != fixed_idx.size() ||
      fit.u_hat.size() != random_idx.size()) {
    return H;
  }

  Eigen::VectorXd x(static_cast<Eigen::Index>(fixed_idx.size()));
  for (std::size_t i = 0; i < fixed_idx.size(); ++i) {
    x[static_cast<Eigen::Index>(i)] = fit.par[i];
  }

  std::vector<double> u = fit.u_hat;
  const double f0 = pollock_joint_objective_at_x_u(
      model, params, fixed_idx, random_idx, x, u);

  for (Eigen::Index i = 0; i < n; ++i) {
    std::vector<double> up = u;
    std::vector<double> um = u;
    up[static_cast<std::size_t>(i)] += eps;
    um[static_cast<std::size_t>(i)] -= eps;

    const double fp = pollock_joint_objective_at_x_u(
        model, params, fixed_idx, random_idx, x, up);
    const double fm = pollock_joint_objective_at_x_u(
        model, params, fixed_idx, random_idx, x, um);

    H(i, i) = (fp - 2.0 * f0 + fm) / (eps * eps);

    for (Eigen::Index j = i + 1; j < n; ++j) {
      std::vector<double> upp = u;
      std::vector<double> upm = u;
      std::vector<double> ump = u;
      std::vector<double> umm = u;

      upp[static_cast<std::size_t>(i)] += eps;
      upp[static_cast<std::size_t>(j)] += eps;

      upm[static_cast<std::size_t>(i)] += eps;
      upm[static_cast<std::size_t>(j)] -= eps;

      ump[static_cast<std::size_t>(i)] -= eps;
      ump[static_cast<std::size_t>(j)] += eps;

      umm[static_cast<std::size_t>(i)] -= eps;
      umm[static_cast<std::size_t>(j)] -= eps;

      const double fpp = pollock_joint_objective_at_x_u(
          model, params, fixed_idx, random_idx, x, upp);
      const double fpm = pollock_joint_objective_at_x_u(
          model, params, fixed_idx, random_idx, x, upm);
      const double fmp = pollock_joint_objective_at_x_u(
          model, params, fixed_idx, random_idx, x, ump);
      const double fmm = pollock_joint_objective_at_x_u(
          model, params, fixed_idx, random_idx, x, umm);

      const double hij = (fpp - fpm - fmp + fmm) / (4.0 * eps * eps);
      H(i, j) = hij;
      H(j, i) = hij;
    }
  }

  return H;
}

void pollock_write_huu_matrix(
    const std::string &path,
    PollockModel &model,
    quadra::ParameterVector params,
    const quadra::OptResult &fit) {
  std::ofstream out(path);
  out << std::setprecision(15);

  const auto random_idx = quadra::build_random_index(params);
  const std::size_t n = random_idx.size();

  out << "row";
  for (std::size_t j = 0; j < n; ++j) {
    out << ",u" << (j + 1);
  }
  out << "\n";

  Eigen::MatrixXd dense = pollock_fd_huu(model, params, fit);

  for (std::size_t i = 0; i < n; ++i) {
    out << "u" << (i + 1);
    for (std::size_t j = 0; j < n; ++j) {
      out << "," << dense(static_cast<Eigen::Index>(i),
                          static_cast<Eigen::Index>(j));
    }
    out << "\n";
  }
}

void pollock_write_huu_sparsity(
    const std::string &path,
    PollockModel &model,
    quadra::ParameterVector params,
    const quadra::OptResult &fit,
    double tol = 1.0e-8) {
  std::ofstream out(path);
  out << "i,j,value,abs_value\n";
  out << std::setprecision(15);

  Eigen::MatrixXd dense = pollock_fd_huu(model, params, fit);

  for (Eigen::Index i = 0; i < dense.rows(); ++i) {
    for (Eigen::Index j = 0; j < dense.cols(); ++j) {
      const double v = dense(i, j);
      if (std::abs(v) > tol) {
        out << (i + 1) << "," << (j + 1) << "," << v << ","
            << std::abs(v) << "\n";
      }
    }
  }
}
#endif





#ifdef WALLEYE_POLLOCK_HUU_PATTERN_COMPARE
void pollock_write_huu_pattern_compare(
    const std::string &path,
    PollockModel &model,
    quadra::ParameterVector params,
    const quadra::OptResult &fit,
    double tol = 1.0e-8) {
  std::ofstream out(path);
  out << "field,value\n";

  const auto fixed_idx = quadra::build_fixed_index(params);
  const auto random_idx = quadra::build_random_index(params);
  const std::size_t n = random_idx.size();

  out << "random_effects," << n << "\n";
  out << "fd_tol," << tol << "\n";
  out << "quadra_pattern_available," << (fit.pattern.available ? "yes" : "no") << "\n";
  out << "quadra_pattern_detected_structure," << fit.pattern.detected_structure << "\n";
  out << "quadra_pattern_nonzeros_reported," << fit.pattern.nonzeros << "\n";

  if (n == 0 || fit.par.size() != fixed_idx.size() ||
      fit.u_hat.size() != n) {
    out << "available,no\n";
    out << "reason,missing random effects or size mismatch\n";
    return;
  }

  Eigen::MatrixXd Hfd = pollock_fd_huu(model, params, fit);

  std::size_t fd_nonzeros_all = 0;
  std::size_t fd_nonzeros_upper = 0;
  std::size_t fd_nonzeros_diag = 0;
  double max_abs_fd = 0.0;
  double min_abs_fd_nonzero = std::numeric_limits<double>::infinity();

  for (Eigen::Index i = 0; i < Hfd.rows(); ++i) {
    for (Eigen::Index j = 0; j < Hfd.cols(); ++j) {
      const double av = std::abs(Hfd(i, j));
      max_abs_fd = std::max(max_abs_fd, av);
      if (av > tol) {
        ++fd_nonzeros_all;
        min_abs_fd_nonzero = std::min(min_abs_fd_nonzero, av);
        if (i <= j) {
          ++fd_nonzeros_upper;
        }
        if (i == j) {
          ++fd_nonzeros_diag;
        }
      }
    }
  }

  const std::size_t fd_nonzeros_offdiag_all =
      fd_nonzeros_all >= fd_nonzeros_diag
          ? fd_nonzeros_all - fd_nonzeros_diag
          : 0;
  const std::size_t fd_nonzeros_offdiag_upper =
      fd_nonzeros_upper >= fd_nonzeros_diag
          ? fd_nonzeros_upper - fd_nonzeros_diag
          : 0;

  out << "available,yes\n";
  out << "fd_nonzeros_all," << fd_nonzeros_all << "\n";
  out << "fd_nonzeros_upper_including_diag," << fd_nonzeros_upper << "\n";
  out << "fd_nonzeros_diag," << fd_nonzeros_diag << "\n";
  out << "fd_nonzeros_offdiag_all," << fd_nonzeros_offdiag_all << "\n";
  out << "fd_nonzeros_offdiag_upper," << fd_nonzeros_offdiag_upper << "\n";
  out << "fd_density_all," << (n == 0 ? 0.0 : static_cast<double>(fd_nonzeros_all) / static_cast<double>(n * n)) << "\n";
  out << "fd_density_upper," << (n == 0 ? 0.0 : static_cast<double>(fd_nonzeros_upper) / static_cast<double>((n * (n + 1)) / 2)) << "\n";
  out << "max_abs_fd," << max_abs_fd << "\n";
  out << "min_abs_fd_nonzero,"
      << (std::isfinite(min_abs_fd_nonzero) ? min_abs_fd_nonzero : 0.0)
      << "\n";
  out << "note,OptPatternInfo does not currently expose individual pattern entries; this compares reported Quadra count to finite-difference numerical sparsity.\n";

  std::ofstream detail(
      "examples/NMFS/afsc_walleye_pollock/outputs/"
      "walleye_pollock_huu_pattern_compare_detail.csv");
  detail << "i,j,fd_nonzero,fd_value,abs_fd_value,band_distance\n";
  detail << std::setprecision(15);

  for (Eigen::Index i = 0; i < Hfd.rows(); ++i) {
    for (Eigen::Index j = 0; j < Hfd.cols(); ++j) {
      const double v = Hfd(i, j);
      const double av = std::abs(v);
      if (av > tol) {
        detail << (i + 1) << "," << (j + 1) << ",yes,"
               << v << "," << av << "," << std::abs(i - j) << "\n";
      }
    }
  }
}
#endif



#ifdef WALLEYE_POLLOCK_HUU_BAND_SUMMARY
void pollock_write_huu_band_summary(
    const std::string &path,
    PollockModel &model,
    quadra::ParameterVector params,
    const quadra::OptResult &fit,
    double tol = 1.0e-8) {
  Eigen::MatrixXd H = pollock_fd_huu(model, params, fit);

  std::ofstream out(path);
  out << "band_distance,count,nonzero_count,mean_abs,max_abs,sum_abs,share_sum_abs,cumulative_share_sum_abs\n";
  out << std::setprecision(15);

  if (H.rows() == 0) {
    return;
  }

  const Eigen::Index n = H.rows();
  std::vector<double> sum_abs(static_cast<std::size_t>(n), 0.0);
  std::vector<double> max_abs(static_cast<std::size_t>(n), 0.0);
  std::vector<std::size_t> count(static_cast<std::size_t>(n), 0);
  std::vector<std::size_t> nonzero_count(static_cast<std::size_t>(n), 0);

  double total_abs = 0.0;

  // Use upper triangle including diagonal so each symmetric pair is counted once.
  for (Eigen::Index i = 0; i < n; ++i) {
    for (Eigen::Index j = i; j < n; ++j) {
      const std::size_t d = static_cast<std::size_t>(j - i);
      const double av = std::abs(H(i, j));

      count[d] += 1;
      sum_abs[d] += av;
      max_abs[d] = std::max(max_abs[d], av);
      total_abs += av;

      if (av > tol) {
        nonzero_count[d] += 1;
      }
    }
  }

  double cumulative = 0.0;
  for (std::size_t d = 0; d < static_cast<std::size_t>(n); ++d) {
    const double mean_abs =
        count[d] > 0 ? sum_abs[d] / static_cast<double>(count[d]) : 0.0;
    const double share =
        total_abs > 0.0 ? sum_abs[d] / total_abs : 0.0;
    cumulative += share;

    out << d << ","
        << count[d] << ","
        << nonzero_count[d] << ","
        << mean_abs << ","
        << max_abs[d] << ","
        << sum_abs[d] << ","
        << share << ","
        << cumulative << "\n";
  }
}
#endif



#ifdef WALLEYE_POLLOCK_HUU_BANDLIMIT_DIAGNOSTIC
void pollock_write_huu_bandlimit_diagnostic(
    const std::string &path,
    PollockModel &model,
    quadra::ParameterVector params,
    const quadra::OptResult &fit) {
  Eigen::MatrixXd H = pollock_fd_huu(model, params, fit);

  std::ofstream out(path);
  out << "band_width,kept_entries,total_entries,kept_entry_share,"
         "retained_abs_share,relative_frobenius_error,"
         "min_eigenvalue,max_eigenvalue,positive_definite,condition_number_abs\n";
  out << std::setprecision(15);

  if (H.rows() == 0) {
    return;
  }

  const Eigen::Index n = H.rows();
  const double full_abs_sum = H.cwiseAbs().sum();
  const double full_frob = H.norm();

  const std::vector<int> bands = {0, 1, 2, 3, 5, 10, 20};

  for (const int bw_raw : bands) {
    const Eigen::Index bw = std::min<Eigen::Index>(
        static_cast<Eigen::Index>(bw_raw), n - 1);

    Eigen::MatrixXd Hb = Eigen::MatrixXd::Zero(n, n);
    std::size_t kept_entries = 0;

    for (Eigen::Index i = 0; i < n; ++i) {
      for (Eigen::Index j = 0; j < n; ++j) {
        if (std::abs(i - j) <= bw) {
          Hb(i, j) = H(i, j);
          ++kept_entries;
        }
      }
    }

    const double retained_abs_share =
        full_abs_sum > 0.0 ? Hb.cwiseAbs().sum() / full_abs_sum : 0.0;
    const double rel_frob_error =
        full_frob > 0.0 ? (H - Hb).norm() / full_frob : 0.0;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(Hb);
    const bool eig_ok = eig.info() == Eigen::Success;
    double min_eval = std::numeric_limits<double>::quiet_NaN();
    double max_eval = std::numeric_limits<double>::quiet_NaN();
    bool pd = false;
    double cond = std::numeric_limits<double>::quiet_NaN();

    if (eig_ok && eig.eigenvalues().size() > 0) {
      min_eval = eig.eigenvalues().minCoeff();
      max_eval = eig.eigenvalues().maxCoeff();
      pd = min_eval > 0.0;
      cond = std::abs(max_eval) / std::max(std::abs(min_eval), 1.0e-300);
    }

    out << bw << ","
        << kept_entries << ","
        << static_cast<std::size_t>(n * n) << ","
        << static_cast<double>(kept_entries) / static_cast<double>(n * n) << ","
        << retained_abs_share << ","
        << rel_frob_error << ","
        << min_eval << ","
        << max_eval << ","
        << (pd ? "yes" : "no") << ","
        << cond << "\n";
  }
}
#endif



#ifdef WALLEYE_POLLOCK_HUU_THRESHOLD_DIAGNOSTIC
void pollock_write_huu_threshold_diagnostic(
    const std::string &path,
    PollockModel &model,
    quadra::ParameterVector params,
    const quadra::OptResult &fit) {
  Eigen::MatrixXd H = pollock_fd_huu(model, params, fit);

  std::ofstream out(path);
  out << "threshold_type,threshold,absolute_threshold,"
         "kept_entries,total_entries,kept_entry_share,"
         "retained_abs_share,relative_frobenius_error,"
         "min_eigenvalue,max_eigenvalue,positive_definite,condition_number_abs\n";
  out << std::setprecision(15);

  if (H.rows() == 0) {
    return;
  }

  const Eigen::Index n = H.rows();
  const double full_abs_sum = H.cwiseAbs().sum();
  const double full_frob = H.norm();
  const double max_abs = H.cwiseAbs().maxCoeff();

  struct ThresholdSpec {
    const char *type;
    double threshold;
    double absolute_threshold;
  };

  std::vector<ThresholdSpec> specs;

  for (double t : {1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6}) {
    specs.push_back({"absolute", t, t});
  }

  for (double r : {1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5}) {
    specs.push_back({"relative_to_max_abs", r, r * max_abs});
  }

  for (const auto &spec : specs) {
    Eigen::MatrixXd Ht = Eigen::MatrixXd::Zero(n, n);
    std::size_t kept_entries = 0;

    for (Eigen::Index i = 0; i < n; ++i) {
      for (Eigen::Index j = 0; j < n; ++j) {
        const double v = H(i, j);
        if (std::abs(v) >= spec.absolute_threshold) {
          Ht(i, j) = v;
          ++kept_entries;
        }
      }
    }

    const double retained_abs_share =
        full_abs_sum > 0.0 ? Ht.cwiseAbs().sum() / full_abs_sum : 0.0;
    const double rel_frob_error =
        full_frob > 0.0 ? (H - Ht).norm() / full_frob : 0.0;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(Ht);
    const bool eig_ok = eig.info() == Eigen::Success;
    double min_eval = std::numeric_limits<double>::quiet_NaN();
    double max_eval = std::numeric_limits<double>::quiet_NaN();
    bool pd = false;
    double cond = std::numeric_limits<double>::quiet_NaN();

    if (eig_ok && eig.eigenvalues().size() > 0) {
      min_eval = eig.eigenvalues().minCoeff();
      max_eval = eig.eigenvalues().maxCoeff();
      pd = min_eval > 0.0;
      cond = std::abs(max_eval) / std::max(std::abs(min_eval), 1.0e-300);
    }

    out << spec.type << ","
        << spec.threshold << ","
        << spec.absolute_threshold << ","
        << kept_entries << ","
        << static_cast<std::size_t>(n * n) << ","
        << static_cast<double>(kept_entries) / static_cast<double>(n * n) << ","
        << retained_abs_share << ","
        << rel_frob_error << ","
        << min_eval << ","
        << max_eval << ","
        << (pd ? "yes" : "no") << ","
        << cond << "\n";
  }
}
#endif




#ifdef WALLEYE_POLLOCK_LAPLACE_STRUCTURE_REPORT
void pollock_write_laplace_structure_report(
    const std::string &path,
    PollockModel &model,
    quadra::ParameterVector params,
    const quadra::OptResult &fit,
    double nonzero_tol = 1.0e-8) {
  const Eigen::MatrixXd H = pollock_fd_huu(model, params, fit);
  const auto report =
      quadra::summarize_laplace_hessian_structure(H, nonzero_tol);

  quadra::write_laplace_structure_report_text(report, path);
  quadra::write_laplace_structure_report_csv(
      report,
      "examples/NMFS/afsc_walleye_pollock/outputs/"
      "walleye_pollock_laplace_structure_report.csv");
}
#endif




#ifdef WALLEYE_POLLOCK_GRADIENT_VOLATILITY
std::vector<int> pollock_fixed_indices(const quadra::ParameterVector &params) {
  std::vector<int> out;
  for (std::size_t i = 0; i < params.params.size(); ++i) {
    if (!params.params[i].is_random) {
      out.push_back(static_cast<int>(i));
    }
  }
  return out;
}

std::vector<int> pollock_random_indices(const quadra::ParameterVector &params) {
  std::vector<int> out;
  for (std::size_t i = 0; i < params.params.size(); ++i) {
    if (params.params[i].is_random) {
      out.push_back(static_cast<int>(i));
    }
  }
  return out;
}

Eigen::VectorXd pollock_fixed_values(const quadra::ParameterVector &params) {
  const auto fixed_idx = pollock_fixed_indices(params);
  Eigen::VectorXd x(static_cast<Eigen::Index>(fixed_idx.size()));
  for (std::size_t i = 0; i < fixed_idx.size(); ++i) {
    x(static_cast<Eigen::Index>(i)) =
        params.params[static_cast<std::size_t>(fixed_idx[i])].value;
  }
  return x;
}

std::vector<double> pollock_profile_gradient_fd_at_x(
    PollockModel &model,
    quadra::ParameterVector params,
    const Eigen::VectorXd &x,
    const std::vector<double> &u_hat,
    const quadra::LaplaceOptions &opts,
    double fd_step = 1.0e-5) {
  const auto fixed_idx = pollock_fixed_indices(params);
  const auto random_idx = pollock_random_indices(params);

  std::vector<double> grad(static_cast<std::size_t>(x.size()), 0.0);

  auto eval = [&](const Eigen::VectorXd &x_eval) -> double {
    had::ADGraph graph;
    auto res = quadra::laplace_eval_at_u_star(
        model, params, fixed_idx, random_idx, x_eval, u_hat, graph, opts);
    return res.value;
  };

  for (Eigen::Index j = 0; j < x.size(); ++j) {
    const double h = fd_step * std::max(1.0, std::abs(x(j)));
    Eigen::VectorXd xp = x;
    Eigen::VectorXd xm = x;
    xp(j) += h;
    xm(j) -= h;
    grad[static_cast<std::size_t>(j)] = (eval(xp) - eval(xm)) / (2.0 * h);
  }

  return grad;
}

quadra::FunctionalGradientVolatilitySummary
pollock_compute_gradient_volatility_fd(
    PollockModel &model,
    quadra::ParameterVector params,
    const quadra::OptResult &fit,
    double perturbation_scale = 1.0e-5,
    double fd_step = 1.0e-5) {
  quadra::FunctionalGradientVolatilitySummary empty;
  if (fit.fixed_gradient.empty()) {
    return empty;
  }

  Eigen::VectorXd x0 = fit.x;
  if (x0.size() == 0) {
    x0 = pollock_fixed_values(params);
  }

  if (x0.size() == 0 ||
      x0.size() != static_cast<Eigen::Index>(fit.fixed_gradient.size())) {
    return empty;
  }

  quadra::LaplaceOptions opts = quadra::default_laplace_options();

  std::vector<std::vector<double>> gradient_samples;
  gradient_samples.reserve(static_cast<std::size_t>(x0.size() * 2));

  for (Eigen::Index j = 0; j < x0.size(); ++j) {
    const double dx =
        perturbation_scale * std::max(1.0, std::abs(x0(j)));

    for (double sign : {-1.0, 1.0}) {
      Eigen::VectorXd xp = x0;
      xp(j) += sign * dx;
      gradient_samples.push_back(
          pollock_profile_gradient_fd_at_x(
              model, params, xp, fit.u_hat, opts, fd_step));
    }
  }

  return quadra::summarize_gradient_volatility(
      gradient_samples, fit.fixed_gradient, fit.fixed_gradient_names,
      perturbation_scale);
}
#endif


#ifdef WALLEYE_POLLOCK_PARAMETER_GEOMETRY
std::vector<int> pollock_parameter_geometry_fixed_indices(
    const quadra::ParameterVector &params) {
  std::vector<int> out;
  for (std::size_t i = 0; i < params.params.size(); ++i) {
    if (!params.params[i].is_random) {
      out.push_back(static_cast<int>(i));
    }
  }
  return out;
}

std::vector<int> pollock_parameter_geometry_random_indices(
    const quadra::ParameterVector &params) {
  std::vector<int> out;
  for (std::size_t i = 0; i < params.params.size(); ++i) {
    if (params.params[i].is_random) {
      out.push_back(static_cast<int>(i));
    }
  }
  return out;
}

Eigen::MatrixXd pollock_parameter_geometry_fd_fixed_hessian(
    PollockModel &model,
    quadra::ParameterVector params,
    const quadra::OptResult &fit,
    const quadra::LaplaceOptions &opts,
    double fd_step = 1.0e-4) {
  Eigen::VectorXd x0 = fit.x;
  if (x0.size() == 0) {
    // Backward-compatible fallback. Prefer fit.x when available.
    const auto fixed_idx = pollock_parameter_geometry_fixed_indices(params);
    x0.resize(static_cast<Eigen::Index>(fixed_idx.size()));
    for (std::size_t i = 0; i < fixed_idx.size(); ++i) {
      x0(static_cast<Eigen::Index>(i)) =
          params.params[static_cast<std::size_t>(fixed_idx[i])].value;
    }
  }

  const Eigen::Index n = x0.size();
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(n, n);

  if (n == 0) {
    return H;
  }

  const auto fixed_idx = pollock_parameter_geometry_fixed_indices(params);
  const auto random_idx = pollock_parameter_geometry_random_indices(params);

  auto eval = [&](const Eigen::VectorXd &x_eval) -> double {
    had::ADGraph graph;
    auto res = quadra::laplace_eval_at_u_star(
        model, params, fixed_idx, random_idx, x_eval, fit.u_hat, graph, opts);
    return res.value;
  };

  for (Eigen::Index i = 0; i < n; ++i) {
    const double hi = fd_step * std::max(1.0, std::abs(x0(i)));

    // Diagonal second derivative.
    {
      Eigen::VectorXd xp = x0;
      Eigen::VectorXd xm = x0;
      xp(i) += hi;
      xm(i) -= hi;

      const double f0 = eval(x0);
      const double fp = eval(xp);
      const double fm = eval(xm);
      H(i, i) = (fp - 2.0 * f0 + fm) / (hi * hi);
    }

    // Mixed second derivatives.
    for (Eigen::Index j = i + 1; j < n; ++j) {
      const double hj = fd_step * std::max(1.0, std::abs(x0(j)));

      Eigen::VectorXd xpp = x0;
      Eigen::VectorXd xpm = x0;
      Eigen::VectorXd xmp = x0;
      Eigen::VectorXd xmm = x0;

      xpp(i) += hi; xpp(j) += hj;
      xpm(i) += hi; xpm(j) -= hj;
      xmp(i) -= hi; xmp(j) += hj;
      xmm(i) -= hi; xmm(j) -= hj;

      const double fpp = eval(xpp);
      const double fpm = eval(xpm);
      const double fmp = eval(xmp);
      const double fmm = eval(xmm);

      const double hij = (fpp - fpm - fmp + fmm) / (4.0 * hi * hj);
      H(i, j) = hij;
      H(j, i) = hij;
    }
  }

  return H;
}
#endif

#ifdef WALLEYE_POLLOCK_FUNCTIONAL_ANALYSIS_REPORT
void pollock_write_functional_analysis_report(
    const std::string &text_path,
    const std::string &csv_path,
    PollockModel &model,
    quadra::ParameterVector params,
    const quadra::OptResult &fit,
    double nonzero_tol = 1.0e-8) {
  const Eigen::MatrixXd H = pollock_fd_huu(model, params, fit);

  quadra::FunctionalOptimizationSummary opt;
  opt.objective_value = fit.value;
  opt.gradient_norm = fit.grad_norm;
  opt.iterations = fit.iterations;
  opt.converged = fit.converged;
  opt.message = fit.message;

  if (!fit.fixed_gradient.empty()) {
    const std::size_t max_i = max_fixed_gradient_index(fit);
    opt.max_gradient_parameter =
        (max_i < fit.fixed_gradient_names.size())
            ? fit.fixed_gradient_names[max_i]
            : ("fixed_" + std::to_string(max_i));
    opt.max_gradient_value = fit.fixed_gradient[max_i];
    opt.max_abs_gradient = std::abs(fit.fixed_gradient[max_i]);
  }

  std::vector<std::string> random_names;
  random_names.reserve(fit.u_hat.size());
  for (std::size_t i = 0; i < fit.u_hat.size(); ++i) {
    random_names.push_back("log_rec_dev_" + std::to_string(i + 1));
  }

  auto report =
      quadra::make_functional_analysis_report(
          opt, H, fit.u_hat, nonzero_tol, random_names);

#ifdef WALLEYE_POLLOCK_PARAMETER_GEOMETRY
  {
    quadra::LaplaceOptions hess_opts = quadra::default_laplace_options();
    const Eigen::MatrixXd Hxx =
        pollock_parameter_geometry_fd_fixed_hessian(model, params, fit, hess_opts);

    report.parameter_geometry =
        quadra::summarize_parameter_geometry(
            Hxx, fit.fixed_gradient, fit.fixed_gradient_names);
  }
#endif

#ifdef WALLEYE_POLLOCK_GRADIENT_VOLATILITY
  {
    report.gradient_volatility =
        pollock_compute_gradient_volatility_fd(
            model, params, fit, 1.0e-5, 1.0e-5);
  }
#endif

  quadra::write_functional_analysis_report_text(report, text_path);
  quadra::write_functional_analysis_report_csv(report, csv_path);
}
#endif


} // namespace

int main()
{
  try
  {
    std::cout << "Synthetic AFSC walleye-pollock-style assessment example\n";
    std::cout << "=======================================================\n\n";
    std::cout << "Synthetic and public-data-safe. Not an official assessment.\n";
  std::cout << "Assessment-scale diagnostic: tolerance is relaxed for synthetic profiling/identifiability checks.\n";
  std::cout << "Recruitment deviations use a fixed AR(1) prior: rho=0.60, sigma=0.15.\n";
#ifdef WALLEYE_POLLOCK_RANDOM_RECRUITMENT_COUNT
    std::cout << "Random recruitment enabled for first "
              << WALLEYE_POLLOCK_RANDOM_RECRUITMENT_COUNT
              << " year(s).\n\n";
#else
    std::cout << "Level 1: fixed-effect index fit with observed-catch removals; random recruitment disabled.\n\n";
#endif

    auto obs = read_obs("examples/NMFS/afsc_walleye_pollock/data/synthetic_walleye_pollock_observations.csv");
    std::cout << "Loaded synthetic rows: " << obs.size() << "\n\n";

    PollockModel model(obs);
    auto params = pollock::make_params(obs.size());
    auto opts = quadra::default_laplace_options();

    auto fit = quadra::optimize_lbfgs(model, params, opts);

    write_summary("examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_fit_summary.csv", fit);
    write_fixed_parameter_estimates(
        "examples/NMFS/afsc_walleye_pollock/outputs/"
        "walleye_pollock_fixed_parameter_estimates.csv",
        fit);
    write_fixed_gradient_diagnostics(
        "examples/NMFS/afsc_walleye_pollock/outputs/"
        "walleye_pollock_fixed_gradient_diagnostics.csv",
        fit);

#ifdef WALLEYE_POLLOCK_FIXED_HESSIAN_DIAGNOSTICS
    {
      quadra::LaplaceOptions hess_opts = quadra::default_laplace_options();
      write_fixed_hessian_diagnostics(
          "examples/NMFS/afsc_walleye_pollock/outputs/"
          "walleye_pollock_fixed_hessian_diagnostics.csv",
          "examples/NMFS/afsc_walleye_pollock/outputs/"
          "walleye_pollock_fixed_hessian_matrix.csv",
          model, params, fit, hess_opts);
    }
#endif


#ifdef WALLEYE_POLLOCK_HUU_DIAGNOSTICS
    pollock_write_huu_diagnostics(
        "examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_huu_diagnostics.csv",
        model, params, fit);
#endif

#ifdef WALLEYE_POLLOCK_HUU_MATRIX_DUMP
    pollock_write_huu_matrix(
        "examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_huu_matrix.csv",
        model, params, fit);
    pollock_write_huu_sparsity(
        "examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_huu_sparsity.csv",
        model, params, fit);
#endif

#ifdef WALLEYE_POLLOCK_HUU_PATTERN_COMPARE
    pollock_write_huu_pattern_compare(
        "examples/NMFS/afsc_walleye_pollock/outputs/"
        "walleye_pollock_huu_pattern_compare.csv",
        model, params, fit);
#endif

#ifdef WALLEYE_POLLOCK_HUU_BAND_SUMMARY
    pollock_write_huu_band_summary(
        "examples/NMFS/afsc_walleye_pollock/outputs/"
        "walleye_pollock_huu_band_summary.csv",
        model, params, fit);
#endif

#ifdef WALLEYE_POLLOCK_HUU_BANDLIMIT_DIAGNOSTIC
    pollock_write_huu_bandlimit_diagnostic(
        "examples/NMFS/afsc_walleye_pollock/outputs/"
        "walleye_pollock_huu_bandlimit_diagnostic.csv",
        model, params, fit);
#endif

#ifdef WALLEYE_POLLOCK_HUU_THRESHOLD_DIAGNOSTIC
    pollock_write_huu_threshold_diagnostic(
        "examples/NMFS/afsc_walleye_pollock/outputs/"
        "walleye_pollock_huu_threshold_diagnostic.csv",
        model, params, fit);
#endif

#ifdef WALLEYE_POLLOCK_LAPLACE_STRUCTURE_REPORT
    pollock_write_laplace_structure_report(
        "examples/NMFS/afsc_walleye_pollock/outputs/"
        "walleye_pollock_laplace_structure_report.txt",
        model, params, fit);
#endif

#ifdef WALLEYE_POLLOCK_FUNCTIONAL_ANALYSIS_REPORT
    pollock_write_functional_analysis_report(
        "examples/NMFS/afsc_walleye_pollock/outputs/"
        "walleye_pollock_functional_analysis_report.txt",
        "examples/NMFS/afsc_walleye_pollock/outputs/"
        "walleye_pollock_functional_analysis_report.csv",
        model, params, fit);
#endif

#ifdef WALLEYE_POLLOCK_MARKDOWN_REPORT
    pollock_example::write_pollock_markdown_report(
        "examples/NMFS/afsc_walleye_pollock/outputs/"
        "walleye_pollock_analysis.md",
        "examples/NMFS/afsc_walleye_pollock/outputs/"
        "walleye_pollock_functional_analysis_report.csv",
        "examples/NMFS/afsc_walleye_pollock/outputs/"
        "walleye_pollock_laplace_structure_report.txt");
#endif

    std::ofstream rec("examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_recruitment_deviations.csv");
    rec << "year,log_rec_dev,ar1_rho,innovation\n";
    const double rec_rho_report = 0.60;
    for (std::size_t i = 0; i < fit.u_hat.size(); ++i)
    {
      const double innovation =
          (i == 0) ? fit.u_hat[i]
                   : (fit.u_hat[i] - rec_rho_report * fit.u_hat[i - 1]);
      rec << (i + 1) << "," << fit.u_hat[i] << ","
          << rec_rho_report << "," << innovation << "\n";
    }

    std::cout << "\nFit diagnostics\n";
    std::cout << "---------------\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "objective          " << fit.value << "\n";
    std::cout << "grad_norm          " << fit.grad_norm << "\n";
    std::cout << "iterations         " << fit.iterations << "\n";
    std::cout << "converged          " << (fit.converged ? "yes" : "no") << "\n";
    std::cout << "message            " << fit.message << "\n";
    if (!fit.fixed_gradient.empty()) {
      const std::size_t max_grad_i = max_fixed_gradient_index(fit);
      const std::string max_grad_name =
          (max_grad_i < fit.fixed_gradient_names.size())
              ? fit.fixed_gradient_names[max_grad_i]
              : ("fixed_" + std::to_string(max_grad_i));
      std::cout << "max_grad_param     " << max_grad_name << "\n";
      std::cout << "max_grad_value     " << fit.fixed_gradient[max_grad_i] << "\n";
      std::cout << "max_abs_grad       "
                << std::abs(fit.fixed_gradient[max_grad_i]) << "\n";
    }

    std::cout << "\nOptimizer structure diagnostics\n";
    std::cout << "-------------------------------\n";
    std::cout << "random effects     " << fit.pattern.random_effect_count << "\n";
    std::cout << "pattern available  " << (fit.pattern.available ? "yes" : "no") << "\n";
    std::cout << "detected structure " << fit.pattern.detected_structure << "\n";
    std::cout << "Hessian nonzeros   " << fit.pattern.nonzeros << "\n";

    std::cout << "\nWrote outputs:\n";
    std::cout << "  examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_fit_summary.csv\n";
    std::cout << "  examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_recruitment_deviations.csv\n";

    return fit.converged ? 0 : 2;
  }
  catch (const std::exception &e)
  {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 1;
  }
}
