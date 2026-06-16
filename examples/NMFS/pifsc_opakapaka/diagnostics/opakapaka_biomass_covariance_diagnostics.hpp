#pragma once

#include "../quadra/opakapaka_model.hpp"

#include "../../../../core/uncertainty/selected_inverse_diagonal.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

namespace opakapaka_example {

inline Eigen::MatrixXd compute_log_b_covariance_submatrix(
  const std::vector<opakapaka_example::Observation> &data,
  const std::vector<double> &u_hat, const Eigen::SparseMatrix<double> &h_uu)
{
const std::size_t n = std::min(data.size(), u_hat.size());
if (n == 0)
{
  return Eigen::MatrixXd();
}

std::vector<int> indices;
indices.reserve(n);
for (std::size_t i = 0; i < n; ++i)
{
  indices.push_back(static_cast<int>(i));
}

const auto log_b_cov =
    quadra::uncertainty::selected_inverse_submatrix_from_spd_hessian(h_uu,
                                                                     indices);

if (!log_b_cov.success)
{
  return Eigen::MatrixXd::Constant(static_cast<Eigen::Index>(n),
                                   static_cast<Eigen::Index>(n),
                                   std::numeric_limits<double>::quiet_NaN());
}

return log_b_cov.covariance;
}


inline Eigen::MatrixXd
log_cov_to_biomass_cov(const Eigen::MatrixXd &log_b_cov,
                     const std::vector<double> &u_hat)
{
const Eigen::Index n = log_b_cov.rows();
Eigen::MatrixXd biomass_cov =
    Eigen::MatrixXd::Constant(n, n, std::numeric_limits<double>::quiet_NaN());

for (Eigen::Index i = 0; i < n; ++i)
{
  const double b_i = std::exp(u_hat[static_cast<std::size_t>(i)]);
  for (Eigen::Index j = 0; j < n; ++j)
  {
    const double b_j = std::exp(u_hat[static_cast<std::size_t>(j)]);
    biomass_cov(i, j) = b_i * b_j * log_b_cov(i, j);
  }
}

return biomass_cov;
}


inline Eigen::MatrixXd covariance_to_correlation(const Eigen::MatrixXd &cov)
{
const Eigen::Index n = cov.rows();
Eigen::MatrixXd corr =
    Eigen::MatrixXd::Constant(n, n, std::numeric_limits<double>::quiet_NaN());

for (Eigen::Index i = 0; i < n; ++i)
{
  for (Eigen::Index j = 0; j < n; ++j)
  {
    const double vii = cov(i, i);
    const double vjj = cov(j, j);
    const double vij = cov(i, j);

    if (std::isfinite(vii) && std::isfinite(vjj) && std::isfinite(vij) &&
        vii > 0.0 && vjj > 0.0)
    {
      double c = vij / std::sqrt(vii * vjj);
      if (c > 1.0 && c < 1.0 + 1.0e-10)
        c = 1.0;
      if (c < -1.0 && c > -1.0 - 1.0e-10)
        c = -1.0;
      corr(i, j) = c;
    }
  }
}

return corr;
}


inline void write_biomass_covariance_diagnostics_csv(
  const std::string &path,
  const std::vector<opakapaka_example::Observation> &data,
  const std::vector<double> &u_hat, const Eigen::SparseMatrix<double> &h_uu)
{
std::ofstream out(path);
out << "metric,value,note\n";

const Eigen::MatrixXd log_b_cov =
    compute_log_b_covariance_submatrix(data, u_hat, h_uu);
const Eigen::MatrixXd biomass_cov = log_cov_to_biomass_cov(log_b_cov, u_hat);
const Eigen::MatrixXd biomass_corr =
    quadra::uncertainty::covariance_to_correlation_matrix(biomass_cov);

const Eigen::Index n = biomass_cov.rows();

bool finite_all = true;
bool positive_diag = true;
double min_diag = std::numeric_limits<double>::infinity();
double max_diag = -std::numeric_limits<double>::infinity();

for (Eigen::Index i = 0; i < n; ++i)
{
  const double v = biomass_cov(i, i);
  if (!std::isfinite(v))
    finite_all = false;
  if (!(v > 0.0))
    positive_diag = false;
  if (std::isfinite(v))
  {
    min_diag = std::min(min_diag, v);
    max_diag = std::max(max_diag, v);
  }

  for (Eigen::Index j = 0; j < n; ++j)
  {
    if (!std::isfinite(biomass_cov(i, j)))
      finite_all = false;
  }
}

double max_abs_asymmetry = 0.0;
if (n > 0)
{
  max_abs_asymmetry =
      (biomass_cov - biomass_cov.transpose()).cwiseAbs().maxCoeff();
}

bool ldlt_success = false;
double min_eigenvalue = std::numeric_limits<double>::quiet_NaN();
double max_eigenvalue = std::numeric_limits<double>::quiet_NaN();

if (n > 0 && finite_all)
{
  Eigen::LDLT<Eigen::MatrixXd> ldlt(biomass_cov);
  ldlt_success = (ldlt.info() == Eigen::Success &&
                  (ldlt.vectorD().array() > -1.0e-10).all());

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(
      0.5 * (biomass_cov + biomass_cov.transpose()));
  if (eig.info() == Eigen::Success)
  {
    min_eigenvalue = eig.eigenvalues().minCoeff();
    max_eigenvalue = eig.eigenvalues().maxCoeff();
  }
}

double mean_nearest_neighbor_corr = std::numeric_limits<double>::quiet_NaN();
double min_nearest_neighbor_corr = std::numeric_limits<double>::quiet_NaN();
double max_nearest_neighbor_corr = std::numeric_limits<double>::quiet_NaN();

if (n > 1)
{
  double sum = 0.0;
  int count = 0;
  min_nearest_neighbor_corr = std::numeric_limits<double>::infinity();
  max_nearest_neighbor_corr = -std::numeric_limits<double>::infinity();

  for (Eigen::Index i = 0; i + 1 < n; ++i)
  {
    const double c = biomass_corr(i, i + 1);
    if (std::isfinite(c))
    {
      sum += c;
      ++count;
      min_nearest_neighbor_corr = std::min(min_nearest_neighbor_corr, c);
      max_nearest_neighbor_corr = std::max(max_nearest_neighbor_corr, c);
    }
  }

  if (count > 0)
  {
    mean_nearest_neighbor_corr = sum / static_cast<double>(count);
  }
}

double mean_lag2_corr = std::numeric_limits<double>::quiet_NaN();
if (n > 2)
{
  double sum = 0.0;
  int count = 0;
  for (Eigen::Index i = 0; i + 2 < n; ++i)
  {
    const double c = biomass_corr(i, i + 2);
    if (std::isfinite(c))
    {
      sum += c;
      ++count;
    }
  }
  if (count > 0)
    mean_lag2_corr = sum / static_cast<double>(count);
}

double mean_lag5_corr = std::numeric_limits<double>::quiet_NaN();
if (n > 5)
{
  double sum = 0.0;
  int count = 0;
  for (Eigen::Index i = 0; i + 5 < n; ++i)
  {
    const double c = biomass_corr(i, i + 5);
    if (std::isfinite(c))
    {
      sum += c;
      ++count;
    }
  }
  if (count > 0)
    mean_lag5_corr = sum / static_cast<double>(count);
}

const bool valid_covariance =
    finite_all && positive_diag && max_abs_asymmetry < 1.0e-8 &&
    ldlt_success && std::isfinite(min_eigenvalue) && min_eigenvalue > -1.0e-8;

auto emit = [&](const std::string &metric, const auto &value,
                const std::string &note)
{
  out << metric << "," << value << "," << note << "\n";
};

emit("n_years", n, "number of fitted biomass states in covariance block");
emit("finite_all", finite_all ? "yes" : "no",
     "all covariance entries finite");
emit("positive_diagonal", positive_diag ? "yes" : "no",
     "all variances positive");
emit("valid_covariance", valid_covariance ? "yes" : "no",
     "finite positive-diagonal symmetric positive-semidefinite check");
emit("ldlt_success", ldlt_success ? "yes" : "no",
     "dense LDLT check on biomass covariance matrix");
emit("max_abs_asymmetry", max_abs_asymmetry,
     "max absolute covariance asymmetry");
emit("min_variance", min_diag, "minimum biomass variance");
emit("max_variance", max_diag, "maximum biomass variance");
emit("min_eigenvalue", min_eigenvalue, "self-adjoint eigenvalue diagnostic");
emit("max_eigenvalue", max_eigenvalue, "self-adjoint eigenvalue diagnostic");
emit("mean_nearest_neighbor_corr", mean_nearest_neighbor_corr,
     "average Corr(B_t,B_tplus1)");
emit("min_nearest_neighbor_corr", min_nearest_neighbor_corr,
     "minimum Corr(B_t,B_tplus1)");
emit("max_nearest_neighbor_corr", max_nearest_neighbor_corr,
     "maximum Corr(B_t,B_tplus1)");
emit("mean_lag2_corr", mean_lag2_corr, "average Corr(B_t,B_tplus2)");
emit("mean_lag5_corr", mean_lag5_corr, "average Corr(B_t,B_tplus5)");
}


inline void write_biomass_covariance_matrix_csv(
  const std::string &path,
  const std::vector<opakapaka_example::Observation> &data,
  const std::vector<double> &u_hat, const Eigen::SparseMatrix<double> &h_uu)
{
std::ofstream out(path);

const std::size_t n = std::min(data.size(), u_hat.size());
if (n == 0)
{
  out << "year\n";
  return;
}

std::vector<int> indices;
indices.reserve(n);
for (std::size_t i = 0; i < n; ++i)
{
  indices.push_back(static_cast<int>(i));
}

const auto log_b_cov =
    quadra::uncertainty::selected_inverse_submatrix_from_spd_hessian(h_uu,
                                                                     indices);

out << "year";
for (std::size_t j = 0; j < n; ++j)
{
  out << ",B_year_" << data[j].year;
}
out << "\n";

for (std::size_t i = 0; i < n; ++i)
{
  out << data[i].year;

  const double b_i = std::exp(u_hat[i]);

  for (std::size_t j = 0; j < n; ++j)
  {
    double cov_biomass = std::numeric_limits<double>::quiet_NaN();

    if (log_b_cov.success &&
        i < static_cast<std::size_t>(log_b_cov.covariance.rows()) &&
        j < static_cast<std::size_t>(log_b_cov.covariance.cols()))
    {
      const double b_j = std::exp(u_hat[j]);
      cov_biomass = b_i * b_j *
                    log_b_cov.covariance(static_cast<Eigen::Index>(i),
                                         static_cast<Eigen::Index>(j));
    }

    out << "," << cov_biomass;
  }

  out << "\n";
}
}


inline void write_biomass_correlation_matrix_csv(
  const std::string &path,
  const std::vector<opakapaka_example::Observation> &data,
  const std::vector<double> &u_hat, const Eigen::SparseMatrix<double> &h_uu)
{
std::ofstream out(path);

const std::size_t n = std::min(data.size(), u_hat.size());
if (n == 0)
{
  out << "year\n";
  return;
}

std::vector<int> indices;
indices.reserve(n);
for (std::size_t i = 0; i < n; ++i)
{
  indices.push_back(static_cast<int>(i));
}

const auto log_b_cov =
    quadra::uncertainty::selected_inverse_submatrix_from_spd_hessian(h_uu,
                                                                     indices);

out << "year";
for (std::size_t j = 0; j < n; ++j)
{
  out << ",B_year_" << data[j].year;
}
out << "\n";

for (std::size_t i = 0; i < n; ++i)
{
  out << data[i].year;

  for (std::size_t j = 0; j < n; ++j)
  {
    double corr = std::numeric_limits<double>::quiet_NaN();

    if (log_b_cov.success &&
        i < static_cast<std::size_t>(log_b_cov.covariance.rows()) &&
        j < static_cast<std::size_t>(log_b_cov.covariance.cols()))
    {
      const double vii = log_b_cov.covariance(static_cast<Eigen::Index>(i),
                                              static_cast<Eigen::Index>(i));
      const double vjj = log_b_cov.covariance(static_cast<Eigen::Index>(j),
                                              static_cast<Eigen::Index>(j));
      const double vij = log_b_cov.covariance(static_cast<Eigen::Index>(i),
                                              static_cast<Eigen::Index>(j));

      if (std::isfinite(vii) && std::isfinite(vjj) && std::isfinite(vij) &&
          vii > 0.0 && vjj > 0.0)
      {
        corr = vij / std::sqrt(vii * vjj);
        if (corr > 1.0 && corr < 1.0 + 1.0e-10)
          corr = 1.0;
        if (corr < -1.0 && corr > -1.0 - 1.0e-10)
          corr = -1.0;
      }
    }

    out << "," << corr;
  }

  out << "\n";
}
}


inline void write_biomass_correlation_decay_csv(
  const std::string &path,
  const std::vector<opakapaka_example::Observation> &data,
  const std::vector<double> &u_hat, const Eigen::SparseMatrix<double> &h_uu)
{
std::ofstream out(path);
out << "lag,count,mean_correlation,min_correlation,max_correlation\n";

const Eigen::MatrixXd log_b_cov =
    compute_log_b_covariance_submatrix(data, u_hat, h_uu);
const Eigen::MatrixXd biomass_cov = log_cov_to_biomass_cov(log_b_cov, u_hat);
const Eigen::MatrixXd biomass_corr =
    quadra::uncertainty::covariance_to_correlation_matrix(biomass_cov);

const Eigen::Index n = biomass_corr.rows();

for (Eigen::Index lag = 0; lag < n; ++lag)
{
  double sum = 0.0;
  double min_corr = std::numeric_limits<double>::infinity();
  double max_corr = -std::numeric_limits<double>::infinity();
  int count = 0;

  for (Eigen::Index i = 0; i + lag < n; ++i)
  {
    const double c = biomass_corr(i, i + lag);
    if (std::isfinite(c))
    {
      sum += c;
      min_corr = std::min(min_corr, c);
      max_corr = std::max(max_corr, c);
      ++count;
    }
  }

  const double mean_corr = count > 0
                               ? sum / static_cast<double>(count)
                               : std::numeric_limits<double>::quiet_NaN();

  out << lag << "," << count << "," << mean_corr << "," << min_corr << ","
      << max_corr << "\n";
}
}


}  // namespace opakapaka_example

// Compatibility aliases for the current Opakapaka driver.
using opakapaka_example::compute_log_b_covariance_submatrix;
using opakapaka_example::log_cov_to_biomass_cov;
using opakapaka_example::covariance_to_correlation;
using opakapaka_example::write_biomass_covariance_diagnostics_csv;
using opakapaka_example::write_biomass_covariance_matrix_csv;
using opakapaka_example::write_biomass_correlation_matrix_csv;
using opakapaka_example::write_biomass_correlation_decay_csv;
