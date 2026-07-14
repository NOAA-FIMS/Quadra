#pragma once

#include "../../../../../core/optimizer.hpp"
#include "../objective/bigeye_quadra_objective.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace pifsc_bigeye_tuna {

struct RecruitmentDiagnosticRow {
  std::size_t index = 0;
  int year = 0;
  double rec_dev = std::numeric_limits<double>::quiet_NaN();
  double recruitment_multiplier = std::numeric_limits<double>::quiet_NaN();
  double prior_z = std::numeric_limits<double>::quiet_NaN();
  double prior_nll = std::numeric_limits<double>::quiet_NaN();
  double abs_rec_dev = std::numeric_limits<double>::quiet_NaN();
};

struct RecruitmentDiagnostics {
  double sigma_rec_dev = 0.35;
  std::size_t n = 0;
  double mean = std::numeric_limits<double>::quiet_NaN();
  double sd = std::numeric_limits<double>::quiet_NaN();
  double max_abs = std::numeric_limits<double>::quiet_NaN();
  int max_abs_year = 0;
  double lag1_correlation = std::numeric_limits<double>::quiet_NaN();
  double roughness = std::numeric_limits<double>::quiet_NaN();
  double total_prior_nll = 0.0;
  std::string interpretation;
  std::vector<RecruitmentDiagnosticRow> rows;
};

inline double sample_sd(const std::vector<double> &x) {
  if (x.size() < 2)
    return std::numeric_limits<double>::quiet_NaN();

  const double mean =
      std::accumulate(x.begin(), x.end(), 0.0) / static_cast<double>(x.size());

  double ss = 0.0;
  for (const double v : x)
    ss += (v - mean) * (v - mean);

  return std::sqrt(ss / static_cast<double>(x.size() - 1));
}

inline double lag1_corr(const std::vector<double> &x) {
  if (x.size() < 3)
    return std::numeric_limits<double>::quiet_NaN();

  std::vector<double> a;
  std::vector<double> b;
  a.reserve(x.size() - 1);
  b.reserve(x.size() - 1);

  for (std::size_t i = 1; i < x.size(); ++i) {
    a.push_back(x[i - 1]);
    b.push_back(x[i]);
  }

  const double ma =
      std::accumulate(a.begin(), a.end(), 0.0) / static_cast<double>(a.size());
  const double mb =
      std::accumulate(b.begin(), b.end(), 0.0) / static_cast<double>(b.size());

  double num = 0.0;
  double va = 0.0;
  double vb = 0.0;

  for (std::size_t i = 0; i < a.size(); ++i) {
    const double da = a[i] - ma;
    const double db = b[i] - mb;
    num += da * db;
    va += da * da;
    vb += db * db;
  }

  const double den = std::sqrt(va * vb);
  return den > 0.0 ? num / den : std::numeric_limits<double>::quiet_NaN();
}

inline double first_difference_roughness(const std::vector<double> &x) {
  if (x.size() < 2)
    return std::numeric_limits<double>::quiet_NaN();

  double ss = 0.0;
  for (std::size_t i = 1; i < x.size(); ++i) {
    const double d = x[i] - x[i - 1];
    ss += d * d;
  }

  return std::sqrt(ss / static_cast<double>(x.size() - 1));
}

inline RecruitmentDiagnostics
make_recruitment_diagnostics(const BigeyeQuadraObjective &objective,
                             const quadra::OptResult &fit,
                             double sigma_rec_dev = 0.35) {
  RecruitmentDiagnostics out;
  out.sigma_rec_dev = sigma_rec_dev;

  const auto years = objective.unique_years();
  const std::size_t n = std::min(years.size(), fit.u_hat.size());
  out.n = n;

  std::vector<double> devs;
  devs.reserve(n);

  for (std::size_t i = 0; i < n; ++i) {
    const double rec_dev = fit.u_hat[i];

    RecruitmentDiagnosticRow row;
    row.index = i;
    row.year = years[i];
    row.rec_dev = rec_dev;
    row.recruitment_multiplier = std::exp(rec_dev);
    row.prior_z = rec_dev / sigma_rec_dev;
    row.prior_nll = 0.5 * row.prior_z * row.prior_z;
    row.abs_rec_dev = std::abs(rec_dev);

    out.total_prior_nll += row.prior_nll;
    out.rows.push_back(row);
    devs.push_back(rec_dev);
  }

  if (!devs.empty()) {
    out.mean = std::accumulate(devs.begin(), devs.end(), 0.0) /
               static_cast<double>(devs.size());
    out.sd = sample_sd(devs);
    out.lag1_correlation = lag1_corr(devs);
    out.roughness = first_difference_roughness(devs);

    auto max_it = std::max_element(out.rows.begin(), out.rows.end(),
                                   [](const RecruitmentDiagnosticRow &a,
                                      const RecruitmentDiagnosticRow &b) {
                                     return a.abs_rec_dev < b.abs_rec_dev;
                                   });

    if (max_it != out.rows.end()) {
      out.max_abs = max_it->abs_rec_dev;
      out.max_abs_year = max_it->year;
    }
  }

  if (std::isfinite(out.sd) && out.sd > 0.75 * sigma_rec_dev) {
    out.interpretation = "Recruitment deviations are using a substantial "
                         "fraction of the prior scale.";
  } else if (std::isfinite(out.sd)) {
    out.interpretation =
        "Recruitment deviations are modest relative to the prior scale.";
  } else {
    out.interpretation = "Insufficient recruitment deviations for summary.";
  }

  return out;
}

inline void
write_recruitment_diagnostics_csv(const RecruitmentDiagnostics &diagnostics,
                                  const std::string &path) {
  std::ofstream out(path);
  if (!out)
    throw std::runtime_error("Could not open recruitment diagnostics CSV: " +
                             path);

  out << "section,metric,target,value,extra\n";
  out << std::setprecision(15);

  out << "summary,sigma_rec_dev,," << diagnostics.sigma_rec_dev << ",\n";
  out << "summary,n,," << diagnostics.n << ",\n";
  out << "summary,mean,," << diagnostics.mean << ",\n";
  out << "summary,sd,," << diagnostics.sd << ",\n";
  out << "summary,max_abs,," << diagnostics.max_abs
      << ",year=" << diagnostics.max_abs_year << "\n";
  out << "summary,lag1_correlation,," << diagnostics.lag1_correlation << ",\n";
  out << "summary,roughness,," << diagnostics.roughness << ",\n";
  out << "summary,total_prior_nll,," << diagnostics.total_prior_nll << ",\n";
  out << "summary,interpretation,,\"" << diagnostics.interpretation << "\",\n";

  for (const auto &row : diagnostics.rows) {
    out << "recruitment,rec_dev," << row.year << "," << row.rec_dev
        << ",index=" << row.index << "\n";
    out << "recruitment,multiplier," << row.year << ","
        << row.recruitment_multiplier << ",index=" << row.index << "\n";
    out << "recruitment,prior_z," << row.year << "," << row.prior_z
        << ",index=" << row.index << "\n";
    out << "recruitment,prior_nll," << row.year << "," << row.prior_nll
        << ",index=" << row.index << "\n";
  }
}

inline void
write_recruitment_diagnostics_text(const RecruitmentDiagnostics &diagnostics,
                                   const std::string &path) {
  std::ofstream out(path);
  if (!out)
    throw std::runtime_error("Could not open recruitment diagnostics text: " +
                             path);

  out << "Recruitment Diagnostics\n";
  out << "=======================\n\n";
  out << std::setprecision(15);

  out << "Summary\n";
  out << "-------\n";
  out << "sigma_rec_dev:          " << diagnostics.sigma_rec_dev << "\n";
  out << "n:                      " << diagnostics.n << "\n";
  out << "mean_rec_dev:           " << diagnostics.mean << "\n";
  out << "sd_rec_dev:             " << diagnostics.sd << "\n";
  out << "max_abs_rec_dev:        " << diagnostics.max_abs << "\n";
  out << "max_abs_year:           " << diagnostics.max_abs_year << "\n";
  out << "lag1_correlation:       " << diagnostics.lag1_correlation << "\n";
  out << "roughness:              " << diagnostics.roughness << "\n";
  out << "total_rec_prior_nll:    " << diagnostics.total_prior_nll << "\n";
  out << "interpretation:         " << diagnostics.interpretation << "\n\n";

  out << "Interpretation notes\n";
  out << "--------------------\n";
  out << "Recruitment deviations should be checked for magnitude, persistence, "
         "and\n";
  out << "roughness. Large persistent runs can indicate that recruitment is "
         "absorbing\n";
  out << "misspecified fleet, selectivity, index, or catch processes rather "
         "than\n";
  out << "representing biological recruitment variation.\n\n";

  out << "Rows\n";
  out << "----\n";
  out << "index,year,rec_dev,recruitment_multiplier,prior_z,prior_nll\n";

  for (const auto &row : diagnostics.rows) {
    out << row.index << "," << row.year << "," << row.rec_dev << ","
        << row.recruitment_multiplier << "," << row.prior_z << ","
        << row.prior_nll << "\n";
  }
}

inline void write_recruitment_diagnostics(
    const std::string &text_path, const std::string &csv_path,
    const BigeyeQuadraObjective &objective, const quadra::OptResult &fit,
    double sigma_rec_dev = 0.35) {
  const auto diagnostics =
      make_recruitment_diagnostics(objective, fit, sigma_rec_dev);

  write_recruitment_diagnostics_text(diagnostics, text_path);
  write_recruitment_diagnostics_csv(diagnostics, csv_path);
}

} // namespace pifsc_bigeye_tuna

using pifsc_bigeye_tuna::RecruitmentDiagnostics;
using pifsc_bigeye_tuna::write_recruitment_diagnostics;
