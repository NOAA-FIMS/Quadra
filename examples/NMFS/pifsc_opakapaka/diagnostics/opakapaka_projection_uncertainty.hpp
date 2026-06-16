#pragma once

#include "../quadra/opakapaka_model.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <string>
#include <vector>

namespace opakapaka_example {

struct ProjectionEnvelopeRow
{
std::string scenario;
int year = 0;
std::string quantity;
double estimate = std::numeric_limits<double>::quiet_NaN();
double mean = std::numeric_limits<double>::quiet_NaN();
double median = std::numeric_limits<double>::quiet_NaN();
double lwr_95 = std::numeric_limits<double>::quiet_NaN();
double upr_95 = std::numeric_limits<double>::quiet_NaN();
double se = std::numeric_limits<double>::quiet_NaN();
std::string note;
};


inline double opakapaka_quantile_sorted(const std::vector<double> &sorted,
                                      double p)
{
if (sorted.empty())
  return std::numeric_limits<double>::quiet_NaN();
if (sorted.size() == 1)
  return sorted.front();

const double x = p * static_cast<double>(sorted.size() - 1);
const std::size_t lo = static_cast<std::size_t>(std::floor(x));
const std::size_t hi = std::min<std::size_t>(lo + 1, sorted.size() - 1);
const double w = x - static_cast<double>(lo);
return (1.0 - w) * sorted[lo] + w * sorted[hi];
}


inline ProjectionEnvelopeRow summarize_projection_samples(
  const std::string &scenario, int year, const std::string &quantity,
  double estimate, std::vector<double> samples, const std::string &note)
{
ProjectionEnvelopeRow row;
row.scenario = scenario;
row.year = year;
row.quantity = quantity;
row.estimate = estimate;
row.note = note;

samples.erase(std::remove_if(samples.begin(), samples.end(),
                             [](double x)
                             { return !std::isfinite(x); }),
              samples.end());

if (samples.empty())
{
  return row;
}

const double sum = std::accumulate(samples.begin(), samples.end(), 0.0);
row.mean = sum / static_cast<double>(samples.size());

double ss = 0.0;
if (samples.size() > 1)
{
  for (double x : samples)
  {
    const double d = x - row.mean;
    ss += d * d;
  }
  row.se = std::sqrt(ss / static_cast<double>(samples.size() - 1));
}
else
{
  row.se = 0.0;
}

std::sort(samples.begin(), samples.end());
row.median = opakapaka_quantile_sorted(samples, 0.50);
row.lwr_95 = opakapaka_quantile_sorted(samples, 0.025);
row.upr_95 = opakapaka_quantile_sorted(samples, 0.975);

return row;
}


inline void write_projection_uncertainty_envelopes_csv(
  const std::string &path,
  const std::vector<opakapaka_example::ProjectionRow>
      &deterministic_projection,
  const std::vector<double> &fitted_log_b, double q_hat,
  double terminal_log_b_variance, int n_samples = 1000,
  unsigned seed = 8675309u)
{
std::ofstream out(path);
out << "scenario,year,quantity,estimate,mean,median,lwr_95,upr_95,se,n_"
       "samples,note\n";

if (deterministic_projection.empty() || fitted_log_b.empty() ||
    !std::isfinite(terminal_log_b_variance) ||
    terminal_log_b_variance < 0.0 || n_samples <= 1)
{
  for (const auto &r : deterministic_projection)
  {
    out << r.scenario << "," << r.year << ",biomass," << r.biomass << ",,,,,,"
        << n_samples
        << ",projection_envelope_unavailable_invalid_terminal_variance\n";
    out << r.scenario << "," << r.year << ",index," << r.index << ",,,,,,"
        << n_samples
        << ",projection_envelope_unavailable_invalid_terminal_variance\n";
  }
  return;
}

const double terminal_log_b_hat = fitted_log_b.back();
const double terminal_sd = std::sqrt(terminal_log_b_variance);

// Infer projection dynamics from deterministic rows. This keeps the envelope
// writer independent of assessment-specific model internals:
//   B_{t+1} = B_t + deterministic_increment_t
// where deterministic_increment_t is read from the point projection.
std::map<std::string, std::vector<opakapaka_example::ProjectionRow>>
    by_scenario;
for (const auto &r : deterministic_projection)
{
  by_scenario[r.scenario].push_back(r);
}

std::mt19937 rng(seed);
std::normal_distribution<double> zdist(0.0, 1.0);

for (auto &kv : by_scenario)
{
  auto &rows = kv.second;
  std::sort(rows.begin(), rows.end(),
            [](const auto &a, const auto &b)
            { return a.year < b.year; });

  std::vector<std::vector<double>> biomass_samples(rows.size());
  std::vector<std::vector<double>> index_samples(rows.size());

  for (int s = 0; s < n_samples; ++s)
  {
    double sampled_b =
        std::exp(terminal_log_b_hat + terminal_sd * zdist(rng));

    for (std::size_t t = 0; t < rows.size(); ++t)
    {
      const double previous_point_b =
          (t == 0) ? std::exp(terminal_log_b_hat) : rows[t - 1].biomass;
      const double deterministic_increment =
          rows[t].biomass - previous_point_b;

      sampled_b = std::max(1.0e-12, sampled_b + deterministic_increment);
      const double sampled_index = q_hat * sampled_b;

      biomass_samples[t].push_back(sampled_b);
      index_samples[t].push_back(sampled_index);
    }
  }

  for (std::size_t t = 0; t < rows.size(); ++t)
  {
    auto b_row = summarize_projection_samples(
        rows[t].scenario, rows[t].year, "biomass", rows[t].biomass,
        biomass_samples[t],
        "terminal_state_parametric_envelope_selected_inverse_delta");
    auto i_row = summarize_projection_samples(
        rows[t].scenario, rows[t].year, "index", rows[t].index,
        index_samples[t],
        "terminal_state_parametric_envelope_selected_inverse_delta");

    auto emit = [&](const ProjectionEnvelopeRow &r)
    {
      out << r.scenario << "," << r.year << "," << r.quantity << ","
          << r.estimate << "," << r.mean << "," << r.median << "," << r.lwr_95
          << "," << r.upr_95 << "," << r.se << "," << n_samples << ","
          << r.note << "\n";
    };

    emit(b_row);
    emit(i_row);
  }
}
}


}  // namespace opakapaka_example

// Compatibility aliases for the current Opakapaka driver.
using opakapaka_example::ProjectionEnvelopeRow;
using opakapaka_example::opakapaka_quantile_sorted;
using opakapaka_example::summarize_projection_samples;
using opakapaka_example::write_projection_uncertainty_envelopes_csv;
