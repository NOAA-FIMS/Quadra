#pragma once

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace pifsc_bigeye_tuna {

constexpr int kAges = 10;

struct Observation {
  int year = 0;
  double catch_mt = 0.0;
  double index = 0.0;
  std::array<double, kAges> age_comp{};
};

struct AgeStructuredParams {
  double log_r0 = std::log(1200.0);
  double log_m = std::log(0.18);
  double log_fbar = std::log(0.025);
  double log_q = std::log(0.00005);
  double sel_a50 = 4.0;
  double sel_slope = 1.2;
};

struct AgeStructuredRow {
  int year = 0;
  double recruitment = 0.0;
  double total_biomass = 0.0;
  double ssb_proxy = 0.0;
  double depletion = 0.0;
  double fbar = 0.0;
  double catch_obs = 0.0;
  double catch_hat = 0.0;
  double index_obs = 0.0;
  double index_hat = 0.0;
};

double logistic_selectivity(double age, double a50, double slope) {
  return 1.0 / (1.0 + std::exp(-slope * (age - a50)));
}

std::array<double, kAges> default_weight_at_age() {
  return {2.0, 8.0, 18.0, 32.0, 48.0, 65.0, 82.0, 98.0, 112.0, 125.0};
}

std::array<double, kAges> default_maturity_at_age() {
  return {0.00, 0.10, 0.35, 0.65, 0.85, 0.95, 1.00, 1.00, 1.00, 1.00};
}

std::vector<std::string> split_csv_line(const std::string &line) {
  std::vector<std::string> out;
  std::stringstream ss(line);
  std::string item;
  while (std::getline(ss, item, ',')) {
    out.push_back(item);
  }
  return out;
}

std::vector<Observation> read_observations(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("Could not open observations CSV: " + path);
  }

  std::string line;
  std::getline(in, line);

  std::vector<Observation> out;
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }

    const auto fields = split_csv_line(line);
    if (fields.size() != 13) {
      throw std::runtime_error("Expected 13 columns in observations CSV");
    }

    Observation obs;
    obs.year = std::stoi(fields[0]);
    obs.catch_mt = std::stod(fields[1]);
    obs.index = std::stod(fields[2]);
    for (int a = 0; a < kAges; ++a) {
      obs.age_comp[static_cast<std::size_t>(a)] = std::stod(fields[3 + a]);
    }

    double age_comp_sum = 0.0;
    for (double v : obs.age_comp) {
      age_comp_sum += v;
    }
    if (age_comp_sum > 0.0) {
      for (double &v : obs.age_comp) {
        v /= age_comp_sum;
      }
    }
    out.push_back(obs);
  }

  return out;
}

double biomass_from_numbers(const std::array<double, kAges> &n,
                            const std::array<double, kAges> &weight) {
  double out = 0.0;
  for (int a = 0; a < kAges; ++a) {
    out += n[static_cast<std::size_t>(a)] * weight[static_cast<std::size_t>(a)];
  }
  return out;
}

double ssb_from_numbers(const std::array<double, kAges> &n,
                        const std::array<double, kAges> &weight,
                        const std::array<double, kAges> &maturity) {
  double out = 0.0;
  for (int a = 0; a < kAges; ++a) {
    out += n[static_cast<std::size_t>(a)] *
           weight[static_cast<std::size_t>(a)] *
           maturity[static_cast<std::size_t>(a)];
  }
  return out;
}

std::array<double, kAges> unfished_equilibrium_numbers(double r0, double m) {
  std::array<double, kAges> n{};
  n[0] = r0;
  for (int a = 1; a < kAges; ++a) {
    n[static_cast<std::size_t>(a)] =
        n[static_cast<std::size_t>(a - 1)] * std::exp(-m);
  }

  // Plus group.
  n[static_cast<std::size_t>(kAges - 1)] /=
      std::max(1.0e-12, 1.0 - std::exp(-m));

  return n;
}

std::vector<AgeStructuredRow> run_deterministic_age_structured_model(
    const std::vector<Observation> &observations,
    const AgeStructuredParams &params) {
  const auto weight = default_weight_at_age();
  const auto maturity = default_maturity_at_age();

  const double r0 = std::exp(params.log_r0);
  const double m = std::exp(params.log_m);
  const double fbar = std::exp(params.log_fbar);
  const double q = std::exp(params.log_q);

  std::array<double, kAges> selectivity{};
  for (int a = 0; a < kAges; ++a) {
    selectivity[static_cast<std::size_t>(a)] = logistic_selectivity(
        static_cast<double>(a + 1), params.sel_a50, params.sel_slope);
  }

  std::array<double, kAges> n = unfished_equilibrium_numbers(r0, m);
  const double unfished_ssb = ssb_from_numbers(n, weight, maturity);

  std::vector<AgeStructuredRow> rows;
  rows.reserve(observations.size());

  for (const auto &obs : observations) {
    const double biomass = biomass_from_numbers(n, weight);
    const double ssb = ssb_from_numbers(n, weight, maturity);

    double catch_hat = 0.0;
    for (int a = 0; a < kAges; ++a) {
      const auto i = static_cast<std::size_t>(a);
      const double f_a = fbar * selectivity[i];
      const double z_a = m + f_a;
      const double harvest_rate =
          z_a > 0.0 ? (f_a / z_a) * (1.0 - std::exp(-z_a)) : 0.0;
      catch_hat += n[i] * weight[i] * harvest_rate;
    }

    AgeStructuredRow row;
    row.year = obs.year;
    row.recruitment = r0;
    row.total_biomass = biomass;
    row.ssb_proxy = ssb;
    row.depletion = ssb / std::max(1.0e-12, unfished_ssb);
    row.fbar = fbar;
    row.catch_obs = obs.catch_mt;
    row.catch_hat = catch_hat;
    row.index_obs = obs.index;
    row.index_hat = q * biomass;
    rows.push_back(row);

    std::array<double, kAges> next{};
    next[0] = r0;

    for (int a = 1; a < kAges; ++a) {
      const auto prev = static_cast<std::size_t>(a - 1);
      const double f_prev = fbar * selectivity[prev];
      const double z_prev = m + f_prev;
      next[static_cast<std::size_t>(a)] = n[prev] * std::exp(-z_prev);
    }

    // Plus group survivor contribution.
    {
      const auto last = static_cast<std::size_t>(kAges - 1);
      const double f_last = fbar * selectivity[last];
      const double z_last = m + f_last;
      next[last] += n[last] * std::exp(-z_last);
    }

    n = next;
  }

  return rows;
}

void write_age_structured_rows(const std::string &path,
                               const std::vector<AgeStructuredRow> &rows) {
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Could not open output CSV: " + path);
  }

  out << "year,recruitment,total_biomass,ssb_proxy,depletion,Fbar,"
      << "catch_obs,catch_hat,index_obs,index_hat\n";

  out << std::fixed << std::setprecision(6);
  for (const auto &row : rows) {
    out << row.year << "," << row.recruitment << "," << row.total_biomass << ","
        << row.ssb_proxy << "," << row.depletion << "," << row.fbar << ","
        << row.catch_obs << "," << row.catch_hat << "," << row.index_obs << ","
        << row.index_hat << "\n";
  }
}

} // namespace pifsc_bigeye_tuna
