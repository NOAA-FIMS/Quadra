#pragma once

#include "../../core/optimizer.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

namespace opakapaka_example {

struct Observation {
  int year = 0;
  double catch_mt = 0.0;
  double index = 0.0;
};

struct ProjectionScenario {
  std::string name;
  double catch_multiplier = 1.0;
};

struct ProjectionRow {
  std::string scenario;
  int year = 0;
  double catch_mt = 0.0;
  double biomass = 0.0;
  double index = 0.0;
};

struct ProjectionOptions {
  int start_year = 2025;
  int years = 10;
  std::vector<ProjectionScenario> scenarios;
};

inline double square(double x) { return x * x; }

inline double safe_log(double x) {
  return std::log(std::max(x, 1.0e-12));
}

inline void add_parameter(quadra::ParameterVector &params,
                          const std::string &name, double value,
                          bool is_random) {
  params.add(quadra::Parameter(name, value,
                               quadra::ParameterTransform::Identity,
                               is_random));
}

inline std::vector<Observation> make_synthetic_opakapaka_data() {
  // Synthetic/public-data-safe data with opakapaka-style scale and trajectory.
  // This is not an official assessment data set.
  std::vector<Observation> data;
  data.reserve(30);

  double biomass = 780.0;
  const double r = 0.34;
  const double K = 950.0;
  const double q = 0.00112;

  for (int i = 0; i < 30; ++i) {
    const int year = 1995 + i;

    // Smooth deterministic catch series. Keep it small enough to avoid
    // pathological toy-model behavior while still forcing a signal.
    const double catch_mt =
        18.0 + 3.0 * std::sin(0.40 * i) + 1.5 * std::cos(0.17 * i);

    // Public-safe synthetic observation noise.
    const double noise = std::exp(0.055 * std::sin(0.73 * i) -
                                  0.035 * std::cos(0.31 * i));
    const double index = q * biomass * noise;

    data.push_back({year, catch_mt, index});

    biomass = biomass + r * biomass * (1.0 - biomass / K) - catch_mt;
    biomass = std::max(20.0, biomass);
  }

  return data;
}

class OpakapakaProjectionModel {
public:
  explicit OpakapakaProjectionModel(std::vector<Observation> observations)
      : data_(std::move(observations)) {
    if (data_.empty()) {
      throw std::invalid_argument("OpakapakaProjectionModel requires data");
    }
  }

  quadra::ParameterVector initial_parameters() const {
    quadra::ParameterVector params;

    // Fixed effects. Keep key biological quantities fixed in this first clean
    // public-API example so the fit is stable and easy to read.
    //
    // log_q is estimated to demonstrate optimize_lbfgs without adding a large
    // identifiability problem to the synthetic example.
    add_parameter(params, "log_q", std::log(0.0010), false);

    // Random effects: latent log-biomass by year.
    for (std::size_t i = 0; i < data_.size(); ++i) {
      const double frac = 0.82 - 0.0015 * static_cast<double>(i);
      add_parameter(params, "log_B_" + std::to_string(i),
                    std::log(950.0 * frac), true);
    }

    return params;
  }

  template <class T> T operator()(const std::vector<T> &par) const {
    const int n = static_cast<int>(data_.size());

    const T log_q = par[0];
    const T q = exp(log_q);

    // Fixed biological parameters for readable example.
    const T r = T(0.34);
    const T K = T(950.0);
    const T sigma_process = T(0.10);
    const T sigma_index = T(0.08);

    T nll = T(0.0);

    // Biomass state prior.
    const T log_B0_expected = log(T(0.82) * K);
    const T log_B0 = par[1];
    nll += T(0.5) * square((log_B0 - log_B0_expected) / T(0.15)) +
           log(T(0.15));

    for (int t = 0; t < n; ++t) {
      const T log_Bt = par[1 + t];
      const T Bt = exp(log_Bt);

      // Index likelihood.
      const T pred_index = q * Bt;
      const T obs_index = T(data_[static_cast<std::size_t>(t)].index);
      nll += T(0.5) *
                 square((log(obs_index) - log(pred_index)) / sigma_index) +
             log(sigma_index);

      // Process equation for next biomass.
      if (t + 1 < n) {
        const T catch_t = T(data_[static_cast<std::size_t>(t)].catch_mt);
        T B_next_pred = Bt + r * Bt * (T(1.0) - Bt / K) - catch_t;
        // Smooth positive guard for the toy projection model. This avoids
        // branching/comparison on AD scalar types in the example code.
        const T guarded_B_next_pred =
            sqrt(B_next_pred * B_next_pred + T(1.0e-8));

        const T log_B_next = par[1 + t + 1];
        nll += T(0.5) *
                   square((log_B_next - log(guarded_B_next_pred)) /
                          sigma_process) +
               log(sigma_process);
      }
    }

    return nll;
  }

  std::vector<ProjectionRow>
  project(const quadra::OptResult &fit,
          const ProjectionOptions &options) const {
    if (fit.u_hat.empty()) {
      throw std::runtime_error("Projection requires fit.u_hat");
    }

    const double log_q = fit.par.at(0);
    const double q = std::exp(log_q);

    const double r = 0.34;
    const double K = 950.0;

    const double terminal_biomass = std::exp(fit.u_hat.back());
    const double recent_catch = data_.back().catch_mt;

    std::vector<ProjectionScenario> scenarios = options.scenarios;
    if (scenarios.empty()) {
      scenarios = {
          {"zero_catch", 0.0},
          {"status_quo", 1.0},
          {"low_catch", 0.75},
          {"high_catch", 1.25},
      };
    }

    std::vector<ProjectionRow> rows;

    for (const auto &scenario : scenarios) {
      double biomass = terminal_biomass;

      for (int y = 0; y < options.years; ++y) {
        const int year = options.start_year + y;
        const double catch_mt = recent_catch * scenario.catch_multiplier;

        biomass = biomass + r * biomass * (1.0 - biomass / K) - catch_mt;
        biomass = std::max(1.0, biomass);

        rows.push_back(
            {scenario.name, year, catch_mt, biomass, q * biomass});
      }
    }

    return rows;
  }

  const std::vector<Observation> &data() const { return data_; }

private:
  std::vector<Observation> data_;
};

inline void write_fit_summary_csv(const std::string &path,
                                  const quadra::OptResult &fit) {
  std::ofstream out(path);
  out << "field,value\n";
  out << "objective," << fit.value << "\n";
  out << "grad_norm," << fit.grad_norm << "\n";
  out << "iterations," << fit.iterations << "\n";
  out << "converged," << (fit.converged ? "yes" : "no") << "\n";
  out << "message," << fit.message << "\n";
  out << "random_effects," << fit.pattern.random_effect_count << "\n";
  out << "detected_structure," << fit.pattern.detected_structure << "\n";
  out << "backend," << fit.pattern.backend << "\n";
  out << "solver," << fit.pattern.solver << "\n";
  out << "complexity," << fit.pattern.complexity << "\n";
  out << "bandwidth," << fit.pattern.bandwidth << "\n";
  out << "hessian_nonzeros," << fit.pattern.nonzeros << "\n";
}

inline void write_projection_csv(const std::string &path,
                                 const std::vector<ProjectionRow> &rows) {
  std::ofstream out(path);
  out << "scenario,year,catch_mt,biomass,index\n";
  out << std::fixed << std::setprecision(6);
  for (const auto &row : rows) {
    out << row.scenario << "," << row.year << "," << row.catch_mt << ","
        << row.biomass << "," << row.index << "\n";
  }
}

} // namespace opakapaka_example
