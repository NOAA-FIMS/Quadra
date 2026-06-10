#pragma once

#include <array>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

namespace sefsc_red_snapper {

struct Observation {
  int year = 0;
  double catch_mt = 0.0;
  double index = 0.0;
  std::array<double, 10> age_comp{};
};

struct ProjectionScenario {
  std::string scenario;
  int projection_year = 0;
  double catch_mt = 0.0;
};

struct DerivedRow {
  int year = 0;
  double biomass = 0.0;
  double ssb_proxy = 0.0;
  double depletion = 0.0;
  double f_proxy = 0.0;
  double index_hat = 0.0;
};

// Level-0 placeholder model:
// This is intentionally minimal. The next patch should replace this with
// Quadra AD/Laplace evaluation and recruitment deviations as random effects.
class RedSnapperModel {
 public:
  explicit RedSnapperModel(std::vector<Observation> obs)
      : observations_(std::move(obs)) {}

  const std::vector<Observation>& observations() const { return observations_; }

  std::vector<DerivedRow> deterministic_trajectory(double log_r0,
                                                    double log_q,
                                                    double log_f) const {
    const double r0 = std::exp(log_r0);
    const double q = std::exp(log_q);
    const double f = std::exp(log_f);

    std::vector<DerivedRow> out;
    out.reserve(observations_.size());

    double biomass = r0;
    const double unfished = r0;

    for (const auto& obs : observations_) {
      biomass = std::max(1.0, biomass + 0.25 * r0 - obs.catch_mt - 0.05 * biomass);
      DerivedRow row;
      row.year = obs.year;
      row.biomass = biomass;
      row.ssb_proxy = 0.35 * biomass;
      row.depletion = biomass / unfished;
      row.f_proxy = f * obs.catch_mt / std::max(1.0, biomass);
      row.index_hat = q * biomass;
      out.push_back(row);
    }

    return out;
  }

 private:
  std::vector<Observation> observations_;
};

}  // namespace sefsc_red_snapper
