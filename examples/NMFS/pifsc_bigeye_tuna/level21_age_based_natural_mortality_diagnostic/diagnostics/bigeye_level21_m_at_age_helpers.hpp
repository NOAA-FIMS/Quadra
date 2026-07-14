#pragma once

#include "../objective/bigeye_quadra_objective.hpp"
#include <array>
#include <cmath>
#include <vector>

namespace pifsc_bigeye_tuna {
namespace level21_m_helpers {

inline constexpr int kBaseFixed = 3;
inline constexpr int kMParamOffset = kBaseFixed;
inline constexpr int kMParams = 2;
inline constexpr int kLonglineSelOffset = kMParamOffset + kMParams;
inline constexpr int kLonglineSelDevs = kAges;
inline constexpr int kInitialDevOffset = kLonglineSelOffset + kLonglineSelDevs;
inline constexpr int kInitialDevs = kAges;
inline constexpr int kPurseSeineSelOffset = kInitialDevOffset + kInitialDevs;
inline constexpr int kPurseSeineSelDevs = kAges;
inline constexpr int kRecruitmentOffset =
    kPurseSeineSelOffset + kPurseSeineSelDevs;

inline std::array<double, kAges>
m_at_age_from_level21_par(const std::vector<double> &par) {
  const double adult_m = 0.45;
  double log_m_young_offset = 0.0;
  double log_m_old_offset = 0.0;

  if (par.size() > static_cast<std::size_t>(kMParamOffset + 0)) {
    log_m_young_offset = par[static_cast<std::size_t>(kMParamOffset + 0)];
  }
  if (par.size() > static_cast<std::size_t>(kMParamOffset + 1)) {
    log_m_old_offset = par[static_cast<std::size_t>(kMParamOffset + 1)];
  }

  const double m_young = adult_m * std::exp(log_m_young_offset);
  const double m_old = adult_m * std::exp(log_m_old_offset);

  std::array<double, kAges> m_at_age{};
  for (int a = 0; a < kAges; ++a) {
    if (a < 2) {
      m_at_age[static_cast<std::size_t>(a)] = m_young;
    } else if (a >= 7) {
      m_at_age[static_cast<std::size_t>(a)] = m_old;
    } else {
      m_at_age[static_cast<std::size_t>(a)] = adult_m;
    }
  }
  return m_at_age;
}

inline double adult_m() { return 0.45; }

} // namespace level21_m_helpers
} // namespace pifsc_bigeye_tuna
