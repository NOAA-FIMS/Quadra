#pragma once

#include "../quadra/bigeye_age_structured.hpp"

namespace pifsc_bigeye_tuna {
namespace level21_layout {

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

} // namespace level21_layout
} // namespace pifsc_bigeye_tuna
