#pragma once

#include "../../common/bigeye_constants.hpp"

#include <array>

namespace bigeye_v2 {

template <typename T>
struct LifeHistoryState {
  std::array<T, kAges> m_at_age{};
  std::array<T, kAges> weight_at_age{};
  std::array<T, kAges> maturity_at_age{};
};

} // namespace bigeye_v2
