#pragma once

namespace bigeye_v2 {

template <typename T> struct LifeHistoryParameters {
  T log_m_young_offset = T(0.0);
  T log_m_old_offset = T(0.0);
};

} // namespace bigeye_v2
