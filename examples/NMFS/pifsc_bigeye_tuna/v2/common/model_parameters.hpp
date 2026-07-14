#pragma once

namespace bigeye_v2 {

template <typename T> struct BigeyeModelParameters {
  T log_m_young_offset = T(0.0);
  T log_m_old_offset = T(0.0);
  T r0 = T(1000.0);
  T sel_a50 = T(5.0);
  T sel_slope = T(1.0);
  T fbar = T(0.2);
  T catch_sigma = T(0.1);
  T q_index = T(0.01);
  T index_sigma = T(0.2);
};

} // namespace bigeye_v2
