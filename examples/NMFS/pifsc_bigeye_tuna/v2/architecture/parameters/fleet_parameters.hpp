#pragma once

namespace bigeye_v2 {

template <typename T>
struct FleetParameters {
  T sel_a50 = T(5.0);
  T sel_slope = T(1.0);
  T fbar = T(0.2);
  T q_index = T(0.01);
  T catch_sigma = T(0.1);
  T index_sigma = T(0.2);
};

} // namespace bigeye_v2
