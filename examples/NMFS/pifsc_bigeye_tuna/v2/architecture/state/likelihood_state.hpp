#pragma once

namespace bigeye_v2 {

template <typename T>
struct LikelihoodState {
  T catch_nll = T(0.0);
  T index_nll = T(0.0);
  T agecomp_nll = T(0.0);
  T total_nll = T(0.0);
};

} // namespace bigeye_v2
