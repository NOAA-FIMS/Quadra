#include <TMB.hpp>

template <class Type> Type objective_function<Type>::operator()() {
  DATA_VECTOR(index_obs);
  DATA_VECTOR(catch_obs);
  DATA_INTEGER(n_years);

  PARAMETER(log_q);
  PARAMETER_VECTOR(log_B);

  const Type r = Type(0.34);
  const Type K = Type(950.0);
  const Type sigma_process = Type(0.10);
  const Type sigma_index = Type(0.08);
  const Type B0 = Type(0.82) * K;

  Type nll = Type(0.0);

  nll -= dnorm(log_B(0), log(B0), sigma_process, true);

  for (int y = 0; y < n_years; ++y) {
    const Type B = exp(log_B(y));
    const Type pred_index = exp(log_q) * B;

    nll -= dnorm(log(index_obs(y)), log(pred_index), sigma_index, true);

    if (y < n_years - 1) {
      const Type surplus = r * B * (Type(1.0) - B / K);
      const Type raw_next = B + surplus - catch_obs(y);

      // Stable positive floor. Avoids log(<=0) under bad trial points.
      const Type floor = Type(1.0);
      const Type B_next_pred =
          CppAD::CondExpGt(raw_next, floor, raw_next, floor);

      nll -= dnorm(log_B(y + 1), log(B_next_pred), sigma_process, true);
    }
  }

  REPORT(log_B);
  REPORT(log_q);
  ADREPORT(log_q);

  return nll;
}
