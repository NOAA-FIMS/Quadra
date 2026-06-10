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
  const Type sigma_initial = Type(0.15);
  const Type B0 = Type(0.82) * K;

  Type nll = Type(0.0);

  nll += Type(0.5) * ((log_B(0) - log(B0)) / sigma_initial) *
         ((log_B(0) - log(B0)) / sigma_initial);

  for (int y = 0; y < n_years; ++y) {
    const Type B = exp(log_B(y));
    const Type pred_index = exp(log_q) * B;

    nll += Type(0.5) * ((log(index_obs(y)) - log(pred_index)) / sigma_index) *
           ((log(index_obs(y)) - log(pred_index)) / sigma_index);

    if (y < n_years - 1) {
      const Type surplus = r * B * (Type(1.0) - B / K);
      const Type raw_next = B + surplus - catch_obs(y);

      // Match Quadra's smooth positive guard for apples-to-apples comparison.
      const Type B_next_pred = sqrt(raw_next * raw_next + Type(1.0e-8));

      nll += Type(0.5) * ((log_B(y + 1) - log(B_next_pred)) / sigma_process) *
             ((log_B(y + 1) - log(B_next_pred)) / sigma_process);
    }
  }

  REPORT(log_B);
  REPORT(log_q);
  ADREPORT(log_q);

  return nll;
}
