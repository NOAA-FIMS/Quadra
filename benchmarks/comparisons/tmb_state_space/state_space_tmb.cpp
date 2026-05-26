#include <TMB.hpp>

template <class Type> Type objective_function<Type>::operator()() {
  DATA_INTEGER(n_state);

  PARAMETER(mu);
  PARAMETER(log_sigma_obs);
  PARAMETER(log_sigma_rw);
  PARAMETER_VECTOR(x);

  Type sigma_obs = exp(log_sigma_obs);
  Type sigma_rw = exp(log_sigma_rw);

  Type nll = Type(0.0);

  for (int t = 0; t < n_state; ++t) {
    Type y_t = Type(1.0) + Type(0.10) * sin(Type(0.03) * Type(t)) +
               Type(0.05) * cos(Type(0.11) * Type(t));

    Type obs_resid = (y_t - (mu + x(t))) / sigma_obs;

    nll += Type(0.5) * obs_resid * obs_resid + log_sigma_obs;

    if (t == 0) {
      Type z = x(t) / sigma_rw;
      nll += Type(0.5) * z * z + log_sigma_rw;
    } else {
      Type z = (x(t) - x(t - 1)) / sigma_rw;
      nll += Type(0.5) * z * z + log_sigma_rw;
    }
  }

  nll += Type(0.001) * mu * mu;
  nll += Type(0.001) * log_sigma_obs * log_sigma_obs;
  nll += Type(0.001) * log_sigma_rw * log_sigma_rw;

  return nll;
}
