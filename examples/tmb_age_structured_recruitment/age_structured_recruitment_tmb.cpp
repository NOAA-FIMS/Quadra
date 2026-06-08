#include <TMB.hpp>

template <class Type> Type logistic(Type x) {
  return Type(1.0) / (Type(1.0) + exp(-x));
}

template <class Type> Type objective_function<Type>::operator()() {
  DATA_VECTOR(index_obs);
  DATA_INTEGER(n_ages);

  PARAMETER(log_R0);
  PARAMETER(log_M);
  PARAMETER(log_q);
  PARAMETER(log_sigma_R);
  PARAMETER(log_sigma_index);
  PARAMETER_VECTOR(x);

  int n_years = index_obs.size();

  Type R0 = exp(log_R0);
  Type M = exp(log_M);
  Type q = exp(log_q);
  Type sigma_R = exp(log_sigma_R);
  Type sigma_index = exp(log_sigma_index);

  vector<Type> N(n_ages);
  vector<Type> sel(n_ages);

  for (int a = 0; a < n_ages; ++a) {
    N(a) = R0 * exp(-M * Type(a));
    sel(a) = logistic((Type(a + 1) - Type(4.0)) / Type(0.8));
  }

  Type nll = Type(0.0);

  for (int y = 0; y < n_years; ++y) {
    Type vulnerable = Type(0.0);
    for (int a = 0; a < n_ages; ++a) {
      vulnerable += sel(a) * N(a);
    }

    nll -=
        dnorm(log(index_obs(y)), log(q) + log(vulnerable), sigma_index, true);
    nll -= dnorm(x(y), Type(0.0), sigma_R, true);

    vector<Type> N_next(n_ages);
    N_next(0) = R0 * exp(x(y));

    for (int a = 1; a < n_ages; ++a) {
      N_next(a) = N(a - 1) * exp(-M);
    }

    N_next(n_ages - 1) += N(n_ages - 1) * exp(-M);

    N = N_next;
  }

  return nll;
}
