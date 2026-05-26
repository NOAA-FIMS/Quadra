#include <TMB.hpp>

template <class Type> Type objective_function<Type>::operator()() {
  DATA_INTEGER(n_obs);

  PARAMETER(mu);
  PARAMETER(log_sigma);
  PARAMETER(u);

  Type sigma = exp(log_sigma);

  Type nll = Type(0.5) * u * u;

  for (int i = 0; i < n_obs; ++i) {
    Type yi = Type(1.0) + Type(0.10) * sin(Type(0.01) * Type(i));

    Type resid = (yi - (mu + u)) / sigma;

    nll += Type(0.5) * resid * resid + log_sigma;
  }

  nll += Type(0.5) * (mu / Type(10.0)) * (mu / Type(10.0));
  nll += Type(0.5) * (log_sigma / Type(2.0)) * (log_sigma / Type(2.0));

  return nll;
}
