#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator()() {
  DATA_VECTOR(catch_observed);
  DATA_VECTOR(index_observed);

  PARAMETER(log_r);
  PARAMETER(log_K);
  PARAMETER(log_q);
  PARAMETER(log_sigma_process);
  PARAMETER(log_sigma_index);
  PARAMETER(logit_B0_frac);

  PARAMETER_VECTOR(u);

  const int n = catch_observed.size();

  Type r = exp(log_r);
  Type K = exp(log_K);
  Type q = exp(log_q);
  Type sigma_process = exp(log_sigma_process);
  Type sigma_index = exp(log_sigma_index);
  Type B0_frac = Type(1.0) / (Type(1.0) + exp(-logit_B0_frac));

  vector<Type> log_biomass(n);
  vector<Type> biomass(n);
  vector<Type> index_predicted(n);

  biomass(0) = B0_frac * K;
  log_biomass(0) = log(biomass(0));

  Type nll = Type(0.0);

  for (int t = 0; t < n; ++t) {
    biomass(t) = exp(log_biomass(t));
    index_predicted(t) = q * biomass(t);

    nll -= dnorm(log(index_observed(t)),
                 log(index_predicted(t)),
                 sigma_index,
                 true);

    if (t < n - 1) {
      const Type production = r * biomass(t) * (Type(1.0) - biomass(t) / K);
      const Type deterministic_next =
          CppAD::CondExpGt(biomass(t) + production - catch_observed(t),
                           Type(1e-9),
                           biomass(t) + production - catch_observed(t),
                           Type(1e-9));

      log_biomass(t + 1) = log(deterministic_next) + u(t);

      nll -= dnorm(u(t), Type(0.0), sigma_process, true);
    }
  }

  ADREPORT(r);
  ADREPORT(K);
  ADREPORT(q);
  ADREPORT(sigma_process);
  ADREPORT(sigma_index);
  ADREPORT(B0_frac);
  ADREPORT(biomass);
  ADREPORT(index_predicted);

  return nll;
}
