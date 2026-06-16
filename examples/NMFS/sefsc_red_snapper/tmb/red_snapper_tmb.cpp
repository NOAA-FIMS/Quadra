#include <TMB.hpp>

template <class Type> Type square(Type x) { return x * x; }

template <class Type>
Type logistic_selectivity(Type age, Type a50, Type slope) {
  return Type(1.0) / (Type(1.0) + exp(-slope * (age - a50)));
}

template <class Type> Type objective_function<Type>::operator()() {
  DATA_VECTOR(catch_obs);
  DATA_VECTOR(index_obs);
  DATA_MATRIX(age_comp_obs);

  PARAMETER(log_r0);
  PARAMETER(log_fbar);
  PARAMETER(log_q);
  PARAMETER(logit_sel_a50);
  PARAMETER(log_sel_slope);
  PARAMETER_VECTOR(log_rec_dev);

  const int n_years = catch_obs.size();
  const int n_ages = age_comp_obs.cols();

  Type r0 = exp(log_r0);
  Type m = Type(0.18);
  Type fbar = exp(log_fbar);
  Type q = exp(log_q);
  Type sel_a50 = Type(1.0) + Type(9.0) * invlogit(logit_sel_a50);
  Type sel_slope = exp(log_sel_slope);

  Type sigma_log_index = Type(0.20);
  Type sigma_log_catch = Type(0.15);
  Type sigma_rec_dev = Type(0.35);
  Type age_comp_effective_n = Type(2.0);
  Type min_positive = Type(1.0e-12);

  vector<Type> weight(n_ages);
  Type weight_values[10] = {Type(0.40), Type(0.85), Type(1.35), Type(1.95),
                            Type(2.60), Type(3.25), Type(3.85), Type(4.35),
                            Type(4.75), Type(5.05)};
  for (int a = 0; a < n_ages; ++a) {
    weight(a) = weight_values[a];
  }

  vector<Type> sel(n_ages);
  for (int a = 0; a < n_ages; ++a) {
    sel(a) = logistic_selectivity(Type(a + 1), sel_a50, sel_slope);
  }

  vector<Type> n(n_ages);
  n(0) = r0;
  for (int a = 1; a < n_ages; ++a) {
    n(a) = n(a - 1) * exp(-m);
  }
  n(n_ages - 1) = n(n_ages - 1) / (Type(1.0) - exp(-m));

  Type nll = Type(0.0);
  Type fixed_prior_nll = Type(0.0);
  Type rec_prior_nll = Type(0.0);
  Type index_nll = Type(0.0);
  Type catch_nll = Type(0.0);
  Type age_comp_nll = Type(0.0);

  fixed_prior_nll +=
      Type(0.5) * square((log_r0 - Type(std::log(1200.0))) / Type(1.0));
  fixed_prior_nll +=
      Type(0.5) * square((log_fbar - Type(std::log(0.025))) / Type(0.75));
  fixed_prior_nll +=
      Type(0.5) * square((log_q - Type(std::log(0.00005))) / Type(1.0));
  fixed_prior_nll += Type(0.5) * square((sel_a50 - Type(4.0)) / Type(0.75));
  fixed_prior_nll +=
      Type(0.5) * square((log_sel_slope - Type(std::log(1.2))) / Type(0.35));

  nll += fixed_prior_nll;

  for (int y = 0; y < n_years; ++y) {
    Type rec_dev = log_rec_dev(y);
    {
      Type term = Type(0.5) * square(rec_dev / sigma_rec_dev);
      rec_prior_nll += term;
      nll += term;
    }

    Type total_biomass = Type(0.0);
    Type catch_hat = Type(0.0);
    Type selected_sum = Type(0.0);
    vector<Type> pred_age_comp(n_ages);

    for (int a = 0; a < n_ages; ++a) {
      Type fa = fbar * sel(a);
      Type za = m + fa;
      total_biomass += n(a) * weight(a);
      catch_hat += n(a) * weight(a) * fa / za * (Type(1.0) - exp(-za));
      pred_age_comp(a) = n(a) * sel(a);
      selected_sum += pred_age_comp(a);
    }

    Type index_hat = q * total_biomass;

    if (index_obs(y) > 0.0) {
      Type z = (log(Type(index_obs(y))) -
                log(CppAD::CondExpGt(index_hat, min_positive, index_hat,
                                     min_positive))) /
               sigma_log_index;
      {
        Type term = Type(0.5) * z * z;
        index_nll += term;
        nll += term;
      }
    }

    if (catch_obs(y) > 0.0) {
      Type z = (log(Type(catch_obs(y))) -
                log(CppAD::CondExpGt(catch_hat, min_positive, catch_hat,
                                     min_positive))) /
               sigma_log_catch;
      {
        Type term = Type(0.5) * z * z;
        catch_nll += term;
        nll += term;
      }
    }

    for (int a = 0; a < n_ages; ++a) {
      pred_age_comp(a) =
          pred_age_comp(a) / CppAD::CondExpGt(selected_sum, min_positive,
                                              selected_sum, min_positive);
      Type obs = age_comp_obs(y, a);
      if (obs > 0.0) {
        {
          Type term = -age_comp_effective_n * obs *
                      log(CppAD::CondExpGt(pred_age_comp(a), min_positive,
                                           pred_age_comp(a), min_positive));
          age_comp_nll += term;
          nll += term;
        }
      }
    }

    vector<Type> next(n_ages);
    next.setZero();
    next(0) = r0 * exp(rec_dev);

    for (int a = 1; a < n_ages; ++a) {
      Type f_prev = fbar * sel(a - 1);
      Type z_prev = m + f_prev;
      next(a) = n(a - 1) * exp(-z_prev);
    }

    Type f_last = fbar * sel(n_ages - 1);
    Type z_last = m + f_last;
    next(n_ages - 1) += n(n_ages - 1) * exp(-z_last);

    n = next;
  }

  REPORT(fixed_prior_nll);
  REPORT(rec_prior_nll);
  REPORT(index_nll);
  REPORT(catch_nll);
  REPORT(age_comp_nll);
  return nll;
}
