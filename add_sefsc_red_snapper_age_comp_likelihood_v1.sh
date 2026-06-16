#!/usr/bin/env bash
set -euo pipefail

python3 - <<'PY'
from pathlib import Path

p = Path("examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit.cpp")
s = p.read_text()

s = s.replace(
"""template <class T>
T logistic_selectivity_t""",
"""template <class T>
T age_comp_nll(const std::array<double, kAges>& observed,
               const std::array<T, kAges>& predicted,
               double effective_n,
               double floor = 1.0e-12) {
  T nll = T(0.0);
  for (int a = 0; a < kAges; ++a) {
    const auto i = static_cast<std::size_t>(a);
    const double obs = std::max(observed[i], 0.0);
    if (obs > 0.0) {
      nll = nll - T(effective_n * obs) * log_t(max_t(predicted[i], floor));
    }
  }
  return nll;
}

template <class T>
T logistic_selectivity_t""",
1)

s = s.replace(
"""    const T sigma_log_catch = T(0.15);
    const double min_positive = 1.0e-12;""",
"""    const T sigma_log_catch = T(0.15);
    const double age_comp_effective_n = 50.0;
    const double min_positive = 1.0e-12;""",
1)

s = s.replace(
"""      if (obs.catch_mt > 0.0) {
        const T z = (log_t(T(obs.catch_mt)) -
                     log_t(max_t(catch_hat, min_positive))) /
                    sigma_log_catch;
        nll = nll + T(0.5) * square_t(z);
      }

      std::array<T, kAges> next{};""",
"""      if (obs.catch_mt > 0.0) {
        const T z = (log_t(T(obs.catch_mt)) -
                     log_t(max_t(catch_hat, min_positive))) /
                    sigma_log_catch;
        nll = nll + T(0.5) * square_t(z);
      }

      std::array<T, kAges> pred_age_comp{};
      T selected_numbers_sum = T(0.0);
      for (int a = 0; a < kAges; ++a) {
        const auto i = static_cast<std::size_t>(a);
        pred_age_comp[i] = n[i] * selectivity[i];
        selected_numbers_sum = selected_numbers_sum + pred_age_comp[i];
      }
      for (int a = 0; a < kAges; ++a) {
        const auto i = static_cast<std::size_t>(a);
        pred_age_comp[i] =
            pred_age_comp[i] / max_t(selected_numbers_sum, min_positive);
      }

      nll = nll + age_comp_nll(obs.age_comp, pred_age_comp,
                               age_comp_effective_n, min_positive);

      std::array<T, kAges> next{};""",
1)

p.write_text(s)
PY

cat > examples/NMFS/sefsc_red_snapper/validation/age_composition_likelihood_checklist.md <<'MD'
# Age-Composition Likelihood Checklist

- [x] predicted selected age composition added
- [x] multinomial-style negative log likelihood added
- [x] fixed effective sample size added
- [ ] selectivity parameters estimated
- [ ] age-composition residuals written
- [ ] Dirichlet-multinomial alternative
MD
