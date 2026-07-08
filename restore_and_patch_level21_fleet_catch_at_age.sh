#!/usr/bin/env bash
set -euo pipefail

OBJ="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/objective/bigeye_quadra_objective.hpp"
BAK="$OBJ.before_fleet_specific_catch_at_age_fix"

if [[ ! -f "$BAK" ]]; then
  echo "ERROR: expected backup not found:"
  echo "  $BAK"
  echo "Show first 40 lines of current file:"
  sed -n '1,40p' "$OBJ"
  exit 1
fi

cp "$BAK" "$OBJ"
echo "Restored Level 21 objective from:"
echo "  $BAK"

python3 - <<'PY'
from pathlib import Path

obj = Path("examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/objective/bigeye_quadra_objective.hpp")
s = obj.read_text()

old_total = """      T total_catch_hat = T(0.0);
      for (int a = 0; a < kAges; ++a) {
        const auto i = static_cast<std::size_t>(a);
        const T total_sel = sel_longline[i] + sel_purse_seine[i];
        const T f_a = fbar * total_sel;
        const T z_a = m_at_age[i] + f_a;
        const T harvest_rate = (f_a / z_a) * (T(1.0) - exp_t(-z_a));
        total_catch_hat = total_catch_hat + n[i] * T(weight[i]) * harvest_rate;
      }
"""

new_total = """      std::array<T, kAges> catch_at_age_longline{};
      std::array<T, kAges> catch_at_age_purse_seine{};
      T longline_catch_hat = T(0.0);
      T purse_seine_catch_hat = T(0.0);

      for (int a = 0; a < kAges; ++a) {
        const auto i = static_cast<std::size_t>(a);
        const T f_longline_a = fbar * sel_longline[i];
        const T f_purse_seine_a = fbar * sel_purse_seine[i];
        const T f_total_a = f_longline_a + f_purse_seine_a;
        const T z_a = m_at_age[i] + f_total_a;
        const T common_catch_factor =
            n[i] * (T(1.0) - exp_t(-z_a)) / max_t(z_a, min_positive);

        catch_at_age_longline[i] = common_catch_factor * f_longline_a;
        catch_at_age_purse_seine[i] = common_catch_factor * f_purse_seine_a;

        longline_catch_hat =
            longline_catch_hat + catch_at_age_longline[i] * T(weight[i]);
        purse_seine_catch_hat =
            purse_seine_catch_hat + catch_at_age_purse_seine[i] * T(weight[i]);
      }
"""

if old_total not in s:
    raise SystemExit("Could not find exact total_catch_hat block in restored Level 21 objective")
s = s.replace(old_total, new_total, 1)

old_age = """        T vulnerable_biomass = T(0.0);
        T selected_numbers_sum = T(0.0);
        std::array<T, kAges> pred_age_comp{};

        for (int a = 0; a < kAges; ++a) {
          const auto i = static_cast<std::size_t>(a);
          vulnerable_biomass = vulnerable_biomass + n[i] * T(weight[i]) * sel[i];
          pred_age_comp[i] = n[i] * sel[i];
          selected_numbers_sum = selected_numbers_sum + pred_age_comp[i];
        }
"""

new_age = """        T vulnerable_biomass = T(0.0);
        T catch_at_age_sum = T(0.0);
        std::array<T, kAges> pred_age_comp{};

        const auto &fleet_catch_at_age =
            is_longline ? catch_at_age_longline : catch_at_age_purse_seine;

        for (int a = 0; a < kAges; ++a) {
          const auto i = static_cast<std::size_t>(a);
          vulnerable_biomass = vulnerable_biomass + n[i] * T(weight[i]) * sel[i];
          pred_age_comp[i] = fleet_catch_at_age[i];
          catch_at_age_sum = catch_at_age_sum + pred_age_comp[i];
        }
"""

if old_age not in s:
    raise SystemExit("Could not find exact selected-numbers age-comp block in restored Level 21 objective")
s = s.replace(old_age, new_age, 1)

old_catch = "        const T catch_hat = total_catch_hat * T(fleet_catch_share(obs.fleet));"
new_catch = "        const T catch_hat = is_longline ? longline_catch_hat : purse_seine_catch_hat;"
if old_catch not in s:
    raise SystemExit("Could not find fleet_catch_share catch_hat line")
s = s.replace(old_catch, new_catch, 1)

old_norm = "          pred_age_comp[i] = pred_age_comp[i] / max_t(selected_numbers_sum, min_positive);"
new_norm = "          pred_age_comp[i] = pred_age_comp[i] / max_t(catch_at_age_sum, min_positive);"
if old_norm not in s:
    raise SystemExit("Could not find selected_numbers_sum normalization line")
s = s.replace(old_norm, new_norm, 1)

obj.write_text(s)
print("Applied surgical Level 21 fleet-specific catch-at-age patch.")
PY

echo
echo "== First 20 lines sanity check =="
sed -n '1,20p' "$OBJ"

echo
echo "== Check old broken patterns are gone =="
grep -n "total_catch_hat\|fleet_catch_share(obs.fleet)\|selected_numbers_sum" "$OBJ" || true

echo
echo "== Check new fleet-specific catch-at-age code =="
grep -n "catch_at_age_longline\|catch_at_age_purse_seine\|f_longline_a\|f_purse_seine_a\|catch_at_age_sum" "$OBJ"

echo
echo "Now run:"
echo "  ./run_bigeye_level21_age_based_m_check.sh"
