from pathlib import Path
import re
import sys

p = Path("examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/objective/bigeye_quadra_objective.hpp")
s = p.read_text()
orig = s

backup = p.with_name(p.name + ".before_level21_baranov_from_current_context")
if not backup.exists():
    backup.write_text(s)

catch_block = (
"      T longline_catch_hat = T(0.0);\n"
"      T purse_seine_catch_hat = T(0.0);\n"
"      std::array<T, kAges> longline_catch_at_age{};\n"
"      std::array<T, kAges> purse_seine_catch_at_age{};\n"
"\n"
"      for (int a = 0; a < kAges; ++a) {\n"
"        const auto i = static_cast<std::size_t>(a);\n"
"        const T f_longline_a = fbar * sel_longline[i];\n"
"        const T f_purse_seine_a = fbar * sel_purse_seine[i];\n"
"        const T z_a = m_at_age[i] + f_longline_a + f_purse_seine_a;\n"
"        const T total_harvest = T(1.0) - exp_t(-z_a);\n"
"\n"
"        const T longline_harvest_rate = (f_longline_a / z_a) * total_harvest;\n"
"        const T purse_seine_harvest_rate = (f_purse_seine_a / z_a) * total_harvest;\n"
"\n"
"        longline_catch_at_age[i] = n[i] * T(weight[i]) * longline_harvest_rate;\n"
"        purse_seine_catch_at_age[i] = n[i] * T(weight[i]) * purse_seine_harvest_rate;\n"
"\n"
"        longline_catch_hat = longline_catch_hat + longline_catch_at_age[i];\n"
"        purse_seine_catch_hat = purse_seine_catch_hat + purse_seine_catch_at_age[i];\n"
"      }\n"
"\n"
)

if "longline_catch_at_age" not in s:
    anchor = "      nll = nll + T(0.5) * square_t(rec_dev / sigma_rec_dev);\n\n"
    if anchor not in s:
        sys.exit("Could not find rec_dev prior anchor")
    s = s.replace(anchor, anchor + catch_block, 1)

# Remove old total catch block if still present.
s = re.sub(
    r'\n      T total_catch_hat = T\(0\.0\);\n      for \(int a = 0; a < kAges; \+\+a\) \{.*?total_catch_hat = total_catch_hat .*?\n      \}\n',
    '\n',
    s,
    count=1,
    flags=re.S,
)

# Replace fleet catch split with fleet-specific catch.
s = re.sub(
    r'const T catch_hat = total_catch_hat \* T\(fleet_catch_share\(obs\.fleet\)\);',
    'const T catch_hat = is_longline ? longline_catch_hat : purse_seine_catch_hat;',
    s,
)
s = re.sub(
    r'const T catch_hat = .*fleet_catch_share\(obs\.fleet\).*;',
    'const T catch_hat = is_longline ? longline_catch_hat : purse_seine_catch_hat;',
    s,
)

# Replace selected-number age comp basis.
s = s.replace(
    "        T selected_numbers_sum = T(0.0);\n",
    "        T fleet_catch_at_age_sum = T(0.0);\n",
)

marker = "        std::array<T, kAges> pred_age_comp{};\n\n"
alias = (
"        const auto &fleet_catch_at_age =\n"
"            is_longline ? longline_catch_at_age : purse_seine_catch_at_age;\n"
"\n"
)
if marker in s and "const auto &fleet_catch_at_age" not in s:
    s = s.replace(marker, marker + alias, 1)

s = re.sub(
    r'pred_age_comp\[i\]\s*=\s*n\[i\]\s*\*\s*sel\[i\]\s*;',
    'pred_age_comp[i] = fleet_catch_at_age[i];',
    s,
)

s = re.sub(
    r'selected_numbers_sum\s*=\s*selected_numbers_sum\s*\+\s*pred_age_comp\[i\]\s*;',
    'fleet_catch_at_age_sum = fleet_catch_at_age_sum + pred_age_comp[i];',
    s,
)

s = re.sub(
    r'pred_age_comp\[i\]\s*=\s*pred_age_comp\[i\]\s*/\s*max_t\(selected_numbers_sum,\s*min_positive\)\s*;',
    'pred_age_comp[i] = pred_age_comp[i] / max_t(fleet_catch_at_age_sum, min_positive);',
    s,
)

required = [
    "longline_catch_at_age",
    "purse_seine_catch_at_age",
    "const T catch_hat = is_longline ? longline_catch_hat : purse_seine_catch_hat;",
    "pred_age_comp[i] = fleet_catch_at_age[i];",
    "fleet_catch_at_age_sum",
]
missing = [x for x in required if x not in s]
if missing:
    sys.exit("Patch validation failed; missing: " + ", ".join(missing))

bad = ["selected_numbers_sum", "total_catch_hat * T(fleet_catch_share"]
bad_found = [x for x in bad if x in s]
if bad_found:
    sys.exit("Patch validation failed; old patterns remain: " + ", ".join(bad_found))

if s == orig:
    sys.exit("No changes made; objective may already be patched")

p.write_text(s)
print("patched", p)
