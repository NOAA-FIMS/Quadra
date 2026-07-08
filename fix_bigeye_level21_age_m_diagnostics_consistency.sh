#!/usr/bin/env bash
set -euo pipefail

L21="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"
DIAG="$L21/diagnostics"

if [[ ! -d "$L21" ]]; then
  echo "ERROR: missing $L21"
  exit 1
fi

python3 - <<'PY'
from pathlib import Path
import re

L21 = Path("examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic")
diag = L21 / "diagnostics"

helper = diag / "bigeye_level21_m_at_age_helpers.hpp"
helper.write_text('''#pragma once

#include <array>
#include <cmath>
#include <vector>
#include "../objective/bigeye_quadra_objective.hpp"

namespace pifsc_bigeye_tuna {
namespace level21_m_helpers {

inline constexpr int kBaseFixed = 3;
inline constexpr int kMParamOffset = kBaseFixed;
inline constexpr int kMParams = 2;
inline constexpr int kLonglineSelOffset = kMParamOffset + kMParams;
inline constexpr int kLonglineSelDevs = kAges;
inline constexpr int kInitialDevOffset = kLonglineSelOffset + kLonglineSelDevs;
inline constexpr int kInitialDevs = kAges;
inline constexpr int kPurseSeineSelOffset = kInitialDevOffset + kInitialDevs;
inline constexpr int kPurseSeineSelDevs = kAges;
inline constexpr int kRecruitmentOffset = kPurseSeineSelOffset + kPurseSeineSelDevs;

inline std::array<double, kAges>
m_at_age_from_level21_par(const std::vector<double> &par) {
  const double adult_m = 0.45;
  double log_m_young_offset = 0.0;
  double log_m_old_offset = 0.0;

  if (par.size() > static_cast<std::size_t>(kMParamOffset + 0)) {
    log_m_young_offset = par[static_cast<std::size_t>(kMParamOffset + 0)];
  }
  if (par.size() > static_cast<std::size_t>(kMParamOffset + 1)) {
    log_m_old_offset = par[static_cast<std::size_t>(kMParamOffset + 1)];
  }

  const double m_young = adult_m * std::exp(log_m_young_offset);
  const double m_old = adult_m * std::exp(log_m_old_offset);

  std::array<double, kAges> m_at_age{};
  for (int a = 0; a < kAges; ++a) {
    if (a < 2) {
      m_at_age[static_cast<std::size_t>(a)] = m_young;
    } else if (a >= 7) {
      m_at_age[static_cast<std::size_t>(a)] = m_old;
    } else {
      m_at_age[static_cast<std::size_t>(a)] = adult_m;
    }
  }
  return m_at_age;
}

inline double adult_m() { return 0.45; }

} // namespace level21_m_helpers
} // namespace pifsc_bigeye_tuna
''')

def backup(p: Path):
    if p.exists():
        b = p.with_name(p.name + ".before_level21_m_diag_fix")
        if not b.exists():
            b.write_text(p.read_text())

def ensure_include(s: str) -> str:
    inc = '#include "bigeye_level21_m_at_age_helpers.hpp"\n'
    if inc in s:
        return s
    if "#pragma once\n" in s:
        return s.replace("#pragma once\n", "#pragma once\n" + inc, 1)
    return inc + s

def patch_file(p: Path):
    if not p.exists():
        return
    backup(p)
    s = p.read_text()
    s = ensure_include(s)

    # Remove duplicated namespace-scope layout constants, if a previous patch inserted them.
    s = re.sub(
        r'\nconstexpr int kBaseFixed = 3;\n'
        r'constexpr int kMParamOffset = kBaseFixed;\n'
        r'constexpr int kMParams = 2;\n'
        r'constexpr int kLonglineSelOffset = kMParamOffset \+ kMParams;\n'
        r'constexpr int kLonglineSelDevs = kAges;\n'
        r'constexpr int kInitialDevOffset = kLonglineSelOffset \+ kLonglineSelDevs;\n'
        r'constexpr int kInitialDevs = kAges;\n'
        r'constexpr int kPurseSeineSelOffset = kInitialDevOffset \+ kInitialDevs;\n'
        r'constexpr int kPurseSeineSelDevs = kAges;\n'
        r'constexpr int kRecruitmentOffset = kPurseSeineSelOffset \+ kPurseSeineSelDevs;\n',
        '\n',
        s
    )

    # Qualify layout constants to the shared helper namespace.
    names = [
        "kMParamOffset", "kMParams", "kLonglineSelOffset", "kLonglineSelDevs",
        "kInitialDevOffset", "kInitialDevs", "kPurseSeineSelOffset",
        "kPurseSeineSelDevs", "kRecruitmentOffset"
    ]
    for name in names:
        s = re.sub(r'(?<!level21_m_helpers::)\b' + name + r'\b',
                   'level21_m_helpers::' + name, s)

    # Make a local m_at_age vector available in diagnostics/report helpers.
    if "m_at_age_from_level21_par(fit.par)" not in s:
        if re.search(r'const double m\s*=', s):
            s = re.sub(
                r'const double m\s*=\s*[^;]+;',
                'const auto m_at_age = level21_m_helpers::m_at_age_from_level21_par(fit.par);\n  const double m = level21_m_helpers::adult_m();',
                s,
                count=1
            )
        elif re.search(r'const double adult_m\s*=\s*0\.45;', s):
            s = re.sub(
                r'const double adult_m\s*=\s*0\.45;',
                'const auto m_at_age = level21_m_helpers::m_at_age_from_level21_par(fit.par);\n  const double adult_m = level21_m_helpers::adult_m();',
                s,
                count=1
            )

    # Convert common scalar-M diagnostic formulas to age-specific M.
    s = s.replace('const double z_a = m + f_a;', 'const double z_a = m_at_age[i] + f_a;')
    s = s.replace('const double z_prev = m + fbar * avg_sel_prev;', 'const double z_prev = m_at_age[prev] + fbar * avg_sel_prev;')
    s = s.replace('const double z_last = m + fbar * avg_sel_last;', 'const double z_last = m_at_age[last] + fbar * avg_sel_last;')
    s = s.replace('const double z_a = adult_m + f_a;', 'const double z_a = m_at_age[i] + f_a;')
    s = s.replace('const double z_prev = adult_m + fbar * avg_sel_prev;', 'const double z_prev = m_at_age[prev] + fbar * avg_sel_prev;')
    s = s.replace('const double z_last = adult_m + fbar * avg_sel_last;', 'const double z_last = m_at_age[last] + fbar * avg_sel_last;')

    # Be conservative with equilibrium initialization replacements.
    s = s.replace(
        'n[static_cast<std::size_t>(a)] =\n          n[static_cast<std::size_t>(a - 1)] * std::exp(-m);',
        'n[static_cast<std::size_t>(a)] =\n          n[static_cast<std::size_t>(a - 1)] * std::exp(-m_at_age[static_cast<std::size_t>(a - 1)]);'
    )
    s = s.replace(
        'n[static_cast<std::size_t>(kAges - 1)] =\n        n[static_cast<std::size_t>(kAges - 1)] / (1.0 - std::exp(-m));',
        'n[static_cast<std::size_t>(kAges - 1)] =\n        n[static_cast<std::size_t>(kAges - 1)] /\n        (1.0 - std::exp(-m_at_age[static_cast<std::size_t>(kAges - 1)]));'
    )

    p.write_text(s)

for filename in [
    "bigeye_age_comp_residual_diagnostics.hpp",
    "bigeye_initial_numbers_diagnostics.hpp",
    "bigeye_longline_prediction_decomposition.hpp",
    "bigeye_purse_seine_prediction_decomposition.hpp",
    "bigeye_level20_parameter_sanity_diagnostics.hpp",
]:
    patch_file(diag / filename)

run = Path("run_bigeye_level21_age_based_m_check.sh")
if run.exists():
    backup(run)
    s = run.read_text()
    marker = "== Level 21 M diagnostic consistency =="
    if marker not in s:
        s += '''
echo
echo "== Level 21 M diagnostic consistency =="
grep -R "m_at_age_from_level21_par\|m_at_age" -n \
  examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/diagnostics \
  | head -80
'''
        run.write_text(s)

print("Installed Level 21 age-M diagnostic consistency patch.")
print("Created diagnostics/bigeye_level21_m_at_age_helpers.hpp")
PY

echo
echo "Run:"
echo "  ./run_bigeye_level21_age_based_m_check.sh"
