#!/usr/bin/env bash
set -euo pipefail

L21="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"
OBJ="$L21/objective/bigeye_quadra_objective.hpp"
DIAG="$L21/diagnostics/bigeye_level21_parameter_sanity_diagnostics.hpp"

if [[ ! -f "$OBJ" ]]; then
  echo "ERROR: missing $OBJ"
  exit 1
fi

ts="$(date +%Y%m%d_%H%M%S)"
cp "$OBJ" "$OBJ.before_level21_m_layout_fix.$ts"
[[ -f "$DIAG" ]] && cp "$DIAG" "$DIAG.before_level21_m_layout_fix.$ts"

python3 - <<'PY'
from pathlib import Path
import re

obj = Path("examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/objective/bigeye_quadra_objective.hpp")
s = obj.read_text()

# Replace whatever Level 20/21 layout currently exists with the intended Level 21 layout.
pattern = re.compile(
    r"    constexpr int kBaseFixed = 3;\n"
    r"(?:    constexpr int kMParamOffset = kBaseFixed;\n"
    r"    constexpr int kMParams = 2;\n)?"
    r"    constexpr int kLonglineSelOffset = [^\n]+;\n"
    r"    constexpr int kLonglineSelDevs = kAges;\n"
    r"    constexpr int kInitialDevOffset = [^\n]+;\n"
    r"    constexpr int kInitialDevs = kAges;\n"
    r"    constexpr int kPurseSeineSelOffset = [^\n]+;\n"
    r"    constexpr int kPurseSeineSelDevs = kAges;\n"
    r"    constexpr int kRecruitmentOffset = [^\n]+;\n"
)
replacement = (
    "    constexpr int kBaseFixed = 3;\n"
    "    constexpr int kMParamOffset = kBaseFixed;\n"
    "    constexpr int kMParams = 2;\n"
    "    constexpr int kLonglineSelOffset = kMParamOffset + kMParams;\n"
    "    constexpr int kLonglineSelDevs = kAges;\n"
    "    constexpr int kInitialDevOffset = kLonglineSelOffset + kLonglineSelDevs;\n"
    "    constexpr int kInitialDevs = kAges;\n"
    "    constexpr int kPurseSeineSelOffset = kInitialDevOffset + kInitialDevs;\n"
    "    constexpr int kPurseSeineSelDevs = kAges;\n"
    "    constexpr int kRecruitmentOffset = kPurseSeineSelOffset + kPurseSeineSelDevs;\n"
)
s, n = pattern.subn(replacement, s, count=1)
if n != 1:
    raise SystemExit("Could not replace Level 21 objective layout constants")

# Ensure error message matches the actual layout.
s = re.sub(
    r'"Level 21 expected[^"]+"',
    '"Level 21 expected 3 base fixed effects, 2 age-based M parameters, longline selectivity logits, initial number deviations, purse-seine age selectivity logits, plus recruitment deviations"',
    s,
    count=1,
)

obj.write_text(s)

diag = Path("examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/diagnostics/bigeye_level21_parameter_sanity_diagnostics.hpp")
if diag.exists():
    s = diag.read_text()

    pattern = re.compile(
        r"  constexpr int kBaseFixed = 3;\n"
        r"(?:  constexpr int kMParamOffset = kBaseFixed;\n"
        r"  constexpr int kMParams = 2;\n)?"
        r"  constexpr int kLonglineSelOffset = [^\n]+;\n"
        r"  constexpr int kLonglineSelDevs = kAges;\n"
        r"  constexpr int kInitialDevOffset = [^\n]+;\n"
        r"  constexpr int kInitialDevs = kAges;\n"
        r"  constexpr int kPurseSeineSelOffset = [^\n]+;\n"
        r"  constexpr int kPurseSeineSelDevs = kAges;\n"
        r"  constexpr int kRecruitmentOffset = [^\n]+;\n"
    )
    replacement = (
        "  constexpr int kBaseFixed = 3;\n"
        "  constexpr int kMParamOffset = kBaseFixed;\n"
        "  constexpr int kMParams = 2;\n"
        "  constexpr int kLonglineSelOffset = kMParamOffset + kMParams;\n"
        "  constexpr int kLonglineSelDevs = kAges;\n"
        "  constexpr int kInitialDevOffset = kLonglineSelOffset + kLonglineSelDevs;\n"
        "  constexpr int kInitialDevs = kAges;\n"
        "  constexpr int kPurseSeineSelOffset = kInitialDevOffset + kInitialDevs;\n"
        "  constexpr int kPurseSeineSelDevs = kAges;\n"
        "  constexpr int kRecruitmentOffset = kPurseSeineSelOffset + kPurseSeineSelDevs;\n"
    )
    s, n = pattern.subn(replacement, s, count=1)
    if n != 1:
        print("WARNING: could not replace diagnostic layout constants; continuing")

    # Ensure M penalty variables exist if the earlier patch partially missed them.
    if "double m_age_penalty = 0.0;" not in s and "double init_penalty = 0.0;" in s:
        s = s.replace("  double init_penalty = 0.0;\n",
                      "  double m_age_penalty = 0.0;\n  double init_penalty = 0.0;\n", 1)

    if "const double sigma_log_m_age_offset = 0.35;" not in s:
        marker = "  for (int a = 0; a < kAges; ++a) {\n    const double v = fit.par[kInitialDevOffset + a];\n"
        if marker in s:
            s = s.replace(
                marker,
                "  const double sigma_log_m_age_offset = 0.35;\n"
                "  const double log_m_young_offset = fit.par[kMParamOffset + 0];\n"
                "  const double log_m_old_offset = fit.par[kMParamOffset + 1];\n"
                "  m_age_penalty += 0.5 * sq(log_m_young_offset / sigma_log_m_age_offset);\n"
                "  m_age_penalty += 0.5 * sq(log_m_old_offset / sigma_log_m_age_offset);\n\n"
                + marker,
                1,
            )

    diag.write_text(s)
PY

echo "== Level 21 objective layout check =="
grep -n "kBaseFixed\|kMParamOffset\|kLonglineSelOffset\|kInitialDevOffset\|kPurseSeineSelOffset\|kRecruitmentOffset\|log_m_young_offset\|log_m_old_offset" "$OBJ" | head -80

echo
echo "Installed Level 21 M layout fix."
echo
echo "Run:"
echo "  ./run_bigeye_level21_age_based_m_check.sh"
