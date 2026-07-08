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
cp "$OBJ" "$OBJ.before_surgical_m_layout_fix.$ts"
[[ -f "$DIAG" ]] && cp "$DIAG" "$DIAG.before_surgical_m_layout_fix.$ts"

python3 - <<'PY'
from pathlib import Path
import re

obj = Path("examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/objective/bigeye_quadra_objective.hpp")
s = obj.read_text()

# Objective: insert M layout constants if missing.
if "constexpr int kMParamOffset" not in s:
    s = s.replace(
        "    constexpr int kBaseFixed = 3;\n",
        "    constexpr int kBaseFixed = 3;\n"
        "    constexpr int kMParamOffset = kBaseFixed;\n"
        "    constexpr int kMParams = 2;\n",
        1,
    )

# Objective: force downstream offsets to account for the two M params.
s = re.sub(
    r"    constexpr int kLonglineSelOffset = [^;\n]+;",
    "    constexpr int kLonglineSelOffset = kMParamOffset + kMParams;",
    s,
    count=1,
)
s = re.sub(
    r"    constexpr int kInitialDevOffset = [^;\n]+;",
    "    constexpr int kInitialDevOffset = kLonglineSelOffset + kLonglineSelDevs;",
    s,
    count=1,
)
s = re.sub(
    r"    constexpr int kPurseSeineSelOffset = [^;\n]+;",
    "    constexpr int kPurseSeineSelOffset = kInitialDevOffset + kInitialDevs;",
    s,
    count=1,
)
s = re.sub(
    r"    constexpr int kRecruitmentOffset = [^;\n]+;",
    "    constexpr int kRecruitmentOffset = kPurseSeineSelOffset + kPurseSeineSelDevs;",
    s,
    count=1,
)

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

    if "constexpr int kMParamOffset" not in s:
        s = s.replace(
            "  constexpr int kBaseFixed = 3;\n",
            "  constexpr int kBaseFixed = 3;\n"
            "  constexpr int kMParamOffset = kBaseFixed;\n"
            "  constexpr int kMParams = 2;\n",
            1,
        )

    s = re.sub(
        r"  constexpr int kLonglineSelOffset = [^;\n]+;",
        "  constexpr int kLonglineSelOffset = kMParamOffset + kMParams;",
        s,
        count=1,
    )
    s = re.sub(
        r"  constexpr int kInitialDevOffset = [^;\n]+;",
        "  constexpr int kInitialDevOffset = kLonglineSelOffset + kLonglineSelDevs;",
        s,
        count=1,
    )
    s = re.sub(
        r"  constexpr int kPurseSeineSelOffset = [^;\n]+;",
        "  constexpr int kPurseSeineSelOffset = kInitialDevOffset + kInitialDevs;",
        s,
        count=1,
    )
    s = re.sub(
        r"  constexpr int kRecruitmentOffset = [^;\n]+;",
        "  constexpr int kRecruitmentOffset = kPurseSeineSelOffset + kPurseSeineSelDevs;",
        s,
        count=1,
    )

    if "double m_age_penalty = 0.0;" not in s and "double init_penalty = 0.0;" in s:
        s = s.replace(
            "  double init_penalty = 0.0;\n",
            "  double m_age_penalty = 0.0;\n"
            "  double init_penalty = 0.0;\n",
            1,
        )

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

    if "age_based_m_prior_nll" not in s and 'txt << "initial_numbers_prior_nll," << init_penalty' in s:
        s = s.replace(
            '  txt << "initial_numbers_prior_nll," << init_penalty << "\\n";\n',
            '  txt << "age_based_m_prior_nll," << m_age_penalty << "\\n";\n'
            '  txt << "m_young," << 0.45 * std::exp(log_m_young_offset) << "\\n";\n'
            '  txt << "m_adult," << 0.45 << "\\n";\n'
            '  txt << "m_old," << 0.45 * std::exp(log_m_old_offset) << "\\n";\n'
            '  txt << "initial_numbers_prior_nll," << init_penalty << "\\n";\n',
            1,
        )

    diag.write_text(s)
PY

echo "== Objective constants =="
grep -n "kBaseFixed\|kMParamOffset\|kMParams\|kLonglineSelOffset\|kInitialDevOffset\|kPurseSeineSelOffset\|kRecruitmentOffset" "$OBJ" | head -80

echo
echo "== M parameters =="
grep -n "log_m_young_offset\|log_m_old_offset\|m_at_age\|sigma_log_m_age_offset" "$OBJ" | head -80

echo
echo "Installed surgical Level 21 M layout fix."
echo
echo "Run:"
echo "  ./run_bigeye_level21_age_based_m_check.sh"
