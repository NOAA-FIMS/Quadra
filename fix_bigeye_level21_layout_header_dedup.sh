#!/usr/bin/env bash
set -euo pipefail

L21="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"
LAYOUT="$L21/diagnostics/bigeye_level21_parameter_layout.hpp"

if [[ ! -d "$L21" ]]; then
  echo "ERROR: missing $L21"
  exit 1
fi

ts="$(date +%Y%m%d_%H%M%S)"
for f in \
  "$L21/diagnostics/bigeye_initial_numbers_diagnostics.hpp" \
  "$L21/diagnostics/bigeye_age_comp_residual_diagnostics.hpp" \
  "$L21/diagnostics/bigeye_purse_seine_prediction_decomposition.hpp" \
  "$L21/diagnostics/bigeye_longline_prediction_decomposition.hpp" \
  "$L21/diagnostics/bigeye_level21_parameter_sanity_diagnostics.hpp" \
  "$L21/reports/bigeye_fit_reports.hpp"; do
  [[ -f "$f" ]] && cp "$f" "$f.before_layout_header_dedup.$ts"
done

cat > "$LAYOUT" <<'HPP'
#pragma once

#include "../quadra/bigeye_age_structured.hpp"

namespace pifsc_bigeye_tuna {
namespace level21_layout {

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

} // namespace level21_layout
} // namespace pifsc_bigeye_tuna
HPP

python3 - <<'PY'
from pathlib import Path
import re

root = Path("examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic")

files = [
    root / "diagnostics/bigeye_initial_numbers_diagnostics.hpp",
    root / "diagnostics/bigeye_age_comp_residual_diagnostics.hpp",
    root / "diagnostics/bigeye_purse_seine_prediction_decomposition.hpp",
    root / "diagnostics/bigeye_longline_prediction_decomposition.hpp",
    root / "diagnostics/bigeye_level21_parameter_sanity_diagnostics.hpp",
    root / "reports/bigeye_fit_reports.hpp",
]

layout_ns_block = re.compile(
    r"\n?namespace level21_layout \{\n"
    r"constexpr int kBaseFixed = 3;\n"
    r"constexpr int kMParamOffset = kBaseFixed;\n"
    r"constexpr int kMParams = 2;\n"
    r"constexpr int kLonglineSelOffset = kMParamOffset \+ kMParams;\n"
    r"constexpr int kLonglineSelDevs = kAges;\n"
    r"constexpr int kInitialDevOffset = kLonglineSelOffset \+ kLonglineSelDevs;\n"
    r"constexpr int kInitialDevs = kAges;\n"
    r"constexpr int kPurseSeineSelOffset = kInitialDevOffset \+ kInitialDevs;\n"
    r"constexpr int kPurseSeineSelDevs = kAges;\n"
    r"constexpr int kRecruitmentOffset = kPurseSeineSelOffset \+ kPurseSeineSelDevs;\n"
    r"\} // namespace level21_layout\n",
    re.M,
)

# Some files got bare namespace-level constants, not inside level21_layout.
bare_block = re.compile(
    r"\n?constexpr int kBaseFixed = 3;\n"
    r"constexpr int kMParamOffset = kBaseFixed;\n"
    r"constexpr int kMParams = 2;\n"
    r"constexpr int kLonglineSelOffset = kMParamOffset \+ kMParams;\n"
    r"constexpr int kLonglineSelDevs = kAges;\n"
    r"constexpr int kInitialDevOffset = kLonglineSelOffset \+ kLonglineSelDevs;\n"
    r"constexpr int kInitialDevs = kAges;\n"
    r"constexpr int kPurseSeineSelOffset = kInitialDevOffset \+ kInitialDevs;\n"
    r"constexpr int kPurseSeineSelDevs = kAges;\n"
    r"constexpr int kRecruitmentOffset = kPurseSeineSelOffset \+ kPurseSeineSelDevs;\n",
    re.M,
)

local_block = re.compile(
    r"\n?\s*constexpr int kBaseFixed = 3;\n"
    r"\s*constexpr int kMParamOffset = kBaseFixed;\n"
    r"\s*constexpr int kMParams = 2;\n"
    r"\s*constexpr int kLonglineSelOffset = kMParamOffset \+ kMParams;\n"
    r"\s*constexpr int kLonglineSelDevs = kAges;\n"
    r"\s*constexpr int kInitialDevOffset = kLonglineSelOffset \+ kLonglineSelDevs;\n"
    r"\s*constexpr int kInitialDevs = kAges;\n"
    r"\s*constexpr int kPurseSeineSelOffset = kInitialDevOffset \+ kInitialDevs;\n"
    r"\s*constexpr int kPurseSeineSelDevs = kAges;\n"
    r"\s*constexpr int kRecruitmentOffset = kPurseSeineSelOffset \+ kPurseSeineSelDevs;\n",
    re.M,
)

for p in files:
    if not p.exists():
        continue
    s = p.read_text()

    s = layout_ns_block.sub("\n", s)
    s = bare_block.sub("\n", s)

    # Remove duplicated function-local constants only if the shared header include will be added.
    s = local_block.sub("\n", s)

    if '#include "bigeye_level21_parameter_layout.hpp"' not in s:
        # Diagnostics are in same dir; reports are in reports and need ../diagnostics
        include = '#include "bigeye_level21_parameter_layout.hpp"\n'
        if "/reports/" in str(p):
            include = '#include "../diagnostics/bigeye_level21_parameter_layout.hpp"\n'
        # Insert after pragma once if present, else after last include.
        if "#pragma once" in s:
            s = s.replace("#pragma once\n", "#pragma once\n\n" + include, 1)
        else:
            includes = list(re.finditer(r"^#include[^\n]*\n", s, re.M))
            if includes:
                pos = includes[-1].end()
                s = s[:pos] + include + s[pos:]
            else:
                s = include + s

    # Qualify layout names if not already in using namespace.
    if "using namespace level21_layout;" not in s:
        if "namespace pifsc_bigeye_tuna {" in s:
            s = s.replace("namespace pifsc_bigeye_tuna {\n", "namespace pifsc_bigeye_tuna {\n\nusing namespace level21_layout;\n", 1)

    p.write_text(s)

# Objective should keep its function-local constants; do not include the diagnostics layout there.

PY

cat > audit_bigeye_level21_layout_dedup.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail
L21="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"

echo "== Shared layout header =="
cat "$L21/diagnostics/bigeye_level21_parameter_layout.hpp"

echo
echo "== Duplicate bare constant definitions scan =="
grep -R "constexpr int kBaseFixed\|constexpr int kMParamOffset\|constexpr int kLonglineSelOffset" -n \
  "$L21/diagnostics" "$L21/reports" "$L21/objective" || true

echo
echo "== Includes layout header =="
grep -R "bigeye_level21_parameter_layout.hpp\|using namespace level21_layout" -n \
  "$L21/diagnostics" "$L21/reports" || true
SH
chmod +x audit_bigeye_level21_layout_dedup.sh

echo "Installed Level 21 layout header de-dup fix."
echo
echo "Run:"
echo "  ./audit_bigeye_level21_layout_dedup.sh"
echo "  ./run_bigeye_level21_age_based_m_check.sh"
