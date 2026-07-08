#!/usr/bin/env bash
set -euo pipefail

L21="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"
OBJ="$L21/objective/bigeye_quadra_objective.hpp"
REPORT="$L21/reports/bigeye_fit_reports.hpp"
INIT_DIAG="$L21/diagnostics/bigeye_initial_numbers_diagnostics.hpp"
AGE_DIAG="$L21/diagnostics/bigeye_age_comp_residual_diagnostics.hpp"
LL_DECOMP="$L21/diagnostics/bigeye_longline_prediction_decomposition.hpp"
PS_DECOMP="$L21/diagnostics/bigeye_purse_seine_prediction_decomposition.hpp"
SANITY="$L21/diagnostics/bigeye_level21_parameter_sanity_diagnostics.hpp"
DRIVER="$L21/quadra/bigeye_level21_age_based_natural_mortality_diagnostic.cpp"

if [[ ! -d "$L21" ]]; then
  echo "ERROR: missing $L21"
  exit 1
fi

ts="$(date +%Y%m%d_%H%M%S)"
for f in "$OBJ" "$REPORT" "$INIT_DIAG" "$AGE_DIAG" "$LL_DECOMP" "$PS_DECOMP" "$SANITY" "$DRIVER"; do
  [[ -f "$f" ]] && cp "$f" "$f.before_level21_index_audit_fix.$ts"
done

python3 - <<'PY'
from pathlib import Path
import re

root = Path("examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic")

layout_block_cpp = """
namespace level21_layout {
constexpr int kBaseFixed = 3;
constexpr int kMParamOffset = kBaseFixed;
constexpr int kMParams = 2;
constexpr int kLonglineSelOffset = kMParamOffset + kMParams;
constexpr int kLonglineSelDevs = kAges;
constexpr int kInitialDevOffset = kLonglineSelOffset + kLonglineSelDevs;
constexpr int kInitialDevs = kAges;
constexpr int kPurseSeineSelOffset = kInitialDevOffset + kInitialDevs;
constexpr int kPurseSeineSelDevs = kAges;
constexpr int kRecruitmentOffset = kPurseSeineSelOffset + kPurseSeineSelDevs;
} // namespace level21_layout

""".lstrip()

def ensure_header_include(path: Path):
    s = path.read_text()
    if "namespace level21_layout" not in s:
        # Insert after namespace pifsc_bigeye_tuna { if present, otherwise after includes.
        if "namespace pifsc_bigeye_tuna {" in s:
            s = s.replace("namespace pifsc_bigeye_tuna {\n", "namespace pifsc_bigeye_tuna {\n\n" + layout_block_cpp, 1)
        else:
            # Fallback after last include.
            m = list(re.finditer(r"^#include[^\n]*\n", s, flags=re.M))
            if m:
                pos = m[-1].end()
                s = s[:pos] + "\nnamespace pifsc_bigeye_tuna {\n\n" + layout_block_cpp + "\n} // namespace pifsc_bigeye_tuna\n\n" + s[pos:]
    return s

def patch_layout_in_file(path: Path):
    if not path.exists():
        return
    s = path.read_text()

    # Add local layout constants if missing and file is in pifsc namespace.
    if "kMParamOffset" not in s:
        if "namespace pifsc_bigeye_tuna {" in s:
            s = s.replace("namespace pifsc_bigeye_tuna {\n", "namespace pifsc_bigeye_tuna {\n\n" + layout_block_cpp, 1)
        else:
            # Add function-local constants after the first function opening if possible.
            s = s.replace("{\n", "{\n  constexpr int kBaseFixed = 3;\n  constexpr int kMParamOffset = kBaseFixed;\n  constexpr int kMParams = 2;\n  constexpr int kLonglineSelOffset = kMParamOffset + kMParams;\n  constexpr int kLonglineSelDevs = kAges;\n  constexpr int kInitialDevOffset = kLonglineSelOffset + kLonglineSelDevs;\n  constexpr int kInitialDevs = kAges;\n  constexpr int kPurseSeineSelOffset = kInitialDevOffset + kInitialDevs;\n  constexpr int kPurseSeineSelDevs = kAges;\n  constexpr int kRecruitmentOffset = kPurseSeineSelOffset + kPurseSeineSelDevs;\n", 1)

    # Replace common stale hard-coded offsets.
    s = re.sub(r"\bfit\.par\[3 \+ a\]", "fit.par[kLonglineSelOffset + a]", s)
    s = re.sub(r"\bfit\.par\[3 \+ static_cast<int>\(a\)\]", "fit.par[kLonglineSelOffset + static_cast<int>(a)]", s)
    s = re.sub(r"\bfit\.par\[3 \+ static_cast<std::size_t>\(a\)\]", "fit.par[kLonglineSelOffset + static_cast<std::size_t>(a)]", s)

    # Level 20 had initial dev offset 13; Level 21 should use kInitialDevOffset.
    s = re.sub(r"\bfit\.par\[13 \+ a\]", "fit.par[kInitialDevOffset + a]", s)
    s = re.sub(r"\bfit\.par\[13 \+ static_cast<int>\(a\)\]", "fit.par[kInitialDevOffset + static_cast<int>(a)]", s)
    s = re.sub(r"\bfit\.par\[13 \+ static_cast<std::size_t>\(a\)\]", "fit.par[kInitialDevOffset + static_cast<std::size_t>(a)]", s)

    # Level 20 had purse-seine selectivity offset 23; Level 21 should use kPurseSeineSelOffset.
    s = re.sub(r"\bfit\.par\[23 \+ a\]", "fit.par[kPurseSeineSelOffset + a]", s)
    s = re.sub(r"\bfit\.par\[23 \+ static_cast<int>\(a\)\]", "fit.par[kPurseSeineSelOffset + static_cast<int>(a)]", s)
    s = re.sub(r"\bfit\.par\[23 \+ static_cast<std::size_t>\(a\)\]", "fit.par[kPurseSeineSelOffset + static_cast<std::size_t>(a)]", s)

    # Common named constants from Level 20.
    s = re.sub(r"constexpr int kLonglineSelOffset = kBaseFixed;", "constexpr int kLonglineSelOffset = kMParamOffset + kMParams;", s)
    s = re.sub(r"constexpr int kInitialDevOffset = kLonglineSelOffset \+ kLonglineSelDevs;", "constexpr int kInitialDevOffset = kLonglineSelOffset + kLonglineSelDevs;", s)
    s = re.sub(r"constexpr int kPurseSeineSelOffset = kInitialDevOffset \+ kInitialDevs;", "constexpr int kPurseSeineSelOffset = kInitialDevOffset + kInitialDevs;", s)

    # Patch stale direct summary extraction in report writers where base fixeds are still assumed.
    s = s.replace("const double log_m = std::log(0.45);", "const double log_m = std::log(0.45);")
    s = s.replace("fit.par[kInitialDevOffset + a]", "fit.par[kInitialDevOffset + a]")

    path.write_text(s)

# Objective layout is authoritative.
obj = root / "objective/bigeye_quadra_objective.hpp"
s = obj.read_text()
if "constexpr int kMParamOffset" not in s:
    s = s.replace(
        "    constexpr int kBaseFixed = 3;\n",
        "    constexpr int kBaseFixed = 3;\n"
        "    constexpr int kMParamOffset = kBaseFixed;\n"
        "    constexpr int kMParams = 2;\n",
        1,
    )
s = re.sub(r"    constexpr int kLonglineSelOffset = [^;\n]+;", "    constexpr int kLonglineSelOffset = kMParamOffset + kMParams;", s, count=1)
s = re.sub(r"    constexpr int kInitialDevOffset = [^;\n]+;", "    constexpr int kInitialDevOffset = kLonglineSelOffset + kLonglineSelDevs;", s, count=1)
s = re.sub(r"    constexpr int kPurseSeineSelOffset = [^;\n]+;", "    constexpr int kPurseSeineSelOffset = kInitialDevOffset + kInitialDevs;", s, count=1)
s = re.sub(r"    constexpr int kRecruitmentOffset = [^;\n]+;", "    constexpr int kRecruitmentOffset = kPurseSeineSelOffset + kPurseSeineSelDevs;", s, count=1)
obj.write_text(s)

for rel in [
    "reports/bigeye_fit_reports.hpp",
    "diagnostics/bigeye_initial_numbers_diagnostics.hpp",
    "diagnostics/bigeye_age_comp_residual_diagnostics.hpp",
    "diagnostics/bigeye_longline_prediction_decomposition.hpp",
    "diagnostics/bigeye_purse_seine_prediction_decomposition.hpp",
    "diagnostics/bigeye_level21_parameter_sanity_diagnostics.hpp",
]:
    patch_layout_in_file(root / rel)

# More targeted repairs for report/diagnostic files.
for rel in [
    "reports/bigeye_fit_reports.hpp",
    "diagnostics/bigeye_initial_numbers_diagnostics.hpp",
    "diagnostics/bigeye_age_comp_residual_diagnostics.hpp",
    "diagnostics/bigeye_longline_prediction_decomposition.hpp",
    "diagnostics/bigeye_purse_seine_prediction_decomposition.hpp",
]:
    p = root / rel
    if not p.exists():
        continue
    s = p.read_text()

    # Ensure M values are read correctly where needed.
    if "const double adult_m = 0.45;" not in s and ("m_at_age" in s or "const double m =" in s):
        s = s.replace("const double m = 0.45;", "const double adult_m = 0.45;\n  const double m = adult_m;", 1)

    # If a diagnostic builds m_at_age, make it use Level 21 M offsets.
    if "m_at_age" in s and "log_m_young_offset" not in s:
        insert = (
            "  const double log_m_young_offset = fit.par[kMParamOffset + 0];\n"
            "  const double log_m_old_offset = fit.par[kMParamOffset + 1];\n"
            "  const double adult_m = 0.45;\n"
            "  const double m_young = adult_m * std::exp(log_m_young_offset);\n"
            "  const double m_old = adult_m * std::exp(log_m_old_offset);\n"
            "  std::array<double, kAges> m_at_age{};\n"
            "  for (int a = 0; a < kAges; ++a) {\n"
            "    m_at_age[static_cast<std::size_t>(a)] =\n"
            "        (a < 3) ? m_young : ((a >= 7) ? m_old : adult_m);\n"
            "  }\n"
        )
        # Insert after fixed quantity declarations near top of function.
        idx = s.find("  const auto weight = default_weight_at_age();")
        if idx != -1:
            s = s[:idx] + insert + s[idx:]

    # Replace stale scalar M transitions in diagnostics with m_at_age when available.
    if "m_at_age" in s:
        s = s.replace("* std::exp(-m);", "* std::exp(-m_at_age[static_cast<std::size_t>(a - 1)]);")
        s = s.replace("/ (1.0 - std::exp(-m));", "/ (1.0 - std::exp(-m_at_age[static_cast<std::size_t>(kAges - 1)]));")
        s = s.replace("const double z_a = m + f_a;", "const double z_a = m_at_age[i] + f_a;")
        s = s.replace("const double z_prev = m + fbar * avg_sel_prev;", "const double z_prev = m_at_age[prev] + fbar * avg_sel_prev;")
        s = s.replace("const double z_last = m + fbar * avg_sel_last;", "const double z_last = m_at_age[last] + fbar * avg_sel_last;")

    p.write_text(s)

# Fix fit summary writer if it prints init devs using stale offset.
report = root / "reports/bigeye_fit_reports.hpp"
if report.exists():
    s = report.read_text()

    # Replace manually indexed output blocks if present.
    s = s.replace("fit.par[3 + a]", "fit.par[kLonglineSelOffset + a]")
    s = s.replace("fit.par[13 + a]", "fit.par[kInitialDevOffset + a]")
    s = s.replace("fit.par[23 + a]", "fit.par[kPurseSeineSelOffset + a]")

    # Add M rows to summary if absent.
    if "log_m_young_offset" not in s and "m_fixed" in s:
        marker = '  out << "m_fixed," << 0.45 << "\\n";\n'
        if marker in s:
            s = s.replace(
                marker,
                marker +
                '  const double log_m_young_offset = fit.par[kMParamOffset + 0];\n'
                '  const double log_m_old_offset = fit.par[kMParamOffset + 1];\n'
                '  out << "log_m_young_offset," << log_m_young_offset << "\\n";\n'
                '  out << "m_young," << 0.45 * std::exp(log_m_young_offset) << "\\n";\n'
                '  out << "log_m_old_offset," << log_m_old_offset << "\\n";\n'
                '  out << "m_old," << 0.45 * std::exp(log_m_old_offset) << "\\n";\n',
                1,
            )
    report.write_text(s)

PY

cat > audit_bigeye_level21_parameter_layout.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

L21="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"

echo "== Objective layout constants =="
grep -n "kBaseFixed\|kMParamOffset\|kMParams\|kLonglineSelOffset\|kInitialDevOffset\|kPurseSeineSelOffset\|kRecruitmentOffset" \
  "$L21/objective/bigeye_quadra_objective.hpp" | head -80

echo
echo "== Report/diagnostic stale hard-coded offset scan =="
grep -R "fit\.par\[[0-9][0-9]* *+ *a\]\|fit\.par\[[0-9][0-9]* *+ *static_cast" -n \
  "$L21/reports" "$L21/diagnostics" || true

echo
echo "== M-at-age usage scan =="
grep -R "m_at_age\|log_m_young_offset\|log_m_old_offset\|m_young\|m_old" -n \
  "$L21/objective" "$L21/reports" "$L21/diagnostics" | head -160 || true
SH
chmod +x audit_bigeye_level21_parameter_layout.sh

echo "Installed Level 21 parameter-layout audit/fix."
echo
echo "Run:"
echo "  ./audit_bigeye_level21_parameter_layout.sh"
echo "  ./run_bigeye_level21_age_based_m_check.sh"
