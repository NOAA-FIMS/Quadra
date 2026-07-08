#!/usr/bin/env bash
set -euo pipefail

L17="examples/NMFS/pifsc_bigeye_tuna/level17_juvenile_mortality_diagnostic"

if [[ ! -d "$L17" ]]; then
  echo "ERROR: missing $L17. Run the Level 17 patch first."
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
for f in \
  "$L17/objective/bigeye_quadra_objective.hpp" \
  "$L17/diagnostics/bigeye_initial_numbers_diagnostics.hpp" \
  "$L17/diagnostics/bigeye_age_comp_residual_diagnostics.hpp"
do
  [[ -f "$f" ]] && cp "$f" "${f}.before_literal_newline_fix.${STAMP}"
done

python3 - <<'PY'
from pathlib import Path
import re

root = Path("examples/NMFS/pifsc_bigeye_tuna/level17_juvenile_mortality_diagnostic")

# 1) Fix objective literal \n fragments and stale M wiring.
p = root / "objective/bigeye_quadra_objective.hpp"
s = p.read_text()

# Convert literal backslash-n sequences that were inserted into source lines.
s = s.replace("\\n    const T log_juvenile_m_multiplier", "\n    const T log_juvenile_m_multiplier")
s = s.replace("\\n    const T juvenile_m", "\n    const T juvenile_m")
s = s.replace("\\n    const T sigma_log_juvenile_m_multiplier", "\n    const T sigma_log_juvenile_m_multiplier")
s = s.replace("\\n    nll = nll + T(0.5)", "\n    nll = nll + T(0.5)")

# Remove stale log_m/m lines left over from Level 16.
s = re.sub(r"\n\s*const T log_m = T\(std::log\(0\.45\)\);[^\n]*", "", s)
s = re.sub(r"\n\s*const T m = exp_t\(log_m\);", "", s)

# Ensure fixed-effect layout is Level 17.
s = s.replace("constexpr int kBaseFixed = 5;", "constexpr int kBaseFixed = 6;")

# Ensure juvenile multiplier exists after slope.
if "const T log_juvenile_m_multiplier = par[5];" not in s:
    s = s.replace(
        "const T log_sel_slope_longline = par[4];",
        "const T log_sel_slope_longline = par[4];\n    const T log_juvenile_m_multiplier = par[5];"
    )

# Ensure adult/juvenile M exists.
if "const T adult_m = T(0.45);" not in s:
    s = s.replace(
        "const T log_q_longline = T(std::log(0.00005));",
        "const T adult_m = T(0.45);\n    const T juvenile_m = adult_m * exp_t(log_juvenile_m_multiplier);\n    const T log_q_longline = T(std::log(0.00005));"
    )

# Ensure sigma prior exists.
if "sigma_log_juvenile_m_multiplier" not in s:
    s = s.replace(
        "const T sigma_ps_sel_dev = T(1.0);",
        "const T sigma_ps_sel_dev = T(1.0);\n    const T sigma_log_juvenile_m_multiplier = T(0.50);"
    )

# Ensure juvenile M prior exists after longline slope prior.
prior_line = "nll = nll + normal_prior(log_sel_slope_longline, std::log(1.2), 0.35);"
prior_add = "    nll = nll + T(0.5) * square_t(log_juvenile_m_multiplier / sigma_log_juvenile_m_multiplier);"
if prior_add not in s:
    s = s.replace(prior_line, prior_line + "\n" + prior_add)

# Ensure m_at_age block exists after weight declaration.
if "std::array<T, kAges> m_at_age" not in s:
    s = s.replace(
        "const auto weight = default_weight_at_age();",
        """const auto weight = default_weight_at_age();

    std::array<T, kAges> m_at_age{};
    for (int a = 0; a < kAges; ++a) {
      m_at_age[static_cast<std::size_t>(a)] =
          (a < 2) ? juvenile_m : adult_m;
    }"""
    )

# Replace any remaining fixed M uses in transitions.
s = s.replace("const T z_a = m + f_a;", "const T z_a = m_at_age[i] + f_a;")
s = s.replace("const T z_prev = m + fbar * avg_sel_prev;", "const T z_prev = m_at_age[prev] + fbar * avg_sel_prev;")
s = s.replace("const T z_last = m + fbar * avg_sel_last;", "const T z_last = m_at_age[last] + fbar * avg_sel_last;")
s = s.replace("] * exp_t(-m);", "] * exp_t(-m_at_age[static_cast<std::size_t>(a - 1)]);")
s = s.replace("/ (T(1.0) - exp_t(-m));", "/ (T(1.0) - exp_t(-m_at_age[static_cast<std::size_t>(kAges - 1)]));")

p.write_text(s)

# 2) Fix initial-numbers diagnostics literal \n and define adult/juvenile variables.
p = root / "diagnostics/bigeye_initial_numbers_diagnostics.hpp"
if p.exists():
    s = p.read_text()
    s = s.replace("\\n  const double juvenile_m", "\n  const double juvenile_m")
    s = s.replace("\\n  const double m = adult_m", "\n  const double m = adult_m")
    s = s.replace("constexpr int kBaseFixed = 5;", "constexpr int kBaseFixed = 6;")

    if "const double adult_m = 0.45;" not in s:
        # Prefer insertion before first use of equilibrium_n.
        s = s.replace(
            "const double log_r0 = fit.par[0];",
            "const double log_r0 = fit.par[0];"
        )
        m_anchor = "  const double m = 0.45;"
        if m_anchor in s:
            s = s.replace(
                m_anchor,
                "  const double adult_m = 0.45;\n  const double juvenile_m = adult_m * std::exp(fit.par[5]);\n  const double m = adult_m;"
            )
        else:
            # Insert after r0 if needed.
            s = s.replace(
                "  const double r0 = std::exp(log_r0);",
                "  const double r0 = std::exp(log_r0);\n  const double adult_m = 0.45;\n  const double juvenile_m = adult_m * std::exp(fit.par[5]);\n  const double m = adult_m;"
            )

    s = s.replace("] * std::exp(-m);", "] * std::exp(-((a - 1) < 2 ? juvenile_m : adult_m));")
    s = s.replace("/ (1.0 - std::exp(-m));", "/ (1.0 - std::exp(-adult_m));")
    p.write_text(s)

# 3) Fix age-comp residual diagnostics literal \n and define adult/juvenile variables.
p = root / "diagnostics/bigeye_age_comp_residual_diagnostics.hpp"
if p.exists():
    s = p.read_text()
    s = s.replace("\\n  const double juvenile_m", "\n  const double juvenile_m")
    s = s.replace("\\n  const double m = adult_m", "\n  const double m = adult_m")
    s = s.replace("constexpr int kBaseFixed = 5;", "constexpr int kBaseFixed = 6;")

    if "const double adult_m = 0.45;" not in s:
        m_anchor = "  const double m = 0.45;"
        if m_anchor in s:
            s = s.replace(
                m_anchor,
                "  const double adult_m = 0.45;\n  const double juvenile_m = adult_m * std::exp(fit.par[5]);\n  const double m = adult_m;"
            )
        else:
            s = s.replace(
                "  const double r0 = std::exp(fit.par[0]);",
                "  const double r0 = std::exp(fit.par[0]);\n  const double adult_m = 0.45;\n  const double juvenile_m = adult_m * std::exp(fit.par[5]);\n  const double m = adult_m;"
            )

    s = s.replace("] * std::exp(-m);", "] * std::exp(-((a - 1) < 2 ? juvenile_m : adult_m));")
    s = s.replace("/ (1.0 - std::exp(-m));", "/ (1.0 - std::exp(-adult_m));")
    s = s.replace("const double z_prev = m + fbar * avg_sel_prev;", "const double z_prev = ((a - 1) < 2 ? juvenile_m : adult_m) + fbar * avg_sel_prev;")
    s = s.replace("const double z_last = m + fbar * avg_sel_last;", "const double z_last = adult_m + fbar * avg_sel_last;")
    p.write_text(s)

PY

echo "Cleaned Level 17 literal newline and juvenile M wiring."
echo
echo "Run:"
echo "  ./inspect_bigeye_level17_juvenile_mortality.sh"
echo "  ./run_bigeye_level17_juvenile_mortality_check.sh"
