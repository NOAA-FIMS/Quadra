#!/usr/bin/env bash
set -euo pipefail

echo "== Add estimated selectivity parameters to SEFSC red snapper Quadra fit =="

python3 - <<'PY'
from pathlib import Path

p = Path("examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit.cpp")
s = p.read_text()

if "invlogit_t" not in s:
    old = """template <class T>
T log_t(const T& x) {
  using std::log;
  return log(x);
}

template <class T>
T max_t"""
    new = """template <class T>
T log_t(const T& x) {
  using std::log;
  return log(x);
}

template <class T>
T invlogit_t(const T& x) {
  return T(1.0) / (T(1.0) + exp_t(-x));
}

template <class T>
T max_t"""
    if old not in s:
        raise SystemExit("Could not find log_t block")
    s = s.replace(old, new, 1)

s = s.replace(
"""    if (par.size() < 3) {
      throw std::runtime_error(
          "RedSnapperQuadraObjective expected parameters: log_r0, log_fbar, log_q");
    }

    const T log_r0 = par[0];
    const T log_fbar = par[1];
    const T log_q = par[2];
""",
"""    if (par.size() < 5) {
      throw std::runtime_error(
          "RedSnapperQuadraObjective expected parameters: log_r0, log_fbar, log_q, logit_sel_a50, log_sel_slope");
    }

    const T log_r0 = par[0];
    const T log_fbar = par[1];
    const T log_q = par[2];
    const T logit_sel_a50 = par[3];
    const T log_sel_slope = par[4];
""",
1)

s = s.replace(
"""    const T q = exp_t(log_q);

    const T sigma_log_index = T(0.20);
""",
"""    const T q = exp_t(log_q);
    const T sel_a50 = T(1.0) + T(9.0) * invlogit_t(logit_sel_a50);
    const T sel_slope = exp_t(log_sel_slope);

    const T sigma_log_index = T(0.20);
""",
1)

s = s.replace(
"          logistic_selectivity_t(T(a + 1), T(4.0), T(1.2));",
"          logistic_selectivity_t(T(a + 1), sel_a50, sel_slope);",
1)

# Only add regularization if not already added.
if "normal_penalty" not in s:
    s = s.replace(
"""    T nll = T(0.0);

    for (const auto& obs : observations_) {
""",
"""    T nll = T(0.0);

    auto normal_penalty = [](const T& x, double mean, double sd) {
      const T z = (x - T(mean)) / T(sd);
      return T(0.5) * z * z;
    };

    nll = nll + normal_penalty(log_r0, std::log(1200.0), 1.0);
    nll = nll + normal_penalty(log_fbar, std::log(0.025), 0.75);
    nll = nll + normal_penalty(log_q, std::log(0.00005), 1.0);
    nll = nll + normal_penalty(sel_a50, 4.0, 2.0);
    nll = nll + normal_penalty(log_sel_slope, std::log(1.2), 1.0);

    for (const auto& obs : observations_) {
""",
1)

if "logit_sel_a50" not in s.split("void write_fit_summary", 1)[1].split("}", 1)[0]:
    s = s.replace(
"""  if (fit.par.size() >= 3) {
    out << "log_r0," << fit.par[0] << "\\n";
    out << "r0," << std::exp(fit.par[0]) << "\\n";
    out << "log_fbar," << fit.par[1] << "\\n";
    out << "fbar," << std::exp(fit.par[1]) << "\\n";
    out << "log_q," << fit.par[2] << "\\n";
    out << "q," << std::exp(fit.par[2]) << "\\n";
  }
""",
"""  if (fit.par.size() >= 3) {
    out << "log_r0," << fit.par[0] << "\\n";
    out << "r0," << std::exp(fit.par[0]) << "\\n";
    out << "log_fbar," << fit.par[1] << "\\n";
    out << "fbar," << std::exp(fit.par[1]) << "\\n";
    out << "log_q," << fit.par[2] << "\\n";
    out << "q," << std::exp(fit.par[2]) << "\\n";
  }

  if (fit.par.size() >= 5) {
    const double sel_a50 = 1.0 + 9.0 / (1.0 + std::exp(-fit.par[3]));
    const double sel_slope = std::exp(fit.par[4]);
    out << "logit_sel_a50," << fit.par[3] << "\\n";
    out << "sel_a50," << sel_a50 << "\\n";
    out << "log_sel_slope," << fit.par[4] << "\\n";
    out << "sel_slope," << sel_slope << "\\n";
  }
""",
1)

# Add selectivity mapping after trajectory params.log_q assignment in all helpers.
s = s.replace(
"""  params.log_r0 = fit.par[0];
  params.log_fbar = fit.par[1];
  params.log_q = fit.par[2];

  const auto rows =
""",
"""  params.log_r0 = fit.par[0];
  params.log_fbar = fit.par[1];
  params.log_q = fit.par[2];
  if (fit.par.size() >= 5) {
    params.sel_a50 = 1.0 + 9.0 / (1.0 + std::exp(-fit.par[3]));
    params.sel_slope = std::exp(fit.par[4]);
  }

  const auto rows =
""")

# Add parameters if missing.
if 'params.add({"logit_sel_a50"' not in s:
    s = s.replace(
"""  params.add({"log_r0", std::log(1200.0), quadra::ParameterTransform::Identity, false});
  params.add({"log_fbar", std::log(0.025), quadra::ParameterTransform::Identity, false});
  params.add({"log_q", std::log(0.00005), quadra::ParameterTransform::Identity, false});
""",
"""  params.add({"log_r0", std::log(1200.0), quadra::ParameterTransform::Identity, false});
  params.add({"log_fbar", std::log(0.025), quadra::ParameterTransform::Identity, false});
  params.add({"log_q", std::log(0.00005), quadra::ParameterTransform::Identity, false});
  params.add({"logit_sel_a50", 0.0, quadra::ParameterTransform::Identity, false});
  params.add({"log_sel_slope", std::log(1.2), quadra::ParameterTransform::Identity, false});
""",
1)

p.write_text(s)
PY

cat > examples/NMFS/sefsc_red_snapper/validation/selectivity_estimation_checklist.md <<'MD'
# Selectivity Estimation Checklist

- [x] estimated selectivity a50 fixed effect added
- [x] estimated selectivity slope fixed effect added
- [x] bounded a50 transform added
- [x] positive slope transform added
- [x] weak selectivity regularization added
- [x] fitted selectivity parameters written to summary
- [ ] age-composition residuals by age/year
- [ ] selectivity-at-age output
- [ ] Dirichlet-multinomial option
MD

echo
echo "Patched selectivity estimation."
echo
echo "Run:"
echo "  ./examples/NMFS/sefsc_red_snapper/run_red_snapper_quadra_fit.sh"
echo "  cat examples/NMFS/sefsc_red_snapper/outputs/quadra_fit_summary.csv"
