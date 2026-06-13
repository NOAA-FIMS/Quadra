#!/usr/bin/env bash
set -euo pipefail

echo "== Add recruitment deviations as random effects to SEFSC red snapper Quadra fit =="

python3 - <<'PY'
from pathlib import Path

p = Path("examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit.cpp")
s = p.read_text()

s = s.replace(
'''      if (par.size() < 5)
      {
        throw std::runtime_error(
            "RedSnapperQuadraObjective expected parameters: log_r0, log_fbar, log_q, logit_sel_a50, log_sel_slope");
      }
''',
'''      if (par.size() < 5 + observations_.size())
      {
        throw std::runtime_error(
            "RedSnapperQuadraObjective expected parameters: log_r0, log_fbar, log_q, logit_sel_a50, log_sel_slope, log_rec_dev[year]");
      }
''',
1)

s = s.replace(
'''      const T sigma_log_index = T(0.20);
      const T sigma_log_catch = T(0.15);
''',
'''      const T sigma_log_index = T(0.20);
      const T sigma_log_catch = T(0.15);
      const T sigma_rec_dev = T(0.35);
''',
1)

old = '''      nll = nll + normal_prior(log_sel_slope, std::log(1.2), 0.35);

      for (const auto& obs : observations_) {
'''
new = '''      nll = nll + normal_prior(log_sel_slope, std::log(1.2), 0.35);

      for (std::size_t t = 0; t < observations_.size(); ++t) {
        const auto& obs = observations_[t];
        const T rec_dev = par[5 + t];
        nll = nll + T(0.5) * square_t(rec_dev / sigma_rec_dev);
'''
if old not in s:
    raise SystemExit("Could not find start of observation loop / prior block")
s = s.replace(old, new, 1)

s = s.replace(
'''      std::array<T, kAges> next{};
      next[0] = r0;
''',
'''      std::array<T, kAges> next{};
      next[0] = r0 * exp_t(rec_dev);
''',
1)

if 'log_rec_dev_' not in s:
    anchor = '''  params.add({"log_sel_slope", std::log(1.2), quadra::ParameterTransform::Identity, false});
'''
    insert = anchor + '''
  for (std::size_t t = 0; t < observations.size(); ++t) {
    params.add({"log_rec_dev_" + std::to_string(t + 1),
                0.0,
                quadra::ParameterTransform::Identity,
                true});
  }
'''
    if anchor not in s:
        raise SystemExit("Could not find log_sel_slope params.add anchor")
    s = s.replace(anchor, insert, 1)

p.write_text(s)
PY

cat > examples/NMFS/sefsc_red_snapper/validation/recruitment_deviation_laplace_checklist.md <<'MD'
# Recruitment-Deviation Laplace Checklist

- [x] one recruitment deviation random effect per fitted year
- [x] Gaussian recruitment-deviation prior
- [x] annual recruitment uses exp(log_rec_dev_t)
- [x] random effects passed through Quadra ParameterVector
- [ ] fitted recruitment deviations written
- [ ] random-effect trajectory written
- [ ] biomass/depletion uncertainty
- [ ] selected inverse diagnostics
MD

echo
echo "Patched recruitment deviations as random effects."
echo
echo "Run:"
echo "  ./examples/NMFS/sefsc_red_snapper/run_red_snapper_quadra_fit.sh"
echo "  cat examples/NMFS/sefsc_red_snapper/outputs/quadra_fit_summary.csv"
