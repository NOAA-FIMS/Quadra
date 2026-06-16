#!/usr/bin/env bash
set -euo pipefail

python3 - <<'PY'
from pathlib import Path

p = Path("examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit.cpp")
s = p.read_text()

if "write_selectivity_at_age" in s:
    print("Already installed")
    raise SystemExit(0)

anchor = """
void write_residual_diagnostics(
"""

helper = r'''
void write_selectivity_at_age(
    const std::string& path,
    const quadra::OptResult& fit)
{
  if (fit.par.size() < 5) {
    return;
  }

  const double a50 =
      1.0 + 9.0 / (1.0 + std::exp(-fit.par[3]));
  const double slope =
      std::exp(fit.par[4]);

  std::ofstream out(path);

  out << "age,selectivity\n";

  for (int age = 1; age <= kAges; ++age) {
    const double sel =
        1.0 / (1.0 + std::exp(-slope * (age - a50)));

    out << age << "," << sel << "\n";
  }
}

'''

if anchor not in s:
    raise SystemExit("Could not find residual diagnostics anchor")

s = s.replace(anchor, helper + "\n" + anchor, 1)

call_anchor = """
  write_residual_diagnostics(
      diagnostics_path,
      observations,
      fit);

  std::cout << "wrote:      " << diagnostics_path << std::endl;
"""

call_block = r'''
  const std::string selectivity_path =
      "examples/NMFS/sefsc_red_snapper/outputs/selectivity_at_age.csv";

  write_selectivity_at_age(
      selectivity_path,
      fit);

  std::cout << "wrote:      "
            << selectivity_path
            << std::endl;

'''

if call_anchor not in s:
    raise SystemExit("Could not find diagnostics call block")

s = s.replace(call_anchor,
              call_anchor + call_block,
              1)

p.write_text(s)
PY

echo
echo "Installed selectivity-at-age output."
echo
echo "Run:"
echo "  ./examples/NMFS/sefsc_red_snapper/run_red_snapper_quadra_fit.sh"
echo "  cat examples/NMFS/sefsc_red_snapper/outputs/selectivity_at_age.csv"
