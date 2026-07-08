#!/usr/bin/env bash
set -euo pipefail

L21="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"
DRIVER="$L21/quadra/bigeye_level21_age_based_natural_mortality_diagnostic.cpp"

if [[ ! -f "$DRIVER" ]]; then
  echo "ERROR: missing $DRIVER"
  exit 1
fi

ts="$(date +%Y%m%d_%H%M%S)"
cp "$DRIVER" "$DRIVER.before_fixed_gradient_diag_fix.$ts"

python3 - <<'PY'
from pathlib import Path
import re

driver = Path("examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/quadra/bigeye_level21_age_based_natural_mortality_diagnostic.cpp")
s = driver.read_text()

if "#include <algorithm>" not in s:
    if "#include <iostream>" in s:
        s = s.replace("#include <iostream>\n", "#include <iostream>\n#include <algorithm>\n", 1)
    else:
        s = "#include <algorithm>\n" + s

# Remove any previously inserted gradient helper block.
start = s.find("struct Level21GradientRow")
if start != -1:
    end = s.find("\nint main(", start)
    if end == -1:
        raise SystemExit("Could not find end of old gradient helper before main()")
    s = s[:start] + s[end+1:]

helper = r'''
struct Level21GradientRow {
  std::string name;
  double value = 0.0;
  double gradient = 0.0;
  double abs_gradient = 0.0;
};

inline std::vector<std::string> level21_fixed_parameter_names() {
  std::vector<std::string> names;
  names.push_back("log_r0");
  names.push_back("log_fbar");
  names.push_back("log_q_purse_seine");
  names.push_back("log_m_young_offset");
  names.push_back("log_m_old_offset");

  for (int a = 0; a < pifsc_bigeye_tuna::kAges; ++a) {
    names.push_back("logit_sel_longline_age_" + std::to_string(a + 1));
  }
  for (int a = 0; a < pifsc_bigeye_tuna::kAges; ++a) {
    names.push_back("init_log_number_dev_age_" + std::to_string(a + 1));
  }
  for (int a = 0; a < pifsc_bigeye_tuna::kAges; ++a) {
    names.push_back("logit_sel_purse_seine_age_" + std::to_string(a + 1));
  }

  return names;
}

inline void write_level21_gradient_by_parameter(
    const std::string &txt_path,
    const std::string &csv_path,
    const quadra::OptResult &fit) {
  const auto names = level21_fixed_parameter_names();

  std::vector<Level21GradientRow> rows;
  const std::size_t n =
      std::min({names.size(), fit.par.size(), fit.fixed_gradient.size()});
  rows.reserve(n);

  for (std::size_t i = 0; i < n; ++i) {
    Level21GradientRow row;
    row.name = names[i];
    row.value = fit.par[i];
    row.gradient = fit.fixed_gradient[i];
    row.abs_gradient = std::abs(row.gradient);
    rows.push_back(row);
  }

  std::sort(rows.begin(), rows.end(),
            [](const auto &a, const auto &b) {
              return a.abs_gradient > b.abs_gradient;
            });

  std::ofstream csv(csv_path);
  csv << "rank,name,value,fixed_gradient,abs_fixed_gradient\n";
  for (std::size_t i = 0; i < rows.size(); ++i) {
    csv << (i + 1) << "," << rows[i].name << "," << rows[i].value << ","
        << rows[i].gradient << "," << rows[i].abs_gradient << "\n";
  }

  std::ofstream txt(txt_path);
  txt << "Level 21 Gradient by Fixed Parameter\n";
  txt << "====================================\n\n";
  txt << "Purpose\n";
  txt << "-------\n";
  txt << "Identify which fixed-effect parameters dominate the final gradient norm.\n";
  txt << "This uses OptResult::fixed_gradient directly.\n\n";

  txt << "Summary\n";
  txt << "-------\n";
  txt << "n_named_fixed_parameters: " << names.size() << "\n";
  txt << "n_values_in_fit_par: " << fit.par.size() << "\n";
  txt << "n_values_in_fixed_gradient: " << fit.fixed_gradient.size() << "\n";
  txt << "n_parameters_checked: " << n << "\n";
  txt << "reported_grad_norm: " << fit.grad_norm << "\n";
  if (!rows.empty()) {
    txt << "max_abs_fixed_gradient: " << rows.front().abs_gradient << "\n";
    txt << "max_abs_fixed_gradient_parameter: " << rows.front().name << "\n";
  }
  txt << "\n";

  txt << "Top gradients\n";
  txt << "-------------\n";
  txt << "rank,name,value,fixed_gradient,abs_fixed_gradient\n";
  const std::size_t top_n = std::min<std::size_t>(25, rows.size());
  for (std::size_t i = 0; i < top_n; ++i) {
    txt << (i + 1) << "," << rows[i].name << "," << rows[i].value << ","
        << rows[i].gradient << "," << rows[i].abs_gradient << "\n";
  }
}
'''

pos = s.find("int main(")
if pos == -1:
    raise SystemExit("Could not find main()")
s = s[:pos] + helper + "\n" + s[pos:]

# Replace any older call variants with the correct signature.
s = re.sub(
    r'write_level21_gradient_by_parameter\(\s*'
    r'"examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_gradient_by_parameter\.txt",\s*'
    r'"examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_gradient_by_parameter\.csv",\s*'
    r'[^;]*?fit\);',
    'write_level21_gradient_by_parameter(\n'
    '        "examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_gradient_by_parameter.txt",\n'
    '        "examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_gradient_by_parameter.csv",\n'
    '        fit);',
    s,
    flags=re.S
)

if "bigeye_level21_gradient_by_parameter.txt" not in s:
    call = r'''
    write_level21_gradient_by_parameter(
        "examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_gradient_by_parameter.txt",
        "examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_gradient_by_parameter.csv",
        fit);
    std::cout << "wrote:      examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_gradient_by_parameter.txt\n";
    std::cout << "wrote:      examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_gradient_by_parameter.csv\n";
'''
    idx = s.rfind("return 0;")
    if idx == -1:
        raise SystemExit("Could not find return 0 insertion anchor")
    s = s[:idx] + call + "\n" + s[idx:]

driver.write_text(s)
PY

cat > run_bigeye_level21_gradient_by_parameter_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

./run_bigeye_level21_age_based_m_check.sh

echo
echo "== Level 21 top fixed-effect gradients =="
grep -n "Summary\|max_abs_fixed_gradient\|Top gradients" -A35 \
  examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_gradient_by_parameter.txt || true

echo
echo "== Compact top 15 fixed gradients =="
awk -F, 'NR==1 || (NR>1 && NR<=16)' \
  examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_gradient_by_parameter.csv || true
SH
chmod +x run_bigeye_level21_gradient_by_parameter_check.sh

echo "Installed Level 21 fixed_gradient-by-parameter diagnostic."
echo
echo "Run:"
echo "  ./run_bigeye_level21_gradient_by_parameter_check.sh"
