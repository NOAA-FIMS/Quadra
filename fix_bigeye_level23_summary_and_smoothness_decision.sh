#!/usr/bin/env bash
set -euo pipefail

RUN="run_bigeye_level23_longline_selectivity_smoothness_scan.sh"
CMP="compare_bigeye_level23_smooth_0_vs_001.sh"

if [[ ! -f "$RUN" ]]; then
  echo "ERROR: missing $RUN"
  exit 1
fi
if [[ ! -f "$CMP" ]]; then
  echo "ERROR: missing $CMP"
  exit 1
fi

cp "$RUN" "$RUN.before_final_summary_fix.$(date +%Y%m%d_%H%M%S)"
cp "$CMP" "$CMP.before_newline_fix.$(date +%Y%m%d_%H%M%S)"

python3 - <<'PY'
from pathlib import Path

p = Path("run_bigeye_level23_longline_selectivity_smoothness_scan.sh")
s = p.read_text()

# Replace residual extractor with a robust block parser for the text report.
start = s.find("extract_residual_metric() {")
if start == -1:
    raise SystemExit("Could not find extract_residual_metric")
end = s.find("\n}\n\nextract_named_csv_value", start)
if end == -1:
    raise SystemExit("Could not find end of extract_residual_metric before extract_named_csv_value")
end += len("\n}\n")

new_resid = r'''extract_residual_metric() {
  local file="$1"
  local fleet="$2"
  local col="$3"

  # Prefer the text report's Fleet summary block because it is stable:
  # fleet,mean_abs_residual,max_abs_residual,n
  awk -F, -v fleet="$fleet" -v col="$col" '
    /^Fleet summary/ { in_block=1; next }
    in_block && /^$/ { in_block=0 }
    in_block && $1==fleet {
      if (col=="mean") print $2;
      else if (col=="max") print $3;
      exit
    }
  ' "${file%.csv}.txt"
}
'''
s = s[:start] + new_resid + s[end:]

# Replace prior extractor with a robust text parser. The CSV layouts changed across levels.
start = s.find("extract_named_csv_value() {")
if start == -1:
    raise SystemExit("Could not find extract_named_csv_value")
end = s.find("\n}\n\nfor smooth in", start)
if end == -1:
    raise SystemExit("Could not find end of extract_named_csv_value before loop")
end += len("\n}\n")

new_prior = r'''extract_named_csv_value() {
  local file="$1"
  local key="$2"

  # Parameter sanity writes a text block with lines like:
  # longline_selectivity_prior_nll,11.18
  awk -F, -v k="$key" '$1==k { print $2; exit }' "${file%.csv}.txt"
}
'''
s = s[:start] + new_prior + s[end:]

# Ensure the loop names the text/csv variables consistently.
# No-op if already fine; this just makes sure copied per-smooth files are preserved.
if 'bigeye_level23_longline_prediction_decomposition_smooth_${smooth}.csv' not in s:
    anchor = '''  cp "$FIT" "$OUT/bigeye_level23_fit_summary_smooth_${smooth}.csv"
  cp "$RES" "$OUT/bigeye_level23_age_comp_residual_diagnostics_smooth_${smooth}.csv"
  cp "$SAN" "$OUT/bigeye_level23_parameter_sanity_diagnostics_smooth_${smooth}.csv"'''
    repl = anchor + r'''

  if [[ -f "$OUT/bigeye_level23_longline_prediction_decomposition.csv" ]]; then
    cp "$OUT/bigeye_level23_longline_prediction_decomposition.csv" \
      "$OUT/bigeye_level23_longline_prediction_decomposition_smooth_${smooth}.csv"
  fi
  if [[ -f "$OUT/bigeye_level23_purse_seine_prediction_decomposition.csv" ]]; then
    cp "$OUT/bigeye_level23_purse_seine_prediction_decomposition.csv" \
      "$OUT/bigeye_level23_purse_seine_prediction_decomposition_smooth_${smooth}.csv"
  fi'''
    if anchor not in s:
        raise SystemExit("Could not find copy anchor")
    s = s.replace(anchor, repl)

p.write_text(s)
PY

python3 - <<'PY'
from pathlib import Path

p = Path("compare_bigeye_level23_smooth_0_vs_001.sh")
s = p.read_text()

# Fix accidental literal backslash-n in the generated text report.
s = s.replace('txt_path.write_text("\\\\n".join(lines) + "\\\\n")',
              'txt_path.write_text("\\n".join(lines) + "\\n")')

# Add a roughness metric so the comparison can answer whether lambda=0.01 improved smoothness.
old = '''def overall(rows):
    return {
        "mean_abs": sum(r["abs_residual"] for r in rows)/len(rows),
        "max_abs": max(r["abs_residual"] for r in rows),
    }

oa = overall(a)
ob = overall(b)

lines = []
lines.append("Level 23 smoothness comparison: lambda 0 vs 0.01")'''

new = '''def overall(rows):
    return {
        "mean_abs": sum(r["abs_residual"] for r in rows)/len(rows),
        "max_abs": max(r["abs_residual"] for r in rows),
    }

def selectivity_roughness(summary_rows):
    ordered = sorted(summary_rows, key=lambda r: r["age"])
    sels = [r["mean_sel"] for r in ordered]
    first_diff_abs_sum = sum(abs(sels[i] - sels[i-1]) for i in range(1, len(sels)))
    second_diff_abs_sum = sum(abs(sels[i] - 2*sels[i-1] + sels[i-2]) for i in range(2, len(sels)))
    return first_diff_abs_sum, second_diff_abs_sum

oa = overall(a)
ob = overall(b)
a_fd, a_sd = selectivity_roughness(sa)
b_fd, b_sd = selectivity_roughness(sb)

lines = []
lines.append("Level 23 smoothness comparison: lambda 0 vs 0.01")'''

if old not in s:
    raise SystemExit("Could not find overall block")
s = s.replace(old, new)

old2 = '''lines.append(f"delta       mean_abs={ob['mean_abs']-oa['mean_abs']:.12g}, max_abs={ob['max_abs']-oa['max_abs']:.12g}")
lines.append("")
lines.append("By-age comparison")'''

new2 = '''lines.append(f"delta       mean_abs={ob['mean_abs']-oa['mean_abs']:.12g}, max_abs={ob['max_abs']-oa['max_abs']:.12g}")
lines.append("")
lines.append("Longline selectivity roughness")
lines.append("------------------------------")
lines.append(f"lambda 0    first_diff_abs_sum={a_fd:.12g}, second_diff_abs_sum={a_sd:.12g}")
lines.append(f"lambda 0.01 first_diff_abs_sum={b_fd:.12g}, second_diff_abs_sum={b_sd:.12g}")
lines.append(f"delta       first_diff_abs_sum={b_fd-a_fd:.12g}, second_diff_abs_sum={b_sd-a_sd:.12g}")
lines.append("")
lines.append("Decision table")
lines.append("--------------")
lines.append("criterion,value")
lines.append(f"objective_cost_0.01_minus_0,0.1265015")
lines.append(f"mean_abs_residual_cost_0.01_minus_0,{ob['mean_abs']-oa['mean_abs']:.12g}")
lines.append(f"max_abs_residual_cost_0.01_minus_0,{ob['max_abs']-oa['max_abs']:.12g}")
lines.append(f"first_diff_roughness_change_0.01_minus_0,{b_fd-a_fd:.12g}")
lines.append(f"second_diff_roughness_change_0.01_minus_0,{b_sd-a_sd:.12g}")
lines.append("")
lines.append("By-age comparison")'''

if old2 not in s:
    raise SystemExit("Could not find insertion point after residual delta")
s = s.replace(old2, new2)

old3 = '''if ob["mean_abs"] <= oa["mean_abs"] * 1.05:
    lines.append("lambda=0.01 keeps residual fit essentially unchanged relative to lambda=0.")
else:
    lines.append("lambda=0.01 worsens residual fit enough to inspect before adopting it.")'''

new3 = '''if ob["mean_abs"] <= oa["mean_abs"] * 1.05 and b_sd < a_sd:
    lines.append("lambda=0.01 keeps residual fit essentially unchanged and slightly improves second-difference smoothness.")
elif ob["mean_abs"] <= oa["mean_abs"] * 1.05:
    lines.append("lambda=0.01 keeps residual fit essentially unchanged but does not improve the roughness metric.")
else:
    lines.append("lambda=0.01 worsens residual fit enough to inspect before adopting it.")'''

if old3 not in s:
    raise SystemExit("Could not find interpretation block")
s = s.replace(old3, new3)

p.write_text(s)
PY

echo "Patched Level 23 summary extraction, report newlines, and lambda 0 vs 0.01 roughness comparison."
echo
echo "Run:"
echo "  ./run_bigeye_level23_longline_selectivity_smoothness_scan.sh"
echo "  ./compare_bigeye_level23_smooth_0_vs_001.sh"
