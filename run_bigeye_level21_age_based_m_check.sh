#!/usr/bin/env bash
set -euo pipefail

echo "== O3 build Bigeye Level 21 age-based natural mortality diagnostic =="
mkdir -p build/examples
c++ ${CXXFLAGS:-} -std=c++17 -O3 -flto \
  -DBIGEYE_LL_SEL_SIGMA=1.75 \
  -DQUADRA_AUDIT_LOGDET_COMPONENTS \
  -DQUADRA_AUDIT_JOINT_LOGDET_ASSEMBLY \
  -DQUADRA_AUDIT_LOGDET_BACKEND_VALUE \
  -I. -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/quadra/bigeye_level21_age_based_natural_mortality_diagnostic.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level21_age_based_m_check

echo
echo "== Run Bigeye Level 21 age-based natural mortality diagnostic =="
./build/examples/pifsc_bigeye_level21_age_based_m_check

echo
echo "== Level 21 fit summary key rows =="
grep -E "objective|grad_norm|converged|log_m_young_offset|log_m_old_offset|sel_longline_age_|init_number_multiplier_age_" \
  examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_fit_summary.csv \
  | head -100 || true

echo
echo "== Level 21 parameter sanity prior blocks =="
grep -n "Prior penalty by block" -A16 \
  examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_parameter_sanity_diagnostics.txt

echo
echo "== Level 21 objective-consistent residual summary =="
grep -n "Fleet summary" -A8 \
  examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_age_comp_residual_diagnostics.txt

echo
echo "== Level 20 sigma=1.75 vs Level 21 compact comparison =="
python3 - <<'PY'
from pathlib import Path
import csv

def read_summary(path):
    d = {}
    if not Path(path).exists():
        return d
    with open(path) as f:
        for row in csv.reader(f):
            if len(row) >= 2:
                d[row[0]] = row[1]
    return d

def fleet_summary(path):
    d = {}
    if not Path(path).exists():
        return d
    lines = Path(path).read_text().splitlines()
    in_block = False
    for line in lines:
        if line.startswith("fleet,mean_abs_residual"):
            in_block = True
            continue
        if in_block:
            if not line.strip():
                break
            parts = line.split(",")
            if len(parts) >= 4:
                d[parts[0]] = {"mean_abs": parts[1], "max_abs": parts[2], "n": parts[3]}
    return d

def prior_blocks(path):
    d = {}
    if not Path(path).exists():
        return d
    lines = Path(path).read_text().splitlines()
    in_block = False
    for line in lines:
        if line.startswith("Prior penalty by block"):
            in_block = True
            continue
        if in_block:
            if not line.strip() or line.startswith("Initial-number"):
                break
            if "," in line:
                k, v = line.split(",", 1)
                d[k] = v
    return d

s20 = {}
f20 = {}
p20 = {}
scan = Path("examples/NMFS/pifsc_bigeye_tuna/level20_longline_selectivity_regularization_scan/outputs/bigeye_level20_ll_selectivity_wide_sigma_scan_summary.csv")
if scan.exists():
    for r in csv.DictReader(scan.open()):
        if r.get("sigma") == "1.75":
            s20 = {"objective": r["objective"], "grad_norm": r["grad_norm"]}
            f20 = {
                "longline": {"mean_abs": r["longline_mean_abs"], "max_abs": r["longline_max_abs"]},
                "purse_seine": {"mean_abs": r["purse_seine_mean_abs"], "max_abs": r["purse_seine_max_abs"]},
            }
            p20 = {
                "longline_selectivity_prior_nll": r["longline_selectivity_prior_nll"],
                "initial_numbers_prior_nll": r["initial_numbers_prior_nll"],
                "recruitment_prior_nll": r["recruitment_prior_nll"],
            }

s21 = read_summary("examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_fit_summary.csv")
f21 = fleet_summary("examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_age_comp_residual_diagnostics.txt")
p21 = prior_blocks("examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_parameter_sanity_diagnostics.txt")

print("level,objective,grad_norm,longline_mean_abs,longline_max_abs,purse_seine_mean_abs,purse_seine_max_abs,m_prior,ll_sel_prior,init_prior,rec_prior")
print(f"level20_sigma_1.75,{s20.get('objective','')},{s20.get('grad_norm','')},{f20.get('longline',{}).get('mean_abs','')},{f20.get('longline',{}).get('max_abs','')},{f20.get('purse_seine',{}).get('mean_abs','')},{f20.get('purse_seine',{}).get('max_abs','')},,{p20.get('longline_selectivity_prior_nll','')},{p20.get('initial_numbers_prior_nll','')},{p20.get('recruitment_prior_nll','')}")
print(f"level21_age_m,{s21.get('objective','')},{s21.get('grad_norm','')},{f21.get('longline',{}).get('mean_abs','')},{f21.get('longline',{}).get('max_abs','')},{f21.get('purse_seine',{}).get('mean_abs','')},{f21.get('purse_seine',{}).get('max_abs','')},{p21.get('age_based_m_prior_nll','')},{p21.get('longline_selectivity_prior_nll','')},{p21.get('initial_numbers_prior_nll','')},{p21.get('recruitment_prior_nll','')}")
PY

echo
echo "== Level 21 M diagnostic consistency =="
grep -R "m_at_age_from_level21_par\|m_at_age" -n   examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/diagnostics   | head -80
