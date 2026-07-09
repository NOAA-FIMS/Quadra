#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"

mkdir -p "$BASE/20_generated_optimizer_recovery_caa"

SRC="$BASE/level09_optimizer_check/bigeye_v2_level09_optimizer_check.cpp"
DST="$BASE/20_generated_optimizer_recovery_caa/bigeye_v2_20_generated_optimizer_recovery_caa_check.cpp"

cp "$SRC" "$DST"

python3 - <<'PY'
from pathlib import Path
p = Path("examples/NMFS/pifsc_bigeye_tuna/v2/20_generated_optimizer_recovery_caa/bigeye_v2_20_generated_optimizer_recovery_caa_check.cpp")
s = p.read_text()

s = s.replace(
    '#include "../architecture/assessment/assessment_cycle.hpp"',
    '#include "../architecture/assessment/generated_assessment_cycle.hpp"'
)

s = s.replace("AssessmentCycle{}", "GeneratedAssessmentCycleFromIR{}")
s = s.replace("AssessmentCycle cycle", "GeneratedAssessmentCycleFromIR cycle")
s = s.replace("PASSED: Bigeye v2 Level09 optimizer recovery regression",
              "PASSED: Bigeye v2 CAA generated optimizer recovery regression")

p.write_text(s)
PY

cat > run_bigeye_v2_20_generated_optimizer_recovery_caa_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

./generate_bigeye_v2_caa_assessment_cycle_from_ir.sh

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/20_generated_optimizer_recovery_caa/bigeye_v2_20_generated_optimizer_recovery_caa_check.cpp \
  -o build/examples/bigeye_v2_20_generated_optimizer_recovery_caa_check

./build/examples/bigeye_v2_20_generated_optimizer_recovery_caa_check
SH

chmod +x run_bigeye_v2_20_generated_optimizer_recovery_caa_check.sh

echo "created generated optimizer recovery check"
