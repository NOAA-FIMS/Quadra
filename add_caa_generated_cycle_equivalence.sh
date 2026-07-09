#!/usr/bin/env bash
set -euo pipefail
echo "Bootstrap installer for generated cycle equivalence test."
mkdir -p examples/NMFS/pifsc_bigeye_tuna/v2/19_generated_cycle_equivalence_caa
cat > run_bigeye_v2_19_generated_cycle_equivalence_caa_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail
./generate_bigeye_v2_caa_assessment_cycle_from_ir.sh
echo "Placeholder: compile and run generated-cycle equivalence test."
SH
chmod +x run_bigeye_v2_19_generated_cycle_equivalence_caa_check.sh
echo "Installed."
