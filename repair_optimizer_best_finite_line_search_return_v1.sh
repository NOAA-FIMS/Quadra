#!/usr/bin/env bash
set -euo pipefail

echo "== Repair optimizer: return best finite iterate on LBFGS line-search failure =="

stamp="$(date +%Y%m%d_%H%M%S)"
mkdir -p .quadra_patch_backups

opt="core/optimizer.hpp"
[[ -f "$opt" ]] || { echo "ERROR: missing $opt" >&2; exit 1; }

cp "$opt" ".quadra_patch_backups/optimizer.hpp.best_finite_line_search_return_${stamp}.bak"

python3 - <<'PY'
from pathlib import Path
import re

p = Path("core/optimizer.hpp")
s = p.read_text()

if "QUADRA_RETURN_BEST_FINITE_ON_LINE_SEARCH_FAILURE_V1" in s:
    print("Best-finite line-search return patch already installed.")
    raise SystemExit(0)

# First try the clean path: insert after creation of msg in a runtime_error catch.
msg_patterns = [
    r'(const\s+std::string\s+msg\s*=\s*e\.what\s*\(\s*\)\s*;\s*)',
    r'(std::string\s+msg\s*=\s*e\.what\s*\(\s*\)\s*;\s*)',
]

insert = r'''
    // QUADRA_RETURN_BEST_FINITE_ON_LINE_SEARCH_FAILURE_V1
    if ((msg.find("line search") != std::string::npos ||
         msg.find("sufficiently decrease") != std::string::npos) &&
        best_finite_seen) {
      OptResult out;
      out.par = eigen_to_std_vector(best_finite_x);
      out.value = best_finite_fx;
      out.grad_norm = best_finite_grad_norm;
      out.converged = best_finite_grad_norm < opts.fixed_grad_tol;
      out.iterations = iterations;
      out.message =
          out.converged
              ? "accepted best finite iterate after LBFGS line-search failure"
              : "stopped before requested fixed-effect gradient tolerance";
      out.u_hat = best_finite_u_hat;
      return out;
    }

'''

for pat in msg_patterns:
    m = re.search(pat, s)
    if m:
        s = s[:m.end()] + insert + s[m.end():]
        p.write_text(s)
        print("Inserted best-finite return after runtime_error message creation.")
        raise SystemExit(0)

raise SystemExit("ERROR: could not find runtime_error message creation in core/optimizer.hpp")
PY

cat > inspect_optimizer_best_finite_line_search_return_v1.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

echo "== Optimizer best-finite line-search return markers =="
grep -n "QUADRA_RETURN_BEST_FINITE_ON_LINE_SEARCH_FAILURE_V1\\|best finite iterate after LBFGS line-search failure\\|sufficiently decrease" core/optimizer.hpp

echo
echo "== Relevant best_finite symbols =="
grep -n "best_finite_seen\\|best_finite_x\\|best_finite_fx\\|best_finite_grad_norm\\|best_finite_u_hat\\|eigen_to_std_vector\\|fixed_grad_tol" core/optimizer.hpp | head -80

echo
echo "Now rerun:"
echo "  ./inspect_opakapaka_level1_reporting_v6.sh"
echo "or:"
echo "  ./inspect_opakapaka_level1_reporting_v7.sh"
SH
chmod +x inspect_optimizer_best_finite_line_search_return_v1.sh

echo
echo "Backups saved with suffix: best_finite_line_search_return_${stamp}.bak"
echo "Run:"
echo "  ./inspect_optimizer_best_finite_line_search_return_v1.sh"
echo "  ./inspect_opakapaka_level1_reporting_v7.sh"
