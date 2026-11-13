# Swordfish operational release contract

The supported release entry point is `python3 scripts/run_release.py`. It builds
and tests an isolated copy of source and all eight supplied inputs, fits the
model, recalculates covariance, checks perturbed starts and projections, generates
figures, and publishes a new result directory only after all acceptance gates
pass. It never updates or overwrites an existing release directory.

```
python3 scripts/run_release.py --output runs/swordfish-release-001
```

Each stage has a 30-minute timeout; timeout/interruption terminates its process
group and retains failure evidence. The supported platforms are POSIX macOS/Linux.

The command requires Python 3.9+, a C++17 compiler (`CXX`, default `c++`), make,
and base R on PATH. The framework's Eigen and LBFGSpp dependencies are bundled.
Use a simple compiler executable in CXX, not a shell command containing flags.
Original inputs must match `reference/source_manifest.json`; they are kept local
and are not uploaded by the runner or committed into CI.

The result contains the copied source, inputs, build outputs, stage logs and
`release_manifest.json`. Manifest status `passed` means the automated numerical
and operational acceptance gates passed for that exact source/input snapshot.
It does not grant independent scientific approval of a stock assessment.

## Publication and failure behavior

A sibling directory lock prevents two cooperative runs from using the same
output name. Work occurs in a uniquely named `.inprogress-...` directory on the
same filesystem as the destination. On success a directory rename publishes the
complete bundle. On failure a `.failed-...` bundle retains logs and `failure.json`,
and the destination is not published. Failed fits do not replace existing fit
CSV outputs. A final manifest records the environment, acceptance gates and stage log paths,
and binds copied source, inputs and all build artifacts by SHA-256. Effective
compiler flags are printed in the retained compilation logs.

Publication protects against process/stage failure and concurrent use of this
runner. It is not a transaction across distributed filesystems or a guarantee
against hostile filesystem modifications, power loss, or a second program that
ignores the destination lock. Stale locks after an uncatchable termination require
operator inspection before removal; no automatic lock stealing is performed.

The low-level `make` targets still use the working directory's `build/` folder.
They are developer tools and may overwrite scratch outputs. Use the isolated
runner for unattended work and release publication. An individual scratch output
is not a completed release until its manifest and acceptance gates exist.

## Covariance identity

`make uncertainty` writes a pending content identity before computation. It binds
all original files, scenario configuration, case source/tests/scripts, framework
headers/vendor files and the complete native fitted parameter CSV. Finishing
requires the same identity and stamps the covariance, coordinate mapping, both
Hessians and curvature diagnostics. Robustness and projection-uncertainty targets
verify this stamp before reuse; mismatches fail with instructions to recalculate.

A model-code change can invalidate covariance even if fitted coordinate values
and file timestamps are unchanged. Source comments currently count as source
changes, deliberately favoring a conservative rebuild over silent reuse.

## Automated gates

`make test` includes biological, population and likelihood SS3 audits, the full
84-coordinate gradient check, baseline dynamics/likelihood tests, equilibrium
benchmark tests, projection stationarity tests, delta-method checks, strict CSV
failure tests and Python input/runtime failure tests. The isolated release adds:

- Five shared-framework regression programs, including the fixed-only optimizer
  path and an existing Laplace optimization case.
- A converged fit and annual SS3 comparison, plus ten management benchmark checks.
- Two-step curvature and selected standard-error comparison, with no silent ridge.
- Three local perturbed starts converging to the same objective.
- Dynamic unfished and historical SPR reproduction, five projection scenarios,
  full correlated projection uncertainty and numerical derivative checks.
- PDF/PNG scientific figures. New document PDFs are separately rendered/reviewed
  before delivery; automated figure generation alone is not a visual sign-off.

The GitHub workflow `swordfish-release.yml` defines Ubuntu/macOS portable and
framework tests plus AddressSanitizer/UndefinedBehaviorSanitizer checks. It does
not include the private raw data. Hosted CI completion cannot be asserted until
the workflow is pushed and executed. The local release evidence records the
actual platform and compiler tested, not an assumed cross-platform certification.

## Validation boundary and remaining scientific work

The implementation is intentionally restricted to the supplied SS3 v3.30.20
swordfish configuration: 1975-2021 history, four quarters, two subseasons, two
sexes, ages 0-15+, one area, 19 catch fleets and eight indices. Projection readers
validate exactly five supported scenario rules and the documented reference
windows. This is a defined support contract, not a general SS3 interpreter.

Scientific limitations persist: published bootstrap forecast configurations and
100 bootstrap model fits are absent; retrospective analysis, broader starts and
structural sensitivities remain future work. Projection uncertainty is a local
parameter-covariance approximation excluding future recruitment variability.
The original report's scenario-description conflicts are retained explicitly.
An independent fisheries-model review is required before operational assessment
advice. Striped marlin requires its own inputs and separate validation.

## Documentation coverage

The guide and API catalog cover all explicit named declarations in the case's
`include`, `src`, `scripts` and `tests`, plus R plotting functions and a focused
appendix for changed Quadra functions. Anonymous lambdas, third-party libraries,
and generated factory bodies are excluded from coverage. Factory generation and
versioned fixture formats are documented in the guide. Rebuild the inventory with
`python3 scripts/document_api.py`; it requires Clang and the generated data
headers. Missing descriptions fail generation rather than creating placeholders.

PDF regeneration uses `docs/requirements.txt` in an isolated environment:

```
python3 -m venv /tmp/quadra-docs
/tmp/quadra-docs/bin/pip install -r docs/requirements.txt
/tmp/quadra-docs/bin/python scripts/build_guide.py
```
