# Swordfish input schema and seasonal dynamics

## Scope and source separation

The parser reads original observations and control settings from the supplied
`swo2023_v004.dat` and `swo2023_v007.ctl`. It does not obtain model observations
from fitted report output. The parser rejects unsupported structural switches,
missing columns, nonfinite input, duplicate observations, and unread trailing
rows. This is deliberately a reader for this assessment's layout, not a general
SS3 importer.

`prepare_model.py` writes the original-input schema to `build/swordfish_inputs.json`
and generates a typed C++ data factory in `build/swordfish_input.hpp`. Recompile
after changing the inputs. Generated files are excluded from git; their generator
and source hashes provide reproducibility. The schema retains source line numbers
and original comments for parameter specifications, but derives canonical parameter
names from their structural positions. Numeric fleet definitions override stale
fleet names in comments.

The same preparation command separately creates `build/ss3_replay.txt` from the
supplied report. That file is strictly a validation fixture, not a model dataset.
`build/preparation_summary.json` records fixed-parameter differences and counts.

## Input contents

| Item | Representation and interpretation |
|---|---|
| Dimensions | 1975–2021; 4 seasons of 3 months; 2 subseasons; 2 sexes; ages 0–15+; 1 population area |
| Catch | 2,223 input rows including equilibrium catches (`year=-999`); preserve fleet-specific biomass/numbers codes |
| Abundance observations | 164 rows; preserve observation month, numeric fleet, observation SE, units and error codes |
| Length observations | 446 samples of 102 input entries (51 bins per sex); 351 included and 95 excluded by negative year |
| Length controls | Per-fleet tail compression, added constants, sex combination, likelihood family and minimum sample size |
| Biological controls | Sex-specific M, von Bertalanffy growth, weight/length, female maturity, sex allocation, first mature age |
| Recruitment | Beverton–Holt parameters, early/main estimation periods, phases, bounds, bias ramp and settlement assignment |
| Selectivity | Length patterns 0, 1, 15 and 24; age patterns 0; fleet mirroring and explicit time-block parameters |
| Parameter specifications | 74 input parameters, retaining 7/14-column bounds, initial values, priors, phases and links |
| Likelihood controls | Q setup, variance adjustments, phase-specific likelihood multipliers and SD-offset setting |

Signed years and original sample sizes are preserved. Do not interpret the two
sex blocks of a combined-sex sample as independent observations. Several input
blocks are duplicates. Conversion to probabilities and SS3
composition processing are implemented for the supplied combined-sex, negative
minimum-tail, multinomial configuration. Other sex/partition or compression
configurations fail explicitly rather than silently using this likelihood.

Survey fleet 25 has index units code 33 (recruitment), not vulnerable biomass.
Survey fleets 24 and 27 have zero likelihood multipliers in the supplied control
file. Keep observations and their use/weight controls separate.

SS3 parameter column 3 is the initial value; it is not necessarily the value
used by the supplied fit, which initialized from `ss.par`. In particular, fixed
male `L_at_Amin_Mal_GP_1` is 83.2 in the control and 83.4 in the report. The input
schema retains 83.2; the report-growth audit explicitly uses 83.4. No silent
replacement occurs.

## State, units and transitions

`State<T>` stores `[sex, age]` in a 32-cell array, with females first. Population
numbers and number-based catches are in thousands of fish; body weights are kg;
their product is metric tonnes. M and fleet-specific F are annual instantaneous
rates. A season has duration `d=0.25` years.

For each sex and age, given beginning-season numbers N:

```
Z = M + sum_f F_f
survivors = N * exp(-Z*d)
catch_f = N * F_f * (1-exp(-Z*d))/Z
natural_deaths = N * M * (1-exp(-Z*d))/Z
```

The integrated-exposure function has a Taylor limit at zero mortality, preserving
both values and derivatives. Survivors plus all mortality removals equal N.
There is no movement in this one-area model. All ages advance on January 1;
oldest-age survivors accumulate with entrants to age 15+. Settlement occurs at
the start of season 3 at age 0, with the supplied female fraction of 0.5. The
actual control assignment is `1 7 1 0` (growth pattern, month, area, age).

The kernel includes Beverton–Holt recruitment:

```
R(B) = 4*h*R0*B / (B0*(1-h) + B*(5*h-1))
```

It is checked at B=0, B=0.2*B0 and B=B0. The historical replay uses report
recruitments rather than this function because independent stock/recruit,
initialization and bias-adjustment matching are handled in a separate, fully
parameter-driven population audit described below.

Mean lengths for ordinary ages use the two-reference-age von Bertalanffy curve,
with linear growth from minimum length to the first reference age. The birth
fraction is 0.5 years. `biology.hpp` implements the aggregate plus-group growth
adjustment, length distributions, maturity and weight integration, and fleet
selectivity. The earlier ordinary-age `mean_length()` helper remains available
for the original audit; the integrated biology uses `length_at_age()` for all
16 ages, including age 15+.

## Explicit spawning-timing compatibility

This comparison uncovered a material SS3-version convention. The report gives
spawn month 5, season 2, and within-season spawn fraction 0.333333. Calendar time
from the start of season 2 is 1/12 year. However, the supplied output is reproduced
by multiplying the annual Z rate directly by 1/3, the fraction of the season.

The pinned [SS3 v3.30.20 population dynamics source](https://github.com/nmfs-stock-synthesis/stock-synthesis/blob/v3.30.20/SS_popdyn.tpl)
uses `spawn_time_seas` directly in the spawning survival exponent, whereas
ordinary seasonal survival multiplies Z by `seasdur(s)`. The numerical replay
confirms the supplied report's behavior across all 47 historical spawning events.
The code therefore exposes two separately named calculations:

- `spawning_biomass`: elapsed calendar years times annual Z.
- `spawning_biomass_ss33020`: elapsed fraction of the spawning season times
  annual Z, for reproducing the supplied SS3 output.

Both sum female numbers times survival times mature weight. These functions
require mature weight at age; the current conditional audit calculates that
vector from the biological parameters and length distributions. `spawning_timing_comparison.csv` records both results using the
same states, rates and mature weights. This is a conditional timing diagnostic,
not a sensitivity refit or a new stock-status estimate.

## What the audit establishes

Starting from the single report state for 1975 season 1, the C++ kernel propagates
all 47 years without resetting to subsequent report states. It uses the report's
fleet-specific apical F and annual recruitment as external drivers. It calculates
F-at-age from apical F and the new selectivity/length integration, then compares:

- 6,016 beginning-season and 6,016 midseason sex/age states;
- 114,304 fleet/season/sex/age catches;
- 188 total and age-1+ biomass values, using independently calculated body weights;
- 47 spawning biomass values using explicit SS3 timing and calculated mature weights;
- 120 beginning-season and 120 midseason mean lengths, excluding aggregate
  plus-group growth;
- mortality mass balance for all 6,016 season/sex/age cells.

The audit allows for the report's limited printed precision. Counts and numbers
use mixed absolute/relative tolerances; the CSV records maximum absolute error
and maximum fraction of tolerance. It fails on any tolerance fraction above 1,
on missing fixture values, or on unread trailing values. Tests independently
exercise zero mortality, plus-group aging, sex allocation, recruitment curve
identities, both spawning conventions, and Quadra automatic derivatives against
analytic derivatives (including the zero-mortality limit).

These checks establish the seasonal kernel and input representation. The
likelihood and fit checks below add objective matching and joint estimation.
Uncertainty and independently calculated reference points remain outstanding.

## Length integration and selectivity

The implementation follows the matching public SS3 v3.30.20 sources:
[biology and growth](https://github.com/nmfs-stock-synthesis/stock-synthesis/blob/v3.30.20/SS_biofxn.tpl),
[age-length key](https://github.com/nmfs-stock-synthesis/stock-synthesis/blob/v3.30.20/SS_ALK.tpl),
[selectivity](https://github.com/nmfs-stock-synthesis/stock-synthesis/blob/v3.30.20/SS_selex.tpl),
and [time/bin setup](https://github.com/nmfs-stock-synthesis/stock-synthesis/blob/v3.30.20/SS_readdata_330.tpl).

- Length-at-age has a normal distribution. Its SD is mean length times CV;
  CV is interpolated using length between the two reference ages. Both CVs
  equal 0.1 in these inputs. Tail probability is accumulated into the first
  and last bins. The CDF uses a scalar-generic approximation with maximum
  absolute error below 7.5e-8; AD differentiates the same approximation.
- Length-bin calculations use midpoints 12.5 through 262.5 cm. Mean weight is
  the probability-weighted sum of `a*length^b`, not weight at mean length.
- For the supplied growth option -999, January plus-group length averages
  linear increments from exact age-15 length toward Linf, with weights
  `exp(-0.2*j)` for j=0 through 15. Subsequent seasonal lengths grow from this
  adjusted January value. This reproduces the reported plus-group sizes.
- Female maturity is the length-based logistic curve. Mature weight integrates
  maturity times fecundity over length, with zero fecundity below age 2.
  SS3 uses the start of the subseason containing the spawn month for this
  integration: quarter 2, subseason 1 (year fraction 0.25) for these inputs.
  This is separate from its spawning-survival exponent convention.
- Length selectivity implements pattern 0 (constant), pattern 1 (logistic with
  95%-width), pattern 24 (double normal with smooth plateau joins), and pattern
  15 (mirror). Pattern 24 currently requires the supplied -999 endpoint settings;
  alternative endpoint/cut-bin settings fail explicitly. No retention or
  sex-specific selectivity offsets are active in the supplied inputs.
- Time-block parameters replace their base values in the supplied model.
  Both ends of each year interval are inclusive. Mirrors use the resulting
  year-specific curve, so a source-fleet block also affects its mirrors.
- Age selectivity, selected biomass, mean selected length and selected weight
  integrate the length-specific curves against length-at-age probabilities.
  No normalization by maximum age selectivity is introduced.

`prepare_biology.py` generates a separate fitted parameter set and records its
mapping in `biology_parameter_mapping.csv`. It checks structural order, labels,
fleet identifiers, parameter bounds, and total count before mapping. Original
control initial values remain unchanged. The comparison uses the report's
printed fitted values; predictions do not read its ALK, selectivity curves,
weights or mature weights. Those appear only as validation targets.

The biology audit performs 194,290 comparisons, including all 13,056 reported
ALK cells, 9,588 length-selectivity values, 162,432 age-selectivity values,
1,598 effective parameters over historical years, and end-year biological
quantities for every age, sex and season. Duplicate report snapshots are retained
as separate comparisons. Mixed absolute/relative tolerances account for printed
parameter and output precision; the largest tolerance fraction must be <=1.
`biology_audit.csv` identifies the worst case for every quantity. Additional tests
cover probability mass and tails, minimum mature age, inclusive block boundaries,
mirrors, constant selection, plus-group growth, and AD derivatives of integrated
predictions with respect to selectivity and growth parameters.

The historical replay now calculates weights, mature weights and F-at-age using
this module. Initial population, annual recruitment and fleet apical F remain report drivers
for this older conditional replay. The separate population audit below computes
those quantities independently.

## Independent initialization, recruitment and hybrid catches

`population.hpp` implements a full 1975-2021 forward population calculation.
The parameter-driven `run_population()` function accepts the original catch
observations, model controls, biological/stock-recruit/selectivity parameters,
initial-F parameter and annual recruitment deviations. It does not read a report,
report-derived population state, annual recruitment estimate or fishing rate.

The population audit supplies the 74 mapped fitted parameters and 57 fitted
recruitment deviations from the report. These are comparison parameter values,
not a new Quadra estimate. `prepare_population.py` maps `Early_InitAge_a` to year
`1975-a`, and early/main recruitment-deviation labels to their named years. It
requires exactly one value for each year 1965-2021 and excludes the forecast
deviation. Original input initial values remain unchanged.

### Initialization

Season-specific biological selection and natural mortality determine annual
survival for each sex and age. For recruitment R settling at the beginning of
quarter 3, January age-1 abundance is its sex allocation times survival through
quarters 3 and 4 at age 0. Each subsequent age equals the previous age times its
full-year survival; age 15+ divides incoming survivors by `1-exp(-annual_Z_15)`.
This closes the infinite tail analytically. The model then advances the
seasonal equilibrium to calculate spawning biomass with the explicit SS3 timing.

First it calculates the unfished equilibrium at `R0=exp(log_R0)` and obtains B0.
Next it applies the initial-F parameter for fleet 6 in season 1 and calculates
spawning biomass per recruit, phiF. With the supplied Beverton-Holt setup,

```
R_initial = (4*h*R0 - (1-h)*B0/phiF) / (5*h-1)
```

The initial equilibrium uses this recruitment. Beginning-1975 ages 1-10 are
then multiplied by their cohort's early recruitment deviation and bias
adjustment. Older ages are outside the early deviation period and retain their
equilibrium abundance. Nonpositive equilibria are rejected explicitly; SS3's
additional near-crash objective penalties are not yet implemented.

### Recruitment and bias

The bias ramp is evaluated from the original control values: zero through
1967.7752, increasing linearly to the maximum 0.9366 at 1995.2499, full through
2020.8972, then declining to zero at 2021.7304. This is evaluated at each integer
year. Recruitment is generated from calculated spawning biomass using
Beverton-Holt and multiplied by `exp(deviation - bias*sigmaR^2/2)` (plus the
supplied regime term, which is zero). The new recruits enter at the beginning
of quarter 3. Fitted deviations are parameters; no reported recruitment outputs
are used as drivers in this calculation.

### Catch conditioning

The solver follows the pinned [SS3 v3.30.20 population dynamics source](https://github.com/nmfs-stock-synthesis/stock-synthesis/blob/v3.30.20/SS_popdyn.tpl)
and [recruitment source](https://github.com/nmfs-stock-synthesis/stock-synthesis/blob/v3.30.20/SS_recruit.tpl).
It preserves the supplied hybrid method 3 and five-pass setting:

1. Calculate a Pope starting estimate using midseason natural survival and
   fleet-specific vulnerability, with a smooth harvest cap near 0.95.
2. Perform four tuning passes. Calculate total Z and integrated exposure, an
   aggregate catch adjustment, and fleet-specific updated F, retaining SS3's
   0.0001 stabilizers and smooth approach to maximum F=4.
3. Recalculate Z and catches from the final apical rates. Those rates determine
   seasonal survival and the following population state.

Each fleet uses its original numbers or biomass catch units. The aggregate
adjustment follows SS3's sum across these mixed units; it is a numerical tuning
heuristic, not a physical total catch. A zero observed catch gives zero F. An
unattainable catch remains visibly unmatched at a finite capped F; the solver
never creates fish to satisfy a target. The forward output records each season's
largest relative catch mismatch.

### Population validation and derivative checks

`make audit-population` computes the entire trajectory before opening its
report-derived validation targets. It compares virgin and initial equilibrium
recruitment, SSB and numbers; beginning-1975 numbers after cohort adjustments;
all historical population states, recruitment, apical F, F-at-age, catch-at-age,
biomass and spawning biomass. It writes `population_audit.csv` and
`population_trajectory.csv`. This is conditional parameter-value reproduction,
not parameter estimation or objective-function equivalence.

Tests cover mixed number/biomass catches generated from known fishing rates,
zero catches, an empty population, negative catches, unattainable targets,
recruitment ramp intervals, and automatic derivatives through initialization
and the finite hybrid algorithm. Derivatives with respect to abundance, initial
F and an early deviation are compared with centered numerical checks. Fixed
biological quantities are cached for this test; their growth/selectivity
derivatives are tested separately. Scalar arrays accumulated by these routines
are explicitly initialized with numerical zeros, including under Quadra AD.

## Likelihood and the first joint fit

The objective follows the pinned [SS3 objective source](https://github.com/nmfs-stock-synthesis/stock-synthesis/blob/v3.30.20/SS_objfunc.tpl),
[observation predictions](https://github.com/nmfs-stock-synthesis/stock-synthesis/blob/v3.30.20/SS_expval.tpl),
[preprocessing](https://github.com/nmfs-stock-synthesis/stock-synthesis/blob/v3.30.20/SS_prelim.tpl),
and [multinomial helper](https://github.com/nmfs-stock-synthesis/stock-synthesis/blob/v3.30.20/SS_miscfxn.tpl).

CPUE predictions integrate selected numbers or biomass at each survey's timing.
The recruitment index uses SS3's recruitment quantity, including its ordering
relative to the year's spawning calculation. Float Q is the inverse-variance
weighted mean of `log(observed)-log(base+1e-6)`. The lognormal negative
log-likelihood includes half the squared standardized log residual and
`sd_offset*log(SE)`, omitting the same constant as SS3. Survey 24's input SE
receives the supplied additive 0.11 adjustment. Zero-lambda surveys retain
predictions and profiled Q but contribute zero to the objective.

Length predictions use season-integrated abundance exposure for fishing fleets
and midseason length distributions. Sexes are combined in predictions; only the
first raw sex block represents the supplied combined-sex observation. Both
observed and predicted proportions receive 0.001 per bin and are normalized
again. The negative tail threshold means no bins are compressed in this model.
The sample-size multiplier defaults to one and the minimum is 0.1. For each
included sample the objective is `N*sum(obs*log(obs/expected))`; the perfect-fit
offset is retained. All 351 included samples have the supplied phase-persistent
length weight of 0.5; the 95 negative-year samples are excluded.

Recruitment contributes `sum(dev^2/(2*sigmaR^2)) + sum(bias)*log(sigmaR)`, with
independent deviations and the supplied SD-offset setting. Catch contributes
SS3's stabilized log-catch residual. The sole positive equilibrium-catch datum
has zero lambda. Ordinary parameter priors are all disabled; selectivity soft
bounds reproduce SS3's symmetric-beta penalty with its 0.0001 stabilizer.
Forecast recruitment, parameter deviations, initial-regime and crash penalties
are zero in the validated configuration. Nonzero ordinary priors, equilibrium
catch likelihood, autocorrelation and active projections are outside this
implementation's current supported configuration and fail explicitly where
applicable. This is not a general SS3 likelihood importer.

`audit-likelihood` computes the population and likelihood before opening the
report targets. It compares 12 printed objective rows, all 164 index fits
(base, Q, predictions, SE, residual and likelihood) and 351 length fits (sample
size, observed/predicted mean and likelihood). At the supplied printed estimates:

| Component | Quadra value | Printed SS3 value |
|---|---:|---:|
| Survey | -164.951027175 | -164.951 |
| Length composition | 481.252427209 | 481.252 |
| Recruitment | -27.018039927 | -27.018 |
| Selectivity soft bounds | 0.016279237 | 0.0162791 |
| Catch | 2.923370e-8 | 2.92335e-8 |
| Total | 289.299639373 | 289.300 |

The model wrapper has 28 estimated control parameters, 20 early deviations and
36 independent main deviations. The last main deviation is minus the sum of
the preceding 36, reproducing SS3's bounded deviation-vector zero-sum constraint.
That gives **84 independent parameters**. Native parameter bounds are represented
with scaled logistic transforms; no change-of-variable Jacobian is added to this
maximum-likelihood target. The dependent deviation is also checked against its
bounds. All deviations are optimized jointly, not integrated as random effects.

`fit-swordfish` uses Quadra's public L-BFGS interface. The first run, starting
from the SS3 report estimates, converged at objective **289.299639305825** with
an L2 gradient norm of **8.9413e-5** in fitting coordinates, below the requested
1e-4. The initial norm was 0.036050. SS3's printed gradient uses its own
parameterization/convention and is not compared numerically to this norm.
This is a common-start reproduction fit; it does not establish global optimality
or robustness from distant starting values. Every full-objective automatic
derivative was compared to a centered numerical derivative; the largest
absolute discrepancy across all 84 was 1.11e-6 (limit 1e-5).

The model caches biology by selectivity-block combination and length keys by
subseason within each evaluation. No AD vertices are retained across evaluations.
Quadra's gradient-only helper now uses its first-order backward pass, and its
public optimizer skips the latent-mode and Hessian/logdet calculations when
there are no random effects. These framework changes preserve the joint target
and have a dedicated fixed-only optimizer regression plus existing scalar,
transformation, stats and random-effect optimizer checks.

Fitted CSVs contain annual and seasonal biomass, SSB, recruitment, fleet F,
annual F, catches, data fits and objective terms. Annual F follows the supplied
**SS3 units 4**: `log(sum(N_Jan*exp(-M))/sum(N_end))` pooling sexes and ages 1-10.
It is not the units-5 arithmetic mean across ages; descriptive prose in the
reference documents should not replace this version-specific definition.

## Remaining implementation

1. Fit robustness from broader starting values and uncertainty estimation.
2. Reference points, dynamic unfished baselines and projection scenarios.
3. Any Laplace-marginalized alternative, reported as a distinct inference target.
4. Striped marlin as a separate model when its inputs are supplied.
