# Bigeye Tuna Diagnostic-Guided Scientific Reasoning Log

This document records the scientific reasoning behind each model-building step.

The goal is not to build a large model quickly. The goal is to build a model
whose complexity is justified by diagnostics.

## Working philosophy

A model feature should earn its place.

For each new feature, we ask:

1. What biological or observation-process hypothesis motivates the feature?
2. What data should identify the feature?
3. What diagnostics changed after adding the feature?
4. Did the feature improve interpretation, stability, or information content?
5. Did the feature introduce confounding, weak curvature, or overparameterization?
6. Should the feature be retained, simplified, pooled, regularized, or removed?

## Diagnostic interpretation rules

### Random-effect Hessian health

If the random-effect Hessian is positive definite and well-conditioned, but the
fixed-effect optimizer stalls, the first suspect is not the Laplace machinery.
The first suspect is fixed-effect scaling, data conflict, or weak
identifiability.

### Fixed-effect gradient health

If the maximum fixed-effect gradient is concentrated on one parameter, that
parameter becomes the first diagnostic target.

Example:

```text
max_gradient_parameter: log_q
```

This suggests catchability/abundance scaling may be weakly identified or
conflicted with other fixed effects.

### Scale confounding

If recruitment scale increases while catchability decreases, the model may be
trading abundance scale against observation scale.

Example pattern:

```text
R0 increases
q decreases
index fit remains similar
```

Interpretation:

The data may not separate absolute abundance from catchability.

### Complexity rule

Do not add movement, regions, tagging, or fleet-specific selectivity when the
current simpler model already shows unresolved q/F/R0 confounding.

## Level 0: single-region baseline

### Model change

The Level 0 model starts with:

- one region
- one aggregate catch series
- one abundance index
- age-structured dynamics
- recruitment deviations as random effects
- fixed effects for recruitment scale, fishing mortality, catchability, and
  selectivity

### Scientific purpose

Establish a baseline diagnostic package before adding fleet structure.

### Diagnostic expectation

Level 0 should answer whether the base model is numerically healthy and whether
the random effects are identifiable under a simple observation process.

### Decision

Level 0 is the baseline. It should not be treated as a tuna assessment. It is a
diagnostic control case.

## Level 1A: aggregated multi-fleet scaffold

### Model change

Two synthetic fleets were introduced in the data file:

- longline
- purse seine

The first scaffold aggregated fleet catches by year and averaged fleet indices
by year so the existing Level 0 biological model could run unchanged.

### Scientific hypothesis

If fleet information is aggregated into one catch series and one index, the
model should behave similarly to Level 0 unless the aggregation creates a scale
conflict.

### Diagnostics observed

The first Level 1A scaffold produced:

```text
objective:                 115.031861
gradient_norm:             0.429520
converged:                 no
max_gradient_parameter:    log_q
max_gradient_value:        0.346158
Huu positive definite:     yes
Huu condition number:      6.20
structural nonzeros:       362 / 400
max random-effect corr:    about 0.21
```

### Scientific interpretation

The random-effect curvature remained healthy. The random-effect Hessian was
positive definite and well-conditioned.

The fixed-effect optimizer stalled, with the largest remaining gradient on
`log_q`.

The aggregated data created a likely scale conflict:

```text
catch = longline catch + purse seine catch
index = average(longline index, purse seine index)
```

This roughly doubled removals while leaving the index scale comparable to a
single-fleet index.

The model responded by increasing the recruitment/abundance scale and reducing
catchability:

```text
R0 increased
q decreased
```

This is consistent with abundance/catchability scale confounding.

### Decision

Level 1A is a useful diagnostic failure.

It should not be treated as a validated multi-fleet model.

The next step is not to add regions, movement, or tagging. The next step is to
make the fleet observation process explicit.

## Level 1B: explicit fleet catchability

### Planned model change

Introduce explicit fleet-level observation terms:

- shared biological population
- fleet-specific catch/index observations
- fleet-specific catchability parameters, initially:
  - `log_q_longline`
  - `log_q_purse_seine`

### Scientific hypothesis

If Level 1A failed because of observation-process aggregation, then allowing
fleet-specific catchability should reduce the `log_q` scaling conflict.

### Diagnostic questions

- Does the model converge?
- Does the maximum fixed-effect gradient decrease?
- Does the maximum gradient remain on a catchability parameter?
- Do fleet-specific q values become identifiable?
- Does Huu remain positive definite?
- Does fixed-effect correlation reveal q/R0/F confounding?
- Does the functional analysis show different influence patterns than Level 1A?

### Decision rule

If explicit fleet q improves convergence and reduces the log_q gradient, then
the problem was likely aggregation-induced observation-process confounding.

If explicit fleet q does not improve convergence, then the issue is more
structural and we should inspect fixed-effect Hessian geometry before adding
additional biological complexity.

## Level 2: fleet-specific selectivity

### Planned model change

Add fleet-specific selectivity only after Level 1B diagnostics show that
fleet-specific catchability is not enough.

### Scientific hypothesis

Different fleets may observe different parts of the age/size distribution.

### Diagnostic questions

- Does fleet-specific selectivity improve fit enough to justify extra
  parameters?
- Which fleet identifies selectivity?
- Are selectivity and q confounded?
- Are selectivity and F confounded?

### Decision rule

Fleet-specific selectivity should be simplified or pooled if diagnostics show
weak curvature or parameter redundancy.

## Level 3: spatial structure

### Planned model change

Add regional abundance states only after fleet observation processes are
diagnostically stable.

### Scientific hypothesis

Regional structure may explain differences among fleets and indices.

### Diagnostic questions

- Are regional abundance states separately identifiable?
- Does q become region-confounded?
- Does the Huu correlation graph show strong regional blocks?

### Decision rule

Do not add movement until regional abundance is at least partially identifiable.

## Level 4: movement

### Planned model change

Add movement among regions.

### Scientific hypothesis

Movement may explain spatial redistribution of abundance.

### Diagnostic questions

- Is movement identifiable separately from regional recruitment?
- Is movement identifiable separately from regional q?
- Does movement reduce residual structure or only absorb noise?

### Decision rule

Movement should be pooled, simplified, or regularized if weak directions show
movement/recruitment/q confounding.

## Level 5: tagging

### Planned model change

Add tagging data only after movement structure has a clear diagnostic need.

### Scientific hypothesis

Tagging may identify movement and reporting rates.

### Diagnostic questions

- Does tagging reduce movement uncertainty?
- Are movement and reporting rate separable?
- Which tag groups contribute information?

### Decision rule

Retain tagging complexity only if it materially improves movement
identifiability.

## Level 2 plan: fleet-specific composition

### Model change

Add synthetic fleet-specific age-composition patterns:

- longline samples older/larger fish
- purse seine samples younger/smaller fish

### Scientific hypothesis

Fleet-specific composition data should provide information about how fleets
sample the population. This may reduce or clarify the R0/F/q confounding seen
in Level 1A and Level 1B.

### Decision rule

Fleet-specific composition earns its place only if diagnostics show improved
separability or interpretable selectivity pressure.

## Level 3 plan: fleet-specific selectivity and vulnerable-biomass indices

### Model change

Add fleet-specific selectivity and connect each fleet index to its own
fleet-vulnerable biomass.

### Scientific hypothesis

Level 2 showed that fleet-specific composition improved optimization but did
not fully resolve the abundance/F/q weak direction. This suggests composition
information was entering the selectivity block but not fully informing the
index likelihood.

If the index is modeled as fleet-vulnerable biomass rather than total biomass,
composition-driven selectivity information can influence the index scale.

### Decision rule

Fleet-specific selectivity and vulnerable-biomass indices earn their place only
if they reduce the primary weak fixed-effect direction or create an interpretable
separation between fleet q, selectivity, F, and R0.

## Level 4 plan: distinct vulnerable-biomass indices

### Model change

Keep the Level 3 model structure but regenerate synthetic fleet indices so that
longline and purse seine observe meaningfully different vulnerable biomass
signals.

### Scientific hypothesis

Level 3 showed that adding fleet-specific selectivity and vulnerable-biomass
indices improved the objective but did not resolve the main weak fixed-effect
geometry. Inspection of the synthetic indices showed that longline and purse
seine index ratios changed only modestly through time, suggesting that the
indices still behaved like a pooled abundance signal.

If the information bottleneck is pooled or nearly pooled index data, then
distinct vulnerable-biomass index trajectories should reduce q/F/R0 confounding.

### Decision rule

If q correlations drop, weak eigenvalues increase, or optimizer behavior
improves, then the diagnostics correctly identified a data-information
bottleneck. If not, the model remains structurally under-identified.

## Wiggle diagnostics

### Purpose

Fixed-effect geometry identifies weak directions through local curvature.
Wiggle diagnostics complement this by perturbing each fixed effect directly and
re-solving random effects through the Laplace machinery.

### Scientific question

Which fixed effects can be moved without meaningfully changing the profiled
objective?

### Interpretation

A small objective change under perturbation suggests that the parameter is
weakly informed or that other model components can compensate for it. A large
change suggests the parameter is locally influential. Asymmetric plus/minus
changes suggest nonlinearity or boundary behavior.

## Level 5 plan: strong selectivity contrast

### Motivation

Level 4 wiggle diagnostics showed that individual q and R0 parameters are
highly influential when perturbed alone, but fleet selectivity parameters can
move with relatively small objective cost. This separates individual parameter
weakness from compensating multi-parameter weakness.

### Model change

Keep the Level 4 model structure but strengthen the synthetic composition
information:

- longline composition dominated by older ages
- purse-seine composition dominated by young ages
- age-composition effective sample size increased

### Scientific hypothesis

If selectivity is weak because composition data are too diffuse or weakly
weighted, then stronger composition contrast should increase selectivity
sensitivity and reduce weak-direction loading on purse-seine selectivity.

### Decision rule

If stronger composition contrast does not improve selectivity geometry, then the
selectivity weakness is structural rather than data-strength related.

## Level 6 plan: purse-seine age-based selectivity

### Motivation

Level 5 successfully reduced the broad R0/F/q confounding, but the remaining
weakest fixed-effect direction was dominated by purse-seine a50. This suggests
the remaining issue is not broad scale identifiability, but an inappropriate
two-parameter logistic selectivity structure for a young-fish purse-seine fleet.

### Model change

Keep longline logistic selectivity, but replace purse-seine logistic selectivity
with a fixed age-based selectivity curve.

### Scientific hypothesis

If the Level 5 weakness is a boundary-like purse-seine a50 problem, then
removing the free purse-seine logistic a50/slope parameters should improve
fixed-effect geometry without requiring more optimizer tuning.

### Decision rule

If the smallest eigenvalue increases, condition number decreases, and wiggle
diagnostics no longer identify purse-seine a50 as weak, then the diagnostics
support using an age-based purse-seine observation process.

## Level 6 result: purse-seine age-based selectivity

### Result

Replacing the purse-seine logistic selectivity curve with a fixed age-based
selectivity pattern improved the model substantially.

The fit moved back toward biologically sensible scale:

- R0 returned near the original recruitment scale.
- q returned near the original catchability scale.
- Fbar dropped from the extreme Level 5 value.
- The maximum fixed-effect gradient became small.

### Diagnostic interpretation

Level 5 successfully reduced broad R0/F/q confounding but exposed a remaining
weak direction dominated by purse-seine a50. Level 6 removed that inappropriate
free logistic purse-seine selectivity parameterization.

This supports the diagnostic-guided conclusion that the purse-seine observation
process in this synthetic example is better represented by age-based
selectivity than by a freely estimated logistic a50/slope curve.

### Diagnostic hardening

Level 6 also exposed a diagnostic-layer issue: some post-fit perturbation
diagnostics can cross regions where the profiled random-effect Hessian is not
factorizable. Safe diagnostics now preserve those failed perturbations as rows
instead of aborting the model run.

## Level 6 follow-up: longline slope geometry scan

### Motivation

Safe wiggle diagnostics showed that decreasing `log_sel_slope_longline` caused a
profiled Laplace/Huu factorization failure while increasing it was stable. This
suggests a one-sided stability boundary around the longline selectivity slope.

### Diagnostic addition

A longline slope geometry scan holds all other fixed effects fixed and sweeps
the longline selectivity slope on the natural scale. For each value, random
effects are re-solved and the profiled objective is evaluated.

### Scientific question

Where does the Laplace approximation become unstable as longline selectivity
slope changes?

### Interpretation

The first failed row below the fitted slope identifies an approximate local
stability boundary. If the boundary is close to the fitted value, the model is
locally sharp and should be treated carefully before adding more structure.

## Level 6 follow-up: recruitment diagnostics

### Motivation

The current recruitment model uses independent annual recruitment deviations
with fixed prior scale `sigma_rec_dev = 0.35`. Recruitment deviations can
represent biological variation, but they can also absorb residual structure from
misspecified selectivity, catch, index, or composition processes.

### Diagnostic addition

Recruitment diagnostics summarize:

- estimated annual recruitment deviations,
- recruitment multipliers,
- prior z-scores,
- prior NLL contribution,
- standard deviation,
- lag-1 correlation,
- first-difference roughness,
- and the year with maximum absolute deviation.

### Scientific question

Are recruitment deviations modest biological variation, or are they carrying
unmodeled structure?

### Interpretation

Large persistent recruitment deviations or strong autocorrelation would suggest
that the independent annual recruitment model is compensating for missing
structure. Modest deviations centered near zero support leaving recruitment as
a simple independent random-effect process for the current diagnostic level.

## Level 7 plan: dual age-based selectivity

### Motivation

Level 6 replaced purse-seine logistic selectivity with fixed age-based
selectivity and substantially improved the model. However, recruitment
diagnostics showed a smooth, persistent recruitment trajectory with high lag-1
correlation. This suggests recruitment may still be absorbing remaining
structural mismatch, possibly from longline logistic selectivity.

### Model change

Level 7 replaces longline logistic selectivity with a fixed age-based
selectivity pattern while retaining purse-seine age-based selectivity. The
fixed-effect dimension is reduced to R0, Fbar, and two fleet-specific q values.

### Scientific hypothesis

If the Level 6 recruitment pattern is compensating for longline logistic
selectivity misspecification, then Level 7 should reduce recruitment
persistence, recruitment prior burden, or age-composition NLL.

### Decision rule

If recruitment lag-1 correlation and rec prior burden drop materially, then
the diagnostic evidence supports the conclusion that recruitment was absorbing
remaining selectivity misspecification. If recruitment persistence remains high,
then the persistent recruitment pattern likely reflects data construction,
observation weighting, or a need for an explicit recruitment process.

## Level 8 plan: tuna weight-at-age

### Motivation

The bigeye tuna examples inherited a small toy/red-snapper-like weight-at-age
schedule that tops out near 5 kg. That is biologically inappropriate for adult
bigeye tuna and affects biomass, vulnerable biomass, catch biomass, q, Fbar,
SSB proxy, and recruitment scaling.

### Model change

Level 8 keeps the Level 6 model structure but replaces the fixed weight-at-age
with a tuna-like schedule:

- age 1: 2 kg
- age 2: 8 kg
- age 3: 18 kg
- age 4: 32 kg
- age 5: 48 kg
- age 6: 65 kg
- age 7: 82 kg
- age 8: 98 kg
- age 9: 112 kg
- age 10: 125 kg

### Scientific question

Were persistent recruitment deviations compensating for unrealistic adult
weight-at-age?

### Interpretation

If recruitment burden, Fbar, q, or objective components shift materially, then
weight-at-age was a hidden structural driver and should be treated as a core
biological input rather than a harmless reporting detail.

## Level 9 plan: estimated natural mortality

### Motivation

Level 8 showed that correcting tuna weight-at-age changed biomass scale but did
not remove the persistent recruitment trajectory. The remaining smooth
recruitment deviations may be compensating for misspecified survival.

### Model change

Level 9 estimates natural mortality as fixed effect `log_m` with a prior
centered at `log(0.18)`. The model otherwise retains the Level 6 structure.

### Scientific question

Can the data support estimated natural mortality, and does allowing M to move
reduce recruitment persistence or age-composition pressure?

### Caution

Estimating M is often weakly identified in stock assessment models. This is a
diagnostic experiment intended to reveal confounding and compensation.

## Level 9 follow-up: objective consistency check

### Motivation

Level 9 produced a large disagreement between the optimizer-reported objective
and the objective-components report. Before interpreting biological results, we
need to verify whether reports are evaluating the same objective and parameter
ordering as the optimizer.

### Diagnostic addition

The objective consistency check constructs the full parameter vector from the
returned fixed effects and random effects, evaluates the objective directly,
and compares that direct joint objective with the profiled Laplace objective.

### Interpretation

The direct joint objective should not equal the profiled Laplace objective,
because the latter includes Laplace adjustment terms. However, it should be on
a plausible scale relative to the reported joint objective components. A large
mismatch indicates stale report code, stale parameter indexing, or inconsistent
fixed assumptions.

## Level 9 follow-up: fixed-effect geometry

### Motivation

After patching Level 9 reports, objective components matched the direct joint
objective. The estimated-M run collapsed recruitment compensation, but it did so
through a biologically implausible ridge involving very large R0, high Fbar,
large q, and nearly flat longline selectivity.

### Diagnostic addition

The Level 9 fixed-effect geometry report computes a finite-difference Hessian
of the profiled Laplace objective across the seven Level 9 fixed effects:

- log_r0
- log_m
- log_fbar
- log_q_longline
- log_q_purse_seine
- logit_sel_a50_longline
- log_sel_slope_longline

### Scientific question

Which parameter combinations define the weak ridge exposed by estimating natural
mortality?

### Interpretation

If `log_m` loads heavily with R0, Fbar, q, or selectivity terms in the weakest
eigen-directions, then estimated M is not independently informed by the data and
is instead part of a compensatory identifiability ridge.

## Level 10 plan: q anchor

Level 9 fixed-effect geometry showed the weakest direction was primarily an
R0-F-q scale ridge. Level 10 fixes longline q at 0.00005 and estimates the
remaining scale terms to test whether one external scale anchor improves the
surface.

## Level 11 plan: fixed M + q anchor

Level 10 fixed longline q but M still escaped high. Level 11 fixes both M and
longline q, leaving R0, Fbar, purse-seine q, and longline selectivity estimated.
This tests whether M was the remaining escape hatch after anchoring q.

## Level 12 plan: tuna life-history correction

Inspection showed that the bigeye levels still carried inherited reef-fish-like
life-history schedules: weight-at-age reached only about 5 kg, and maturity was
near complete by age 5. Level 12 clones Level 11 and keeps fixed M and fixed
longline q, but replaces weight-at-age and maturity-at-age with provisional
bigeye-tuna-like diagnostic schedules. These are placeholders and should be
replaced by stock-specific assessment inputs before production use.

## Level 13 plan: estimated initial numbers

Level 12 converged cleanly after tuna life-history correction, but recruitment
deviations remained strongly positive in the early modeled years. Level 13 tests
whether those recruitment deviations are compensating for an incorrect initial
age structure by estimating initial log-number deviations by age, with a prior
to keep the experiment diagnostic rather than unconstrained.

## Level 14 plan: M sensitivity with estimated initial numbers

Level 13 showed that recruitment deviations were compensating for a bad initial
age structure, especially excessive age-10+ biomass. Level 14 keeps estimated
initial numbers but sweeps fixed M values to test whether higher M reduces
plus-group pressure and the need to crush the age-10+ initial state.

## Level 15 plan: fixed M=0.45 best diagnostic model

Level 14 M sensitivity found the best objective near M=0.45. Level 15 freezes
that setting as the current best diagnostic model and adds age-composition
residual diagnostics. This should identify the remaining fleet/year/age
structure that is driving compensation.

## Level 16 plan: age-specific purse-seine selectivity

Level 15 residual diagnostics showed purse-seine age-1 was greatly overpredicted
every year. Level 16 replaces the fixed purse-seine selectivity vector with
estimated age-specific selectivity logits.

## Level 17 plan: juvenile mortality diagnostic

Level 16 improved objective after estimating purse-seine age selectivity, but
the repeated purse-seine age-1 residual persisted. Level 17 adds a single
juvenile mortality multiplier for ages 1-2, keeping adult M fixed at 0.45.

## Level 16 plan: age-specific purse-seine selectivity

Level 15 residual diagnostics showed purse-seine age-1 was greatly overpredicted
every year. Level 16 replaces the fixed purse-seine selectivity vector with
estimated age-specific selectivity logits.

## Level 18: longline selectivity diagnostic

Level 16 resolved the apparent purse-seine age-composition mismatch once the
diagnostics were corrected to use the objective-consistent prediction path. The
remaining age-composition residuals are concentrated in the longline fleet,
especially older ages. Level 18 clones Level 16 and adds a longline prediction
decomposition report:

- numbers-at-age
- longline selectivity
- selected numbers
- predicted composition
- observed composition
- residual

This level is diagnostic only. Its purpose is to determine whether the remaining
mismatch is caused by the logistic longline selectivity curve being too
restrictive, especially for old ages and the plus group.

## Level 19: flexible longline age selectivity

Level 18 showed that the remaining age-composition mismatch is concentrated in
longline old ages. The logistic longline selectivity curve was essentially fully
selected at both age 8 and age 10, which prevented the model from shifting
selected composition from the plus group back into age 8.

Level 19 replaces the logistic longline selectivity curve with age-specific
longline selectivity logits under an informative template prior. This is the
longline analogue of the successful Level 16 purse-seine age-specific
selectivity diagnostic.

## Level 19 parameter sanity diagnostics

After Level 19 reduced the longline age-composition mismatch, the next question
is whether the improvement came from a defensible selectivity shape or from
excessive tradeoffs with initial numbers. This diagnostic adds block-level prior
penalties and age-specific parameter summaries for initial numbers, longline
selectivity, purse-seine selectivity, and recruitment deviations.

## Level 20: longline selectivity regularization scan

Level 19 showed that flexible longline age selectivity greatly improved the
longline age-composition fit while leaving initial numbers and recruitment
well behaved. Level 20 scans the longline selectivity prior width to find the
smallest amount of selectivity flexibility that preserves most of the Level 19
improvement.

Scanned values:

- sigma_ll_sel_dev = 0.50
- sigma_ll_sel_dev = 0.75
- sigma_ll_sel_dev = 1.00
- sigma_ll_sel_dev = 1.25

Decision rule: prefer the tightest prior that retains most of the objective and
residual improvement, rather than automatically choosing the loosest fit.

## Level 20 wide sigma scan

The first Level 20 scan improved monotonically through sigma_ll_sel_dev = 1.25,
so a wider scan was added for 1.25, 1.50, 1.75, and 2.00. This checks whether
there is an objective/residual elbow or whether the model continues to request
more longline selectivity flexibility.

## Level 21: age-based natural mortality diagnostic

Age-based natural mortality is now on the explicit model-development list. Level
21 starts from the Level 20 sigma_ll_sel_dev = 1.75 diagnostic compromise and
adds a restrained age-based M structure:

- ages 1–3 share a young-M multiplier
- ages 4–7 use the adult fixed-M anchor
- ages 8–10 share an old-M multiplier

Both M multipliers are penalized on the log scale with sigma = 0.35. This tests
whether age-based M improves the remaining fit without letting M become a fully
free age-specific selectivity surrogate.
