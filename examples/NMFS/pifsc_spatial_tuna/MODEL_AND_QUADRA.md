# The complete Quadra spatial tuna assessment

## Purpose and scientific scope

This document is the specification of the model implemented in this repository.
It is deliberately more detailed than a normal README. A reader should be able
to reconstruct the state transitions, joint probability model, optimization,
Laplace approximation, posterior simulation, management calculations, and file
outputs without guessing what “the tuna model” or “Quadra” means.

The program is a seasonal, age-structured, multi-region, multi-fleet statistical
catch-at-age assessment. It estimates population scale, stock–recruit
productivity, fishery catchability, survey/index catchability, selectivity,
retention, observation error, movement, optional availability multipliers, and
annual recruitment deviations. Recruitment deviations are Gaussian random
effects. Quadra supplies automatic differentiation, optimization, and the
Laplace marginalization machinery. The assessment executable supplies the
biology, observation models, workflow, exact posterior target, samplers,
diagnostics, reference points, and projections.

This is currently a worked end-to-end assessment and testing platform, not a
completed operational assessment of a named tuna stock. The executable creates
a 12-year, 6-age example data object in
`src/advanced_tuna_spatial_assessment_example.cpp`. Configuration controls the
workflow and model switches; it does not yet provide a general external data
schema. The `data/opal_raw` audit path is separate from the active example-data
construction. Results must not be represented as management advice for a real
stock until real data ingestion, units, fleet definitions, and scientific
review are completed.

## Code map: where each idea lives

| Concern | Authoritative implementation |
|---|---|
| Data dimensions, storage order, validation | `include/tuna/tuna_spatial_data.hpp` |
| Parameters, transforms, dynamics, likelihoods, reports | `include/tuna/tuna_spatial_assessment_model.hpp` |
| Phased fitting, diagnostics, simulation/recovery | `include/tuna/tuna_assessment_acceptance.hpp` |
| Equilibrium reference points and projections | `include/tuna/tuna_reference_points.hpp` |
| Runtime configuration and end-to-end orchestration | `src/advanced_tuna_spatial_assessment_example.cpp` |
| Native nonlinear transport map | `include/tuna/dependency_free_transport_flow.hpp` |
| Native flow training and QFLOW serialization | `src/train_dependency_free_transport_flow.cpp` |
| Build/run/test entry points | `Makefile` |
| Default runtime choices | `config/tuna_assessment.conf` |

The model objective intentionally starts as an assessment narrative:
`evaluate_impl` selects the observation processes for the current fitting
phase, validates the objective inputs, and calls `evaluate_assessment_story`.
Within that story the code reads parameters, transforms them, adds priors and
smoothing penalties, initializes the population, advances years and seasons,
predicts observations, evaluates likelihoods, recruits and ages fish, and
reports derived quantities. Small member functions name the biological pieces:
`bounded_age50`, `logistic_age`, `normalized_logistic_age`,
`beverton_holt_recruitment`, `initialize_state`, and `spawning_biomass`.

## Indices, arrays, and timing

Use $y=0,\ldots,Y-1$ for years, $s=0,\ldots,S-1$ for seasons,
$r,j=0,\ldots,R-1$ for regions, $a=0,\ldots,A-1$ for age array positions,
and $f=0,\ldots,F-1$ for fleets. Displayed biological age is $a+1$. The last
age is a plus group.

The population state $N_{y,s,r,a}$ is numbers at age immediately before the
fishing and natural mortality calculation in season $s$. Within a season the
order is:

1. calculate fleet fishing rates from the pre-season state;
2. calculate capture and survival during the season;
3. evaluate observations predicted for that fleet-season-region;
4. move survivors among regions.

After the final season, recruitment for the next year is calculated, survivors
age, and the oldest fish accumulate in the plus group. The code records annual
spawning biomass before the first seasonal mortality calculation. The spawning
biomass function nevertheless applies survival to `spawning_fraction`; this is
an explicit within-year spawning-timing adjustment.

Flattened input orders are part of the interface:

| Quantity | Storage order |
|---|---|
| weight at age | `[year, age]` |
| movement | `[season, from_region, to_region]` |
| effort and scalar observations | `[fleet, year, season, region]` |
| catch numbers | `[fleet, year, season, region, age]` |
| fixed availability | `[fleet, season, region, age]` |

`TunaSpatialAssessmentData::validate()` checks dimensions, finite/range
constraints, recruitment proportions and movement rows. Missing scalar
observations are represented by non-positive values and omitted from their
likelihoods. A composition is omitted when its total observed count is zero.

## Inputs and fixed biological quantities

The model treats these as data, not estimated parameters:

- annual natural mortality $M_a$;
- maturity $m_a$;
- year-specific weight $w_{y,a}$;
- regional recruitment shares $\rho_r$, which sum to one;
- effort $E_{f,y,s,r}$;
- optional availability surface $A_{f,s,r,a}$, otherwise one;
- the spawning fraction $\tau$, constrained to $0<\tau<1$;
- observed indices, retained biomass, discarded biomass, catch-at-age counts,
  and optional total catch used for catch conditioning.

The supplied movement matrix initializes movement logits. It is not a fixed
movement matrix once movement parameters are estimated.

## Estimated parameters and transformations

All optimizer and sampler coordinates are unconstrained. The biological values
are:

$$
R_0=\exp(\ell_{R_0}),\qquad
h=0.2+0.8\mathrm{logit}^{-1}(\eta_h),\qquad
\sigma_R=\exp(\ell_{\sigma_R}).
$$

Thus steepness lies strictly between 0.2 and 1. For every fleet:

$$
q_f=e^{\ell_{q_f}},\quad q^I_f=e^{\ell_{q^I_f}},\quad
\sigma^I_f=e^{\ell_{\sigma^I_f}},\quad
\theta_f=e^{\ell_{\theta_f}},\quad
\sigma^{ret}_f=e^{\ell_{\sigma^{ret}_f}},\quad
\sigma^{dis}_f=e^{\ell_{\sigma^{dis}_f}}.
$$

The raw age-50 transform is

$$
a_{50}(x)=1+(A-1)\mathrm{logit}^{-1}(x),
$$

so selectivity and retention age-50 values stay between ages 1 and $A$.
Slopes are exponentiated. Fleet selectivity is an ascending logistic normalized
to equal one at the oldest modeled age:

$$
s_{f,a}=\frac{[1+e^{-k^s_f((a+1)-a^s_{50,f})}]^{-1}}
{[1+e^{-k^s_f(A-a^s_{50,f})}]^{-1}+10^{-12}}.
$$

Retention is an unnormalized ascending logistic,

$$
p^{ret}_{f,a}=[1+e^{-k^{ret}_f((a+1)-a^{ret}_{50,f})}]^{-1}.
$$

If availability scales are estimated, their exponentiated value is either one
per fleet or one per fleet-season-region. Effective availability is

$$
V_{f,s,r,a}=A_{f,s,r,a}\exp(\ell^V_{f,s,r}).
$$

For each season and origin, $R-1$ movement logits are free and the last
destination has reference logit zero. Softmax gives a row-stochastic matrix:

$$
P_{s,r,j}=\frac{e^{\gamma_{s,r,j}}}
{1+\sum_{k=0}^{R-2}e^{\gamma_{s,r,k}}},\quad j<R-1,
\qquad
P_{s,r,R-1}=\frac1{1+\sum_{k=0}^{R-2}e^{\gamma_{s,r,k}}}.
$$

Movement may be shared across seasons. Annual recruitment deviations
$\epsilon_y$ are the model's random effects; their indices are exposed to
Quadra separately from fixed effects.

## Initial population and unfished spawning biomass

At the start of year 1, recruitment in region $r$ is bias-corrected lognormal:

$$
N_{0,0,r,0}=R_0\rho_r\exp(\epsilon_0-\sigma_R^2/2).
$$

Older ages are generated by applying annual natural survival successively.
The plus group divides incoming survivors by $1-e^{-M_{A-1}}$. This is an
equilibrium-style initialization under natural mortality, with the first
recruitment deviation applied to the entire initialized cohort scale through
the age recursion.

For any state $N$, year-$y$ spawning biomass is

$$
SSB_y=\sum_{r,a}N_{r,a}w_{y,a}m_a e^{-M_a\tau}.
$$

`ssb0` is this quantity from the initialized state using year-1 weights. It is
the $SSB_0$ appearing in the stock–recruit curve; it is therefore a derived
quantity conditional on the initialized state, not an independently estimated
parameter.

## Seasonal fishing, catch, survival, and movement

In the effort-driven model, instantaneous fleet fishing mortality is

$$
F_{f,y,s,r,a}=q_f E_{f,y,s,r}s_{f,a}V_{f,s,r,a},
\qquad Z_{y,s,r,a}=M_a/S+\sum_fF_{f,y,s,r,a}.
$$

The Baranov capture equation is

$$
C_{f,y,s,r,a}=N_{y,s,r,a}F_{f,y,s,r,a}
\frac{1-e^{-Z_{y,s,r,a}}}{Z_{y,s,r,a}},
$$

and survivors are $N^S_{r,a}=N_{r,a}e^{-Z_{r,a}}$. Predicted retained and
discarded biomass are

$$
B^{ret}_{f,r}=\sum_a C_{f,r,a}p^{ret}_{f,a}w_{y,a},\qquad
B^{dis}_{f,r}=\sum_a C_{f,r,a}(1-p^{ret}_{f,a})w_{y,a}.
$$

The index's vulnerable biomass is

$$
B^V_{f,r}=\sum_aN_{r,a}w_{y,a}s_{f,a}.
$$

Notice that effective availability enters fishing mortality but does not enter
this vulnerable-biomass sum. Index availability is represented through the
fleet index catchability, not the fishing availability surface.

With catch conditioning, a fleet-season-region removal scale is observed total
catch divided by vulnerable numbers or biomass (according to fleet catch
units). Capture becomes $C=N F$, without the Baranov exploitation multiplier.
Survival uses a smooth positive approximation to $1-\sum_fF_f$, multiplied by
natural survival. This makes the conditioned rate a removal fraction rather
than instantaneous fishing mortality. The distinction is retained in
reference-point calculations.

Movement occurs after survival:

$$
N'_{j,a}=\sum_rN^S_{r,a}P_{s,r,j}.
$$

Consequently movement conserves post-mortality numbers, subject to floating
point error, because each origin row of $P$ sums to one.

## Recruitment and annual aging

After all seasons, Beverton–Holt expected recruitment is

$$
\bar R(SSB)=\frac{4hR_0SSB}
{SSB_0(1-h)+SSB(5h-1)+10^{-12}}.
$$

For the next year:

$$
N_{y+1,0,r}=\rho_r\bar R(SSB_y)
\exp(\epsilon_{y+1}-\sigma_R^2/2).
$$

All post-season survivors age one class. Fish already in the final class and
fish aging into it are added together. Recruitment deviations have density
$\epsilon_y\sim N(0,\sigma_R)$. The subtraction of $\sigma_R^2/2$ makes the
multiplicative recruitment deviation mean one on the arithmetic scale.

## Observation models

Positive index observations follow a lognormal model expressed as a normal
density on log observations:

$$
\log I_{f,y,s,r}\sim N(\log(q^I_fB^V_{f,y,s,r}+10^{-12}),\sigma^I_f).
$$

Retained and discarded biomass use analogous log-scale normals centered on
their predictions, with fleet-specific standard deviations. Because the
normal density is evaluated on the logged datum without an explicit Jacobian,
the objective is the likelihood for log-transformed observations; for fixed
data the omitted Jacobian is constant in parameters.

Catch-at-age proportions are predicted by normalizing capture numbers within a
fleet-year-season-region. Counts use a Dirichlet-multinomial distribution with
fleet concentration parameter $\theta_f$. The entire contribution is
multiplied by `composition_likelihood_weight`. Zero-total rows contribute
nothing.

Fitting phases change which observation components are active:

| Phase | Index | Composition | Retained | Discard |
|---|---:|---:|---:|---:|
| InitializeRecruitment | configured | off | off | off |
| InitializeCatchability | configured | off | configured | configured |
| InitializeMovement | configured | configured | configured | configured |
| Full | configured | configured | configured | configured |

Parameter locking during phases is implemented in the fitting layer, separately
from this likelihood switch table.

## Priors, regularization, and complete joint objective

When priors are enabled, independent normals are placed on log catchabilities,
raw age-50 parameters, log slopes, log observation standard deviations,
availability logs, and movement logits. Their means and standard deviations are
the `TunaAssessmentControls` values. In particular the retention log-slope mean
is -0.2 and observation log-standard-deviation mean is -1. The composition
concentration, $R_0$, steepness, and recruitment standard deviation do not
receive corresponding explicit fixed-effect priors in the current objective.

Quadratic first-difference penalties smooth availability across adjacent
seasons and movement logits across adjacent parameter seasons:

$$
Q_V=\frac{w_V}{2}\sum(\ell^V_{s}-\ell^V_{s-1})^2,\qquad
Q_P=\frac{w_P}{2}\sum(\gamma_s-\gamma_{s-1})^2.
$$

The joint negative log objective is exactly

$$
J(\theta,\epsilon)=NLL_{prior}+Q_V+Q_P+NLL_I+NLL_{comp}
+NLL_{ret}+NLL_{dis}+NLL_R.
$$

The reported `nll_decomposition_total` is the sum above and is intended to
match the returned objective. Each component is also reported separately.

## What Quadra does

Quadra is the C++ differentiation and inference engine used by this assessment.
The tuna model is templated on scalar type, so the same objective can run with
ordinary doubles or Quadra automatic-differentiation scalars. That supplies
gradients, Hessian information, and the higher directional derivatives needed
by the exact derivative of the Laplace objective.

Let $\theta$ denote fixed effects and $u=\epsilon$ recruitment random effects.
For each fixed-effect vector Quadra solves

$$
\hat u(\theta)=\arg\min_u J(\theta,u)
$$

with Newton iterations using the random-effect Hessian
$H_{uu}=\partial^2J/\partial u\partial u'$. The Laplace negative log marginal
objective is, up to the configured normalizing constant,

$$
L(\theta)=J(\theta,\hat u)+\tfrac12\log|H_{uu}(\theta,\hat u)|
-\tfrac{n_u}{2}\log(2\pi).
$$

The structure-aware machinery is active: a persistent random-Hessian tape
captures the objective's derivative topology and is reused across inner Newton
iterations. The default rebuild interval is one marginal evaluation because
longer cross-evaluation reuse exposed stale frozen-topology reads in this model.
This still retains the important reuse inside each mode solve. Mixed derivatives
and third directional derivatives provide the exact Laplace gradient when a
gradient-based consumer requests it. Density-only samplers omit that work.

Fixed effects sampled in the posterior exclude the configured locked parameters
(`logit_steepness` and one anchor fleet's index catchability by default). The
partition, parameter names/order, data, controls, fitted point, and geometry are
fingerprinted to prevent stale checkpoint or transport reuse.

## Optimization and phased fitting

The workflow can fit sequential recruitment, catchability, movement, and full
phases. Each phase exposes an intended parameter subset and carries its result
forward. L-BFGS settings, iteration limits, starts, seeds, and whether to run
directly in the full phase are runtime configuration. Multistart fitting is an
explicit comprehensive-analysis option. Checkpoints record values and metadata;
loading is rejected when fingerprints or parameter identity disagree.

The initial point and fitted point are evaluated through the same model. Output
contains objective/gradient information, likelihood decomposition, parameter
estimates, biomass trajectory, movement, and phase diagnostics. Acceptance
checks are scientific workflow gates, not additions to the probability model.

## Simulation and recovery testing

Two related concepts must not be confused:

1. `data.model_consistent=true` constructs the worked example using predictions
   intended to be internally consistent with the model.
2. The simulation-estimation loop takes a complete truth data object, perturbs
   its already-populated observations, refits, and compares estimated with truth
   terminal depletion.

Indices and biomass are multiplied by independent mean-corrected lognormal
noise. Catch-at-age cells are independently Poisson sampled with their existing
cell count as the mean (bounded below by one for the Poisson mean). The current
`simulate_observed_data` helper therefore does not draw a new latent population
trajectory or multinomial/Dirichlet-multinomial composition conditional on a
fixed row total. It is an observation-perturbation recovery experiment. Seeds
are deterministic. Summaries include mean, median, 10th/90th percentile
depletion bias, finite replicate count, and the rate estimated below the low
depletion threshold.

This distinction is important: successful recovery tests demonstrate behavior
under this particular perturbation experiment, not universal identifiability or
frequentist coverage.

## Posterior target and samplers

The retained fixed-effect target is the Laplace-marginal posterior proportional
to $\exp[-L(\theta)]$, with locked coordinates omitted. Available kernels are:

- AD-NUTS, which requests exact marginal gradients and adapts during warmup;
- pCN in frozen marginal-Hessian-whitened coordinates;
- a frozen full-covariance symmetric random walk;
- experimental Gaussian-mixture, KDE, and polynomial-flow independence
  proposals with exact Metropolis correction;
- native QFLOW RealNVP independence Metropolis;
- native QFLOW mixture i-SIR, the preferred nonlinear-proposal sampler after
  its flow ensemble has been trained and validated.

The optional `transport_isir` transition includes the current state among $K$
candidates, draws $K-1$ candidates from an equal-weight mixture of frozen
flows, computes exact weights

$$
w_k\propto \frac{\exp[-L(\theta_k)]}{q_{mix}(\theta_k)},
$$

and resamples one candidate proportional to those weights. Including the
current state makes this an exact invariant kernel for the same target; the
flows improve proposals but never replace the assessment likelihood.

QFLOW is a dependency-free C++/Eigen RealNVP format. Training uses masked affine
coupling layers, whole-chain holdout validation, deduplication, early stopping,
gradient clipping, AdamW, and deterministic seeds. Manifests bind artifacts to
parameter names/order, assessment and geometry fingerprints, source hashes,
architecture, validation likelihood, and round-trip/parity checks. A mismatch
fails closed.

Sampler diagnostics include chain movement, acceptance/kernel diagnostics,
R-hat, bulk ESS, tail ESS, target-evaluation cost, and profiling of mode solving
versus exact-gradient work. `sampler_summary.csv` is authoritative;
`nuts_summary.csv` is only a compatibility copy.

## Reconstructing recruitment effects

Laplace sampling produces fixed-effect draws, not joint draws. For every
retained fixed-effect draw, the workflow re-solves the conditional recruitment
mode and draws

$$
u\mid\theta,y\ \approx\ N(\hat u(\theta),H_{uu}^{-1}).
$$

The mode and Gaussian draw are both saved. Failures are counted. Reference
points and projections use a reproducibly thinned subset of these reconstructed
joint draws; they do not silently hold recruitment deviations at the fitted
mode.

## Reference points and projections

A fishing pattern is built from fitted recent fishing. In the effort model it
is instantaneous $F$-at-age; in catch-conditioned mode it is explicitly
tagged as a removal-fraction pattern. For each fishing multiplier, the reference
engine repeats seasonal fishing, survival, movement, Beverton–Holt recruitment,
and aging until the maximum relative state change is below $10^{-10}$, after
at least 20 and at most 1000 iterations. Terminal-year weights are used.

The unfished equilibrium gives $B_0$. A fishing-multiplier grid is searched
for maximum equilibrium total yield (retained plus discard), producing MSY,
$B_{MSY}$, and the $F_{MSY}$ multiplier. Boundary solutions are flagged.
Reported status includes terminal $B/B_{MSY}$ and status-quo
$F/F_{MSY}$. Projection scenarios apply specified fishing multipliers and
save annual spawning biomass, depletion, retained yield, discard yield, and
total yield for each management draw.

These are deterministic equilibrium and projection calculations conditional
on a parameter/random-effect draw. They do not include implementation error,
future environmental covariates, or a separate management procedure unless a
scenario explicitly supplies such behavior.

## Configuration and reproducibility

Runtime configuration is dependency-free `key = value` text. Blank lines,
comments, and section labels are accepted. Precedence is

```text
compiled defaults < configuration file < QUADRA_TUNA_* environment variables
```

The resolved configuration is saved with results. Important sections are
`run`, `output`, `data`, `model`, `fit`, `sampling`, and `simulation`. The
default file is `config/tuna_assessment.conf`; another file can be passed with
`--config` or `make ... CONFIG=...`. Configuration makes workflow choices
runtime choices, but, as noted above, the complete biological data object is
still constructed in C++.

Primary commands are:

```bash
make fit-advanced-tuna
make sample-advanced-tuna
make run-comprehensive-analysis
make test-fast
make test
make train-transport-flow
make check-transport-gradient
make check-transport-reproducibility
make check-posterior-assessment
```

The fit and sampling commands generate the report after the model run. The full
simulation recovery test is cached after success; `force-simulation-recovery`
forces a rerun.

## Output contract

Machine-readable files belong under `build/assessment_outputs/data`; rendered
Markdown and figures belong under `build/assessment_outputs/report`. The exact
set depends on whether fitting, sampling, sensitivity, retrospective, or
simulation work was requested. Core products include:

| Product | Meaning |
|---|---|
| `effective_configuration.csv` | resolved runtime settings |
| fit/checkpoint files | fitted coordinates and identity/fingerprint metadata |
| `likelihood_decomposition.csv` | objective components defined above |
| biomass/movement/parameter CSVs | fitted derived state and estimates |
| `posterior_draws.csv` | retained Laplace-marginal fixed-effect draws |
| `sampler_summary.csv` | convergence and efficiency by parameter |
| `sampler_identity.csv` | kernel identity and exact evaluation accounting |
| `marginal_sampler_profile.csv` | mode, gradient, and overhead timing |
| `marginal_whitening_cache.csv` | validated frozen proposal geometry |
| `posterior_random_effect_draws.csv` | conditional modes and reconstructed effects |
| `posterior_reconstruction_summary.csv` | reconstruction success/failure counts |
| `posterior_reference_points.csv` | draw-wise equilibrium management quantities |
| `posterior_projection_draws.csv` | scenario trajectories by draw |

`make check-posterior-assessment` requires convergence, finite and complete
latent reconstruction, valid reference points, four complete biologically
admissible projection scenarios, and the CSV/report directory separation.

## Diagnostics and interpretation

Reported `r0`, steepness, recruitment SD, SSB, depletion, movement, terminal
numbers, regional biomass, catchabilities, and likelihood components are
derived from the exact objective evaluation. Observation predictions can be
enabled for row-level residual auditing. Identifiability deserves particular
attention where catchability and availability multiply, and where survey
catchability and population scale trade off; anchoring and priors are part of
the current strategy, not proof that all contrasts are data-identified.

R-hat and ESS assess sampled marginal chains. They do not validate biological
assumptions, the Laplace approximation, the observation model, or real-data
quality. Simulation recovery is complementary but limited as described above.
The likelihood decomposition, phase behavior, retrospective patterns,
sensitivity runs, posterior predictive behavior, and raw-data provenance all
need review before scientific use.

## Known assumptions and deliberate limitations

- Natural mortality, maturity, weights, recruitment allocation, and the fixed
  availability surface have no estimated uncertainty.
- Selectivity and retention are ascending logistic curves; domed shapes and
  time variation are absent.
- Movement is age invariant and either seasonal or shared across seasons.
- Recruitment deviations are independent with a common SD; there is no AR(1)
  process or environmental covariate.
- Sexes are not modeled separately.
- Discard mortality is not separated from capture: captured fish contribute to
  fishing removal before being classified retained/discarded.
- Scalar observation errors are independent log-scale normals.
- Composition overdispersion is fleet-specific but otherwise shared across
  years, seasons, and regions.
- The assessment $SSB_0$ is derived from its initialized state; reference-point
  code independently iterates equilibrium dynamics.
- Laplace integration and conditional Gaussian reconstruction are local
  approximations in recruitment-effect space.
- The worked example is not yet a general data-driven stock-assessment package.

## Symbol glossary

| Symbol | Definition |
|---|---|
| $N$ | numbers at age by region |
| $R_0$ | unfished recruitment scale |
| $h$ | Beverton–Holt steepness |
| $SSB,SSB_0$ | spawning biomass and initialized unfished spawning biomass |
| $M_a$ | annual natural mortality at age |
| $F_f,Z$ | fleet fishing rate and total seasonal mortality |
| $q_f,q^I_f$ | fishing and index catchability |
| $s_{f,a}$ | normalized selectivity |
| $p^{ret}_{f,a}$ | retention probability |
| $V_{f,s,r,a}$ | effective fishing availability |
| $P_{s,r,j}$ | movement probability from $r$ to $j$ |
| $\epsilon_y,\sigma_R$ | recruitment deviation and its SD |
| $\theta$ | fixed effects collectively; also contextually composition concentration when subscripted by fleet |
| $u$ | random effects collectively |
| $J$ | joint negative log objective |
| $L$ | Laplace-marginal negative log objective |
| $H_{uu}$ | random-effect Hessian/conditional precision |

## End-to-end assessment story in one paragraph

The executable resolves and records configuration, constructs and validates the
data, defines unconstrained fixed and random parameters, fits them in controlled
phases, and evaluates a population initialized near unfished equilibrium. In
each year it computes spawning biomass; in each season it exposes fish to
fleet-, age-, region-, and availability-specific fishing plus natural mortality,
predicts catch composition, index, retained biomass, and discard biomass,
moves survivors, then after the last season generates Beverton–Holt recruitment
with a bias-corrected random deviation and ages survivors into a plus group.
Quadra differentiates the resulting joint negative log likelihood, conditionally
optimizes annual recruitment effects, and supplies a Laplace-marginal target for
optimization and exact-correction posterior kernels. The workflow diagnoses the
fit and chains, reconstructs conditional recruitment effects for every retained
draw, computes equilibrium reference points and scenario projections, validates
the products, writes CSVs separately from the rendered report, and records the
configuration and fingerprints needed to determine exactly what was run.
