# Manuscript outline

## Working title

**Building a Model-Structure-Aware Statistical Inference Engine: Domain-Expert,
AI-Guided Development**

Alternative titles:

1. **From Statistical Model to Inference Engine: Domain-Expert, AI-Guided
   Development of Quadra**
2. **Quadra: A Model-Structure-Aware Engine for Exact Laplace Inference**
3. **Building Scientific Inference Software with a Domain Expert and an AI
   Coding Agent**
4. **Model Structure as a Computational Asset: Building an Exact Laplace
   Inference Engine with AI-Guided Development**

The recommended title keeps the technical contribution first and treats the
development process as a substantive secondary contribution. “Domain-expert”
is preferred over “domain-export.”

## Central thesis

A statistical inference engine can obtain substantial performance,
transparency, and reliability gains by treating model structure as an explicit
computational asset. Quadra demonstrates this through structure-aware automatic
differentiation, random-effect Hessian analysis, exact Laplace gradients, and
reusable derivative machinery. Its development also provides a concrete case
study of domain-expert, AI-guided scientific software engineering in which the
human defines scientific intent and acceptance criteria while an AI coding
agent accelerates implementation, diagnosis, testing, documentation, and
iteration.

## Intended contribution

The manuscript should make three contributions and keep them clearly
distinguished:

1. **Computational:** describe Quadra's model-structure-aware differentiation
   and exact Laplace inference architecture.
2. **Scientific-software:** show how invariants, numerical audits, benchmarks,
   and end-to-end examples turn a prototype into an inspectable inference
   workflow.
3. **Development-method:** document what AI guidance accelerated, what remained
   the responsibility of the domain expert, and which verification practices
   prevented plausible but incorrect changes from being accepted.

The paper should not claim that AI independently established statistical or
scientific validity.

## Target audience and possible venues

- Statistical computing and automatic-differentiation researchers.
- Fisheries stock-assessment and ecological-modeling practitioners.
- Scientific-software engineers.
- Researchers studying human–AI collaboration in technical development.

Potential venue classes include statistical software journals, computational
statistics journals, scientific-computing journals, and fisheries-methods
journals. The final balance among theory, software, fisheries application, and
AI-development process should be adjusted to the selected venue.

## Proposed manuscript structure

### 1. Abstract

- Problem: hierarchical scientific models require reliable and efficient
  inference over fixed and random effects.
- Gap: generic differentiation and default Laplace pathways may fail to exploit
  stable model structure or expose enough information for diagnosis.
- Approach: Quadra combines model-structure discovery, reusable derivative
  machinery, conditional random-effect optimization, and exact Laplace
  gradients.
- Demonstration: a spatial, multi-fleet, multi-season tuna assessment exercises
  the complete fit, simulation-estimation, sampling, projection, and reporting
  workflow.
- Development result: domain-expert direction paired with AI-assisted
  implementation shortened the diagnose–modify–validate loop, while numerical
  audits and scientific acceptance tests remained authoritative.
- Main results: report validated numerical agreement and measured changes in
  runtime, evaluation counts, and memory—not qualitative speed claims alone.

### 2. Introduction

#### 2.1 Motivation

- Random effects are central to contemporary ecological and fisheries models.
- Laplace approximation is attractive but its derivatives and inner mode solves
  can dominate runtime and complicate debugging.
- Model topology is often stable even when parameter values change; inference
  software should exploit that fact explicitly.

#### 2.2 The practical software problem

- A model may compile yet remain unusable because of excessive evaluations,
  opaque optimizer behavior, stale derivative state, or unstable nested solves.
- Low memory use is not necessarily evidence of efficiency; the relevant
  questions include redundant graph construction, derivative order, data
  locality, parallelism, and evaluation counts.

#### 2.3 Research questions

1. How can stable model structure be discovered and reused safely?
2. How can an exact Laplace gradient be evaluated without constructing
   unnecessary dense derivative objects?
3. Which tests distinguish a faster calculation from an incorrect one?
4. How can a domain expert and an AI coding agent collaborate on inference
   software without outsourcing scientific judgment?

#### 2.4 Contributions and paper organization

- State the three contributions above.
- Explicitly separate software capability, empirical evidence, and development
  process.

### 3. Statistical and computational background

#### 3.1 Hierarchical objective

- Define fixed effects $\theta$, random effects $u$, and joint negative log
  objective $J(\theta,u)$.
- Define the conditional mode $\hat u(\theta)$.

#### 3.2 Laplace objective

- Present

  $$
  L(\theta)=J(\theta,\hat u)+\frac{1}{2}\log|H_{uu}|-
  \frac{n_u}{2}\log(2\pi).
  $$

- Explain the roles of the inner Newton solve, Hessian factorization, and outer
  optimizer.

#### 3.3 Exact marginal gradient

- Separate the envelope contribution, implicit mode response, and derivative of
  the log determinant.
- Explain why mixed second and third derivatives appear.
- Contrast exact differentiation with finite-difference validation and with
  approximate alternatives.

#### 3.4 Model structure

- Define graph topology, Hessian sparsity pattern, random/fixed partition, and
  stable versus value-dependent structure.
- Explain when structure reuse is valid and when a tape or plan must be rebuilt.

### 4. Quadra architecture

#### 4.1 Design goals

- Model code templated on scalar type.
- Explicit parameter metadata and fixed/random partition.
- Inspectable inference stages rather than a monolithic optimizer call.
- Exactness first, with optimized pathways required to demonstrate parity.

#### 4.2 Automatic differentiation and reporting

- Scalar types and derivative orders.
- Model evaluation and named report values.
- Separation between model probability calculations and workflow/reporting
  products.

#### 4.3 Structure discovery and persistent derivative machinery

- Detection of the random-effect Hessian pattern.
- Persistent structure-aware tapes or graph plans.
- Reuse within a conditional mode solve.
- Fingerprinting, rebuild intervals, and safeguards against stale topology.

#### 4.4 Conditional random-effect optimization

- Newton iterations, factorization, damping, and backtracking.
- Full Newton trial steps with safeguards.
- Convergence and failure diagnostics.

#### 4.5 Exact Laplace-gradient contraction

- Let $T=D^3J$ and $H_{uu}^{-1}=LL'$.
- Present

  $$
  \mathrm{tr}(H_{uu}^{-1}\dot H_{uu}[d])
  =\sum_k T[d,l_k,l_k].
  $$

- Explain cubic polarization and the reduction from entry-wise dense
  construction to $2n_u+1$ third-directional evaluations per direction.
- Discuss dense and sparse regimes without implying that the current tuna
  example establishes performance for every sparsity pattern.

#### 4.6 Outer optimization

- Phased fitting and parameter masks.
- L-BFGS objective-only line searches.
- Exact-gradient evaluation at accepted iterates and accepted-objective reuse.
- Iteration and evaluation budgets.
- Human-readable stacked progress and convergence coloring as operational
  observability rather than an algorithmic contribution.

#### 4.7 Complete workflow

- Fit → simulation-estimation → report generation.
- Checkpoints, fingerprints, failure propagation, and output contracts.
- Optional posterior kernels and conditional random-effect reconstruction.

### 5. Domain-expert, AI-guided development method

#### 5.1 Roles and decision authority

- Domain expert: scientific requirements, model interpretation, acceptable
  approximations, validation thresholds, and final decisions.
- AI coding agent: repository inspection, hypothesis generation,
  implementation, test construction, documentation, and repetitive comparison.
- Version control and executable tests as the shared external memory.

#### 5.2 Iterative development loop

1. Observe a concrete failure or performance symptom.
2. Form a testable computational hypothesis.
3. Establish a numerical or behavioral baseline.
4. Make the smallest relevant implementation change.
5. Compare objectives, gradients, trajectories, and runtime.
6. Reject changes that are faster but numerically inconsistent.
7. Preserve accepted behavior in regression tests and documentation.

#### 5.3 Representative episodes

- Resolving compilation/API mismatch for optimizer iteration reporting.
- Diagnosing apparent infinite optimization through evaluation-level
  observability and budgets.
- Moving exact-gradient work out of rejected line-search trials.
- Increasing the inner Newton trial step while retaining safeguards.
- Replacing dense entry-wise third-order contraction with the exact factorized
  identity.
- Testing and rejecting a reverse-replay optimization when gradient parity did
  not hold.
- Extending the workflow from fitting to simulations, reports, and interactive
  scientific visualization.

#### 5.4 Why the process is scientifically defensible

- AI suggestions are hypotheses, not evidence.
- Independent finite differences, parity tests, simulation recovery, and
  end-to-end output checks provide evidence.
- Failed experiments should be reported where they reveal important boundary
  conditions or validation lessons.

#### 5.5 Reproducibility of AI involvement

- Record model/tool versions and repository commits.
- Summarize prompts or development episodes without treating conversational
  transcripts as sufficient provenance.
- Identify which code was AI-assisted and how it was reviewed.
- Discuss confidentiality, licensing, and authorship policy for the target
  venue.

### 6. Case study: spatial multi-fleet tuna assessment

#### 6.1 Why this model is a useful stress test

- Multiple fleets, regions, seasons, and ages.
- Seasonal movement and fleet-specific availability/selectivity.
- Annual recruitment random effects.
- Several observation likelihoods and management-derived quantities.

#### 6.2 Model specification

- Refer to the auditable model document for complete equations.
- Summarize population dynamics, observation processes, priors, and parameter
  transforms needed to understand the computational experiment.

#### 6.3 Inference workflow

- Phased fixed-effect optimization and conditional recruitment modes.
- Exact Laplace objective and gradient.
- Posterior sampling or exact-correction kernels where included in results.
- Simulation-estimation, reference points, and projections.

#### 6.4 Output and visualization

- Machine-readable output contract.
- Interactive spatial fishery pulse map: biomass, movement, fleet catch, and
  discards.
- Management-strategy race as a communication layer over deterministic
  projections.
- Make clear that schematic region placement is not geographic inference.

### 7. Verification and experimental design

#### 7.1 Correctness hierarchy

- Unit tests for derivative primitives and trace contractions.
- Exact-versus-reference parity tests.
- Finite-difference audits of marginal gradients.
- Optimizer trajectory and convergence checks.
- Model-level acceptance, reference-point, projection, and recovery tests.
- Full workflow and failure-propagation tests.

#### 7.2 Performance protocol

- Specify hardware, compiler, flags, thread counts, and software revision.
- Separate compilation time, wall time, CPU time, peak resident memory, and
  objective/gradient evaluation counts.
- Warm-up and repetition strategy.
- Use identical data, starting values, tolerances, and iteration caps.
- Report numerical agreement alongside every speed comparison.

#### 7.3 Comparators and ablations

- Non-reused versus reused structure/tape.
- Entry-wise versus factorized exact trace contraction.
- Gradient-based versus objective-only line-search trials.
- Conservative versus full Newton initial steps with identical safeguards.
- Serial versus parallel directional evaluation where appropriate.
- Default Laplace pathway or another established engine, if a fair and
  reproducible comparator can be configured.

#### 7.4 Metrics

- Objective and gradient agreement.
- Inner Newton iterations.
- Outer iterations and objective/gradient evaluations.
- Time per marginal evaluation and per accepted optimizer iteration.
- Peak memory.
- End-to-end workflow time.
- Simulation recovery and posterior diagnostic metrics where applicable.

### 8. Results

#### 8.1 Numerical validation

- Table of maximum absolute and relative errors for gradient and contraction
  tests.
- Optimizer endpoint and trajectory comparisons.

#### 8.2 Computational performance

- Baseline-to-final waterfall showing the incremental effect of each accepted
  optimization.
- Runtime, evaluation count, and peak-memory table.
- Scaling with random-effect dimension and Hessian structure, using synthetic
  benchmarks in addition to the tuna example.

#### 8.3 End-to-end tuna workflow

- Fit status and core scientific outputs.
- Simulation/recovery results.
- Reference points and projection products.
- Example frames or stills from the spatial animation.

#### 8.4 Development-process findings

- Time or iteration savings where records support them.
- Types of tasks well suited to AI assistance.
- Errors caught by tests, including the rejected reverse-replay attempt.
- Areas where domain expertise was indispensable.

### 9. Discussion

#### 9.1 Technical implications

- Stable model structure can be exploited without weakening exactness.
- Evaluation accounting and observability are first-class performance tools.
- Higher memory consumption is not itself a goal; computation should allocate
  memory only when caching or parallel work produces measured value.

#### 9.2 Implications for scientific software development

- AI can compress implementation and diagnostic cycles.
- Strong acceptance tests make rapid iteration safer and more informative.
- Human scientific ownership remains essential.

#### 9.3 Generalizability

- Other hierarchical models with repeated conditional mode solves.
- Limits imposed by dynamic control flow, changing sparsity, very large random
  fields, or unavailable third-order operations.

#### 9.4 Limitations

- The tuna example is synthetic and has two schematic regions.
- A single case study cannot establish broad performance superiority.
- Laplace approximation has statistical limitations independent of derivative
  correctness.
- AI-development observations are partly retrospective unless prospective
  process metrics are collected.
- Comparisons may be sensitive to compiler, hardware, and tuning choices.

#### 9.5 Future work

- Sparse inverse-factor contractions and larger spatial random fields.
- Broader benchmark suite and external engine comparisons.
- Automated structural-validity checks for longer tape reuse.
- Prospective measurement of human–AI development episodes.
- Additional real-data scientific-center examples.

### 10. Conclusion

- Reiterate that model structure is both a statistical concept and a
  computational resource.
- Summarize how exactness, performance, and inspectability were advanced
  together.
- Close with the development lesson: AI assistance is most useful when paired
  with domain ownership, explicit invariants, and executable evidence.

## Proposed figures

1. **Quadra inference architecture:** model evaluation → structure discovery →
   conditional mode solve → exact Laplace objective/gradient → optimizer or
   sampler.
2. **Nested optimization sequence:** accepted outer iterates, objective-only
   line-search trials, and inner random-effect Newton solves.
3. **Factorized trace contraction:** comparison of dense entry-wise and
   inverse-factor directional calculations.
4. **Development loop:** domain requirement → AI-assisted implementation →
   numerical audit → accept/reject → regression test.
5. **Performance waterfall:** cumulative runtime improvement from each accepted
   optimization.
6. **Scaling results:** runtime and memory versus random-effect dimension and
   Hessian structure.
7. **Tuna model schematic:** ages, seasons, regions, fleet removals, movement,
   and recruitment.
8. **Spatial fishery pulse stills:** selected historical frames and management
   scenarios.

## Proposed tables

1. Quadra capabilities and associated validation evidence.
2. Benchmark environment and reproducibility settings.
3. Exact-gradient and third-order contraction agreement.
4. Runtime, evaluation count, and peak-memory ablations.
5. Tuna assessment configuration and dimensions.
6. Scientific workflow outputs and acceptance status.
7. Human and AI responsibilities during representative development episodes.
8. Failed or rejected optimizations and the tests that rejected them.

## Supplementary material

- Complete tuna model equations and configuration.
- Derivation of the exact factorized trace contraction.
- Benchmark source code and raw timing data.
- Gradient-audit and parity-test outputs.
- Simulation-recovery protocol.
- Full output schema.
- AI-assistance disclosure and reproducibility statement.

## Evidence still needed before drafting results

- A clean, versioned benchmark matrix run on identified hardware.
- Peak resident-memory measurements, not observations from an activity monitor.
- Repeated timings with variability estimates.
- Scaling cases beyond the 12 recruitment random effects in the tuna example.
- A fair comparator configuration if comparative claims are retained.
- Archived numerical parity and gradient-audit outputs for the final revision.
- A concise chronology of AI-guided development episodes and rejected changes.
- Target-journal authorship and generative-AI disclosure requirements.

## Recommended drafting order

1. Statistical background and Quadra architecture.
2. Verification and benchmark protocol.
3. Tuna case study.
4. Results after the benchmark matrix is frozen.
5. AI-guided development method and discussion.
6. Introduction, abstract, and conclusion last.
