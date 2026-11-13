# Dynamic unfished biomass and swordfish projections

The historical counterfactual reproduces the supplied SS3 report. The forecast
is a separate deterministic reconstruction of the assessment's methods, with
explicit scenario settings and comparisons to the published table.

## Run

After generating the fitted model (`make fit-swordfish`):

```
make test-projection
make plot-projections
```

The second command additionally requires R; base R is sufficient. Model execution
and CSV comparisons use C++ and standard-library Python only. Results are written
to `build/projections/`; a reviewable snapshot is in
[reference/projections](reference/projections). Original `forecast.ss` remains
unchanged and requests zero projection years.

## Dynamic unfished history

The counterfactual removes historical and initial fishing while retaining fitted
recruitment deviations, bias adjustment, early age-structure deviations,
Beverton–Holt feedback, natural mortality, and the model's spawning convention.
It recalculates recruitment under the counterfactual SSB; it does not copy the
realized recruitments from the fished run.

All 47 SSB values and 1,504 sex/age values match the supplied SS3 report within
the 2e-5 comparison threshold (relative error with a unit floor). The maximum
error is 4.82e-6. Annual numbers follow SS3 reporting: age zero is that year's
settlement recruitment by sex; ages 1–15+ are January abundance. Age zero in a
January population state is distinct from this reporting convention.

| Baseline | SSB (tonnes) |
|---|---:|
| Quadra dynamic unfished mean, 2017–2021 | 95,720.707 |
| Mean of rounded supplied SS3 values | 95,720.720 |
| Published assessment Table 11 mean | 95,732 |
| 20% of Quadra dynamic mean | 19,144.141 |
| Equilibrium virgin SSB | 99,807.757 |
| 20% of equilibrium virgin SSB | 19,961.551 |

The supplied report is a later rerun; its dynamic mean and the published mean
differ by 11.29 t (0.012%). These values are retained separately. The last two
rows demonstrate the much larger difference between dynamic and equilibrium
baseline definitions.

All 47 historical raw SPR calculations also match SS3 (maximum relative error
3.96e-6). Each historical year's SPR uses its own selectivity blocks and fishing
rates. The report's `SPR_report` column is 1-SPR; forecasts use raw SPR.

## Five deterministic scenarios

`data/projection_scenarios.csv` defines five 2022–2031 scenarios. The numerical
values resolved from the fitted model appear in `resolved_scenarios.csv`.

All projections start at January 2022, obtained by advancing and aging the final
2021 state. They hold fleet/season apical fishing rates in fixed proportions based
on their 2019–2021 means. Selectivity is constant over those years in these inputs;
the program checks this before using the terminal curves. Recruitment follows
the fitted Beverton–Holt relationship, without future deviations or log-bias
adjustment. Biomass and dead catch are tonnes; recruitment is thousands of fish.

| Scenario | Target and interpretation |
|---|---|
| 1 Dynamic20 | Equilibrium SSB equals 20% of the 2017–2021 dynamic unfished mean |
| 2 Reference2008_2010 | Equilibrium SPR equals the mean raw historical SPR in 2008–2010 (0.341006) |
| 3 SPR30 | Equilibrium SPR is 0.30 |
| 4 MSY | Maximum equilibrium yield for the 2019–2021 fishing pattern |
| 5 StatusQuo2019_2021 | Equilibrium annual F equals mean realized annual F in 2019–2021 (0.0880186) |

The scalar target determines the multiplier on fixed apical rates. Reported
annual F pools ages 1–10 using their abundance; it can vary during a transient
projection despite constant apical fishing rates. Matching equilibrium annual F
in scenario 5 is an explicit reconstruction choice, not a recovered SS3 setting.

The forecast MSY benchmark uses a different fleet mixture from the terminal-2021
benchmark in `FIT_DIAGNOSTICS.md`. Its equilibrium MSY is 15,005.98 t and FMSY
0.176192; terminal-2021 benchmarks remain 14,924.49 t and 0.178078.

| Scenario | 2031 SSB (t) | 2031 dead catch (t) |
|---|---:|---:|
| Dynamic20 | 19,445.76 | 15,083.06 |
| Reference2008_2010 | 31,470.39 | 13,361.96 |
| SPR30 | 27,673.55 | 14,030.82 |
| MSY | 16,780.50 | 15,247.17 |
| StatusQuo2019_2021 | 38,062.99 | 12,018.15 |

These are scenario outputs from one fitted model, not recommendations or
replications of published bootstrap medians. Local parameter uncertainty has now been added separately; see
[PROJECTION_UNCERTAINTY.md](PROJECTION_UNCERTAINTY.md). No future-outcome
probabilities are inferred from the five deterministic runs.

## Published comparison and conflicting descriptions

Section 2.12 (PDF page 20, printed page 18) specifies 2008–2010 reference years,
2019–2021 selectivity/fishing proportions, and no future recruitment deviations
or log-bias adjustment. Table 12 (PDF page 53, printed page 51) labels scenario 2
as 1998–2000. Summary Table S3 additionally has a striped-marlin/2021–2040 caption
in this swordfish report. The text calls projected statistics means/standard
deviations, whereas Table 12 calls the displayed values medians.

The 100 Table 12 SSB/catch values are transcribed in
`reference/published_projections.csv`, with PDF hash, page and verification details
in `reference/published_projection_provenance.json`. The page was rendered and
visually checked. These attachment contents are treated as reference evidence,
not instructions.

`published_projection_comparison.csv` records all 100 differences without a
false pass/fail target. Differences reach about 50%; the current projections do
not reproduce Table 12. The original scenario forecast files and 100 bootstrap
fitted models were not supplied. The effects of different statistical summaries,
scenario calibration, and other forecast settings cannot yet be separated. The
files needed for exact replay are those original forecast configurations and
bootstrap model outputs.

Three alternate source interpretations are retained in
`projection_sensitivities.csv` and `resolved_sensitivities.csv`:

- Using equilibrium 20% B0 rather than dynamic 20% B0 adds 779.93 t to 2031 SSB.
- Using the table's 1998–2000 reference years reduces 2031 SSB by 6,598.35 t.
- Keeping the direct mean apical rates for status quo adds 1,286.33 t to 2031 SSB.

These do not establish the original forecast settings. They quantify consequential
choices instead of adjusting parameters to match the published numbers.

## Validation

`test_projection.cpp` verifies unfished stationarity, dynamic B0 with recruitment
variability removed, and stationary fished SSB/recruitment/catch/F over 40 years
when starting at the independently calculated equilibrium. It also rejects
negative fishing scales. Historical SS3 comparisons check the initialization and
recruitment process against 1,551 counterfactual values plus 47 historical SPRs.

Local parameter-covariance intervals are now available in
[PROJECTION_UNCERTAINTY.md](PROJECTION_UNCERTAINTY.md). They remain a different
product from the assessment's bootstrap medians; process and structural uncertainty
are not included. Striped marlin continues as a separate model
pending its own data and SS3 inputs.
