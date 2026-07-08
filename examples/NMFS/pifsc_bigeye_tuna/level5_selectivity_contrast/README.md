# Level 5: Strong Fleet Selectivity Contrast

Level 4 showed that distinct vulnerable-biomass index trajectories improved
optimization behavior but did not collapse the dominant weak direction. Wiggle
diagnostics showed that individual q and R0 parameters were highly influential
when moved alone, while fleet selectivity parameters could move with relatively
small objective cost.

This level tests whether the selectivity weakness is caused by insufficient
composition information.

## Model change

Compared with Level 4:

- keep the same model structure
- keep distinct vulnerable-biomass index trajectories
- make fleet age compositions strongly separated
- increase age-composition effective sample size in the objective

## Synthetic contrast

Longline:
- primarily older fish
- ages 6-10 dominate

Purse seine:
- primarily young fish
- ages 1-3 dominate

## Scientific hypothesis

If selectivity was weak because composition information was too diffuse or too
weakly weighted, then strong composition contrast should:

- increase sensitivity of selectivity wiggles
- reduce weak eigenvector loading on purse-seine a50
- improve optimization behavior
- clarify whether q/F/R0 confounding is separate from selectivity weakness

## Decision rule

If stronger composition information does not improve selectivity geometry, then
the weakness is structural rather than data-strength related.
