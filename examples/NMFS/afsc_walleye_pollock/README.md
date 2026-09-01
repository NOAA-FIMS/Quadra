# Synthetic AFSC Walleye-Pollock-Style Example

This is a synthetic, public-data-safe Quadra example inspired by Alaska
walleye pollock assessment structure. It is not an official assessment.

## Model scope

- catch, index, and age-composition observations
- 5 fixed effects
- recruitment deviations as random effects
- Laplace fit through Quadra

The default runner currently executes the fast Level-1 configuration, which
disables random recruitment effects. Use the scaling diagnostic to exercise
configurations with up to 20 random effects.

## Build and run

From anywhere in this repository checkout:

```bash
./examples/NMFS/afsc_walleye_pollock/run_walleye_pollock_example.sh
```

The runner requires Bash and a C++17 compiler. Set `CXX` to override the
compiler. It writes the executable to `build/examples/afsc_walleye_pollock`
and CSV diagnostics to `examples/NMFS/afsc_walleye_pollock/outputs/`.

For the longer random-effect scaling sweep:

```bash
./examples/NMFS/afsc_walleye_pollock/run_pollock_random_effect_scaling.sh
```

That sweep runs six independently compiled configurations (`0, 1, 2, 5, 10,
20` recruitment effects), replaces the Pollock diagnostic files between cases,
and writes `random_effect_scaling_summary.csv` plus one log per case. It is a
scaling/identifiability diagnostic, not a strict optimizer acceptance test.

## Expected success

The default run exits zero, reports `converged yes`, and lists the generated
CSV files. A compiler error involving Eigen usually means the repository was
not cloned with its vendored dependencies intact or an incompatible compiler
was forced through `CXX`.
