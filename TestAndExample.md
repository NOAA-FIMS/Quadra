
# Quadra Test and Example Scaffold

This scaffold adds a lightweight testing and example layout for Quadra.

## Layout

```text
tests/
  test_curvature_depends_on_theta.cpp
  test_poisson_random_effect.cpp
  test_ar1_random_walk.cpp
  test_hdot_validation.cpp

examples/
  fisheries_random_year_effects.cpp
  fisheries_age_selectivity_random_walk.cpp
  fisheries_index_cpue_laplace.cpp

docs/
  TESTING_PLAN.md
```

## Build examples

From the Quadra repository root, copy `tests/` and `examples/` into your project,
then compile a test with something like:

```bash
clang++ -std=c++17 -O3 -flto \
  -Icore/eigen \
  -o test_curvature tests/test_curvature_depends_on_theta.cpp
```

If your headers live in `core/`, the tests assume:

```cpp
#include "../core/optimizer.hpp"
#include "../core/laplace.hpp"
#include "../model/parameter.hpp"
```

Adjust include paths if your local layout differs.

## Validation build

To validate exact `Hdot` against the finite-difference fallback on small models:

```bash
clang++ -std=c++17 -O3 -DQUADRA_VALIDATE_HDOT \
  -Icore/eigen \
  -o test_hdot tests/test_hdot_validation.cpp
```

Do not use `QUADRA_VALIDATE_HDOT` for large random-effect tests.
