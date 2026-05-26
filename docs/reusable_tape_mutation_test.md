# Reusable tape mutation behavior test

This test checks whether Quadra/had_quadra supports reusing an existing graph
after mutating independent variable values in-place.

The test intentionally:

1. Builds the graph once
2. Constructs:
   ```text
   y = x^2
   ```
3. Runs reverse mode
4. Mutates:
   ```cpp
   x.val = 3.0;
   ```
5. Runs reverse mode again WITHOUT rebuilding the graph

Expected results:

| x | y | dy/dx |
|---|---|---|
| 2 | 4 | 4 |
| 3 | 9 | 6 |

If the second pass updates correctly, Quadra can likely support reusable tapes
without rebuilding graphs every optimizer iteration.
