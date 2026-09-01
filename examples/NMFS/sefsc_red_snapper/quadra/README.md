# Quadra Implementation

This directory contains the Quadra drivers and age-structured model code. Use
the shell runners in the parent directory instead of compiling these files by
hand; the runners supply the required include paths and AD-graph translation
unit and place binaries under the repository's `build/examples/` directory.

The primary entry point is `../run_red_snapper_quadra_fit.sh`. The Level-0,
objective, and deterministic age-structured runners provide shorter diagnostic
stages.
