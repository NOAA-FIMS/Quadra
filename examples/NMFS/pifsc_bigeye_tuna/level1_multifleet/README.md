# Level 1: Multi-Fleet Bigeye Tuna Prototype

This level adds multiple fleets to the Level 0 single-region model.

The first implementation is intentionally conservative. It introduces a
multi-fleet data file and aggregates fleets into a single biological time series
so the diagnostic pipeline remains runnable. The next step is to make the fleet
likelihood explicit.

## Purpose

Use the diagnostic suite to ask whether adding fleets improves the model or
introduces confounding.

## Initial diagnostic questions

- Which fleet appears to inform abundance?
- Does splitting data into fleets improve or degrade curvature?
- Does Huu remain positive definite?
- Does fleet structure create q/F/recruitment confounding?

## Status

Scaffolded. Fleet observations are aggregated for the first runnable baseline.
