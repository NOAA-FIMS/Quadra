# Level 3: Fleet-Specific Selectivity and Vulnerable-Biomass Indices

Level 2 added fleet-specific composition data. It improved optimization, but the
primary weak fixed-effect direction still involved abundance scale, fishing
mortality, and both q parameters.

The likely issue is that fleet-specific composition information was not fully
connected to the fleet-specific index likelihood. The Level 2 index model still
treated both fleet indices as observations of the same total biomass trajectory.

## Model change

Compared with Level 2:

- keep one shared biological population
- keep fleet-specific q
- add fleet-specific selectivity
- model each fleet index as an observation of fleet-vulnerable biomass

## Scientific hypothesis

Fleet-specific composition should identify fleet-specific selectivity. Once
indices observe fleet-vulnerable biomass, that selectivity information should
flow into the index likelihood and may reduce the abundance/q/F confounding.
