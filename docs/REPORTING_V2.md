# Report registry v2

Quadra's report registry is explicit and macro-free.

## Hierarchical names

```cpp
out.estimate("biomass/SSB_2025", ssb);
out.estimate("reference_points/MSY", msy);
out.add("diagnostics/max_gradient", max_grad);
```

The group is inferred from the path before the last slash.

## Metadata

```cpp
quadra::ReportMetadata meta;
meta.units = "metric tons";
meta.description = "Spawning stock biomass";

out.estimate("biomass/SSB", ssb, meta);
```

## CSV export

```cpp
auto report = evaluate_report_with_uncertainty(model, params, cov);
report.to_csv("report.csv");
```

CSV columns:

```text
path,group,name,value,std_error,estimate_uncertainty,units,description
```

## Design principle

Optimization does not touch reporting. Reporting is evaluated post-fit only.
