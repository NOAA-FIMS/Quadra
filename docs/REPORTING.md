# Explicit report registry

Quadra uses an explicit report registry instead of TMB-style macros.

## Model pattern

```cpp
template <typename T>
T objective(const std::vector<T>& p) const;

template <typename T, typename ReportLike>
void report(const std::vector<T>& p, ReportLike& out) const;
```

Optimization calls only the objective. Reporting is evaluated post-fit.

## Report API

```cpp
out.add("theta", theta);          // value only
out.estimate("exp_theta", exp(theta)); // value + delta-method SE
```

## Post-fit use

```cpp
auto fit = optimize_lbfgs(model, params, lopts);
auto cov = estimate_fixed_covariance(model, params, fit, lopts, copts);
auto report = evaluate_report_with_uncertainty(model, params, cov);
report.print();
```
