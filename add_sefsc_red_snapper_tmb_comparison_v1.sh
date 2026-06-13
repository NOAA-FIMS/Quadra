#!/usr/bin/env bash
set -euo pipefail

echo "== Add SEFSC red snapper Quadra vs TMB comparison scaffold =="

mkdir -p examples/NMFS/sefsc_red_snapper/tmb
mkdir -p examples/NMFS/sefsc_red_snapper/outputs

cat > examples/NMFS/sefsc_red_snapper/tmb/red_snapper_tmb.cpp <<'CPP'
#include <TMB.hpp>

template <class Type>
Type square(Type x) { return x * x; }

template <class Type>
Type invlogit(Type x) {
  return Type(1.0) / (Type(1.0) + exp(-x));
}

template <class Type>
Type logistic_selectivity(Type age, Type a50, Type slope) {
  return Type(1.0) / (Type(1.0) + exp(-slope * (age - a50)));
}

template <class Type>
objective_function<Type>::operator()() {
  DATA_VECTOR(catch_obs);
  DATA_VECTOR(index_obs);
  DATA_MATRIX(age_comp_obs);

  PARAMETER(log_r0);
  PARAMETER(log_fbar);
  PARAMETER(log_q);
  PARAMETER(logit_sel_a50);
  PARAMETER(log_sel_slope);
  PARAMETER_VECTOR(log_rec_dev);

  const int n_years = catch_obs.size();
  const int n_ages = age_comp_obs.cols();

  Type r0 = exp(log_r0);
  Type m = Type(0.18);
  Type fbar = exp(log_fbar);
  Type q = exp(log_q);
  Type sel_a50 = Type(1.0) + Type(9.0) * invlogit(logit_sel_a50);
  Type sel_slope = exp(log_sel_slope);

  Type sigma_log_index = Type(0.20);
  Type sigma_log_catch = Type(0.15);
  Type sigma_rec_dev = Type(0.35);
  Type age_comp_effective_n = Type(2.0);
  Type min_positive = Type(1.0e-12);

  vector<Type> weight(n_ages);
  for (int a = 0; a < n_ages; ++a) {
    weight(a) = Type(0.35) * pow(Type(a + 1), Type(2.8));
  }

  vector<Type> sel(n_ages);
  for (int a = 0; a < n_ages; ++a) {
    sel(a) = logistic_selectivity(Type(a + 1), sel_a50, sel_slope);
  }

  vector<Type> n(n_ages);
  n(0) = r0;
  for (int a = 1; a < n_ages; ++a) {
    n(a) = n(a - 1) * exp(-m);
  }
  n(n_ages - 1) = n(n_ages - 1) / (Type(1.0) - exp(-m));

  Type nll = Type(0.0);

  nll += Type(0.5) * square((log_r0 - Type(std::log(1200.0))) / Type(1.0));
  nll += Type(0.5) * square((log_fbar - Type(std::log(0.025))) / Type(0.75));
  nll += Type(0.5) * square((log_q - Type(std::log(0.00005))) / Type(1.0));
  nll += Type(0.5) * square((sel_a50 - Type(4.0)) / Type(0.75));
  nll += Type(0.5) * square((log_sel_slope - Type(std::log(1.2))) / Type(0.35));

  for (int y = 0; y < n_years; ++y) {
    Type rec_dev = log_rec_dev(y);
    nll += Type(0.5) * square(rec_dev / sigma_rec_dev);

    Type total_biomass = Type(0.0);
    Type catch_hat = Type(0.0);
    Type selected_sum = Type(0.0);
    vector<Type> pred_age_comp(n_ages);

    for (int a = 0; a < n_ages; ++a) {
      Type fa = fbar * sel(a);
      Type za = m + fa;
      total_biomass += n(a) * weight(a);
      catch_hat += n(a) * weight(a) * fa / za * (Type(1.0) - exp(-za));
      pred_age_comp(a) = n(a) * sel(a);
      selected_sum += pred_age_comp(a);
    }

    Type index_hat = q * total_biomass;

    if (index_obs(y) > 0.0) {
      Type z = (log(Type(index_obs(y))) - log(CppAD::CondExpGt(index_hat, min_positive, index_hat, min_positive))) / sigma_log_index;
      nll += Type(0.5) * z * z;
    }

    if (catch_obs(y) > 0.0) {
      Type z = (log(Type(catch_obs(y))) - log(CppAD::CondExpGt(catch_hat, min_positive, catch_hat, min_positive))) / sigma_log_catch;
      nll += Type(0.5) * z * z;
    }

    for (int a = 0; a < n_ages; ++a) {
      pred_age_comp(a) = pred_age_comp(a) / CppAD::CondExpGt(selected_sum, min_positive, selected_sum, min_positive);
      Type obs = age_comp_obs(y, a);
      if (obs > 0.0) {
        nll -= age_comp_effective_n * obs * log(CppAD::CondExpGt(pred_age_comp(a), min_positive, pred_age_comp(a), min_positive));
      }
    }

    vector<Type> next(n_ages);
    next.setZero();
    next(0) = r0 * exp(rec_dev);

    for (int a = 1; a < n_ages; ++a) {
      Type f_prev = fbar * sel(a - 1);
      Type z_prev = m + f_prev;
      next(a) = n(a - 1) * exp(-z_prev);
    }

    Type f_last = fbar * sel(n_ages - 1);
    Type z_last = m + f_last;
    next(n_ages - 1) += n(n_ages - 1) * exp(-z_last);

    n = next;
  }

  return nll;
}
CPP

cat > examples/NMFS/sefsc_red_snapper/tmb/run_red_snapper_tmb_fit.R <<'RSCRIPT'
#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(TMB))

data_candidates <- c(
  "examples/NMFS/sefsc_red_snapper/data/red_snapper_synthetic_observations.csv",
  "examples/NMFS/sefsc_red_snapper/data/synthetic_red_snapper_observations.csv",
  "examples/NMFS/sefsc_red_snapper/data/red_snapper_observations.csv"
)

data_path <- data_candidates[file.exists(data_candidates)][1]
if (is.na(data_path)) {
  stop("Could not find SEFSC red snapper observation CSV")
}

obs <- read.csv(data_path)

catch_col <- grep("catch", names(obs), value = TRUE)[1]
index_col <- grep("index", names(obs), value = TRUE)[1]
age_cols <- grep("^age_|^age[0-9]+|comp", names(obs), value = TRUE)
age_cols <- age_cols[sapply(obs[age_cols], is.numeric)]

if (is.na(catch_col) || is.na(index_col) || length(age_cols) == 0) {
  stop("Could not infer catch/index/age-composition columns from data CSV")
}

catch_obs <- as.numeric(obs[[catch_col]])
index_obs <- as.numeric(obs[[index_col]])
age_comp_obs <- as.matrix(obs[, age_cols, drop = FALSE])
age_comp_obs <- age_comp_obs / rowSums(age_comp_obs)

cpp <- "examples/NMFS/sefsc_red_snapper/tmb/red_snapper_tmb.cpp"
TMB::compile(cpp)
dyn.load(TMB::dynlib(sub("\\.cpp$", "", cpp)))

parameters <- list(
  log_r0 = log(1200.0),
  log_fbar = log(0.025),
  log_q = log(0.00005),
  logit_sel_a50 = 0.0,
  log_sel_slope = log(1.2),
  log_rec_dev = rep(0.0, length(catch_obs))
)

obj <- MakeADFun(
  data = list(catch_obs = catch_obs, index_obs = index_obs, age_comp_obs = age_comp_obs),
  parameters = parameters,
  random = "log_rec_dev",
  DLL = "red_snapper_tmb",
  silent = TRUE
)

fit <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1000, iter.max = 1000))
pl <- obj$env$parList(obj$env$last.par.best)

summary_path <- "examples/NMFS/sefsc_red_snapper/outputs/tmb_fit_summary.csv"
out <- data.frame(
  field = c("objective", "convergence", "message", "log_r0", "r0",
            "log_fbar", "fbar", "log_q", "q", "logit_sel_a50",
            "sel_a50", "log_sel_slope", "sel_slope", "random_effects"),
  value = c(fit$objective, fit$convergence, fit$message,
            pl$log_r0, exp(pl$log_r0),
            pl$log_fbar, exp(pl$log_fbar),
            pl$log_q, exp(pl$log_q),
            pl$logit_sel_a50,
            1.0 + 9.0 / (1.0 + exp(-pl$logit_sel_a50)),
            pl$log_sel_slope, exp(pl$log_sel_slope),
            length(pl$log_rec_dev))
)
write.csv(out, summary_path, row.names = FALSE, quote = FALSE)

rec_path <- "examples/NMFS/sefsc_red_snapper/outputs/tmb_recruitment_deviations.csv"
write.csv(data.frame(year = seq_along(pl$log_rec_dev),
                     log_rec_dev = as.numeric(pl$log_rec_dev),
                     rec_multiplier = exp(as.numeric(pl$log_rec_dev))),
          rec_path, row.names = FALSE, quote = FALSE)

cat("wrote:", summary_path, "\n")
cat("wrote:", rec_path, "\n")
RSCRIPT
chmod +x examples/NMFS/sefsc_red_snapper/tmb/run_red_snapper_tmb_fit.R

cat > examples/NMFS/sefsc_red_snapper/compare_quadra_tmb_fit.py <<'PY'
#!/usr/bin/env python3
from pathlib import Path
import csv
import math

out = Path("examples/NMFS/sefsc_red_snapper/outputs")

def read_summary(path):
    d = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            try:
                d[row["field"]] = float(row["value"])
            except Exception:
                d[row["field"]] = row["value"]
    return d

q = read_summary(out / "quadra_fit_summary.csv")
t = read_summary(out / "tmb_fit_summary.csv")

fields = ["objective", "r0", "fbar", "q", "sel_a50", "sel_slope", "random_effects"]
path = out / "quadra_vs_tmb_fit_comparison.csv"

with open(path, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["field", "quadra", "tmb", "difference", "relative_difference"])
    for field in fields:
        qv = q.get(field, "")
        tv = t.get(field, "")
        diff = ""
        rel = ""
        if isinstance(qv, float) and isinstance(tv, float):
            diff = qv - tv
            rel = diff / tv if tv != 0 and math.isfinite(tv) else ""
        w.writerow([field, qv, tv, diff, rel])

print(f"wrote: {path}")
PY
chmod +x examples/NMFS/sefsc_red_snapper/compare_quadra_tmb_fit.py

cat > examples/NMFS/sefsc_red_snapper/run_quadra_vs_tmb_comparison.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

./examples/NMFS/sefsc_red_snapper/run_red_snapper_quadra_fit.sh
Rscript examples/NMFS/sefsc_red_snapper/tmb/run_red_snapper_tmb_fit.R
python3 examples/NMFS/sefsc_red_snapper/compare_quadra_tmb_fit.py
cat examples/NMFS/sefsc_red_snapper/outputs/quadra_vs_tmb_fit_comparison.csv
SH
chmod +x examples/NMFS/sefsc_red_snapper/run_quadra_vs_tmb_comparison.sh

echo
echo "Installed Quadra vs TMB comparison scaffold."
echo
echo "Run:"
echo "  ./examples/NMFS/sefsc_red_snapper/run_quadra_vs_tmb_comparison.sh"
