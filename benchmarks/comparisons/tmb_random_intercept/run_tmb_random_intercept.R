#!/usr/bin/env Rscript

if (!requireNamespace("TMB", quietly = TRUE)) {
  stop("TMB is not installed. Install TMB before running this benchmark.")
}

library(TMB)

out_dir <- "benchmarks/comparisons/tmb_random_intercept/comparison_outputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

template <- "benchmarks/comparisons/tmb_random_intercept/random_intercept_tmb.cpp"

TMB::compile(template)
dyn.load(TMB::dynlib(sub("\\.cpp$", "", template)))

run_one <- function(n_obs) {
  data <- list(n_obs = as.integer(n_obs))

  parameters <- list(
    mu = 0.0,
    log_sigma = log(0.5),
    u = 0.0
  )

  random <- "u"

  start <- proc.time()[["elapsed"]]
  obj <- TMB::MakeADFun(
    data = data,
    parameters = parameters,
    random = random,
    DLL = "random_intercept_tmb",
    silent = TRUE
  )

  opt <- stats::nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr
  )

  elapsed_ms <- 1000.0 * (proc.time()[["elapsed"]] - start)

  data.frame(
    engine = "tmb",
    n_obs = n_obs,
    elapsed_ms = elapsed_ms,
    objective = opt$objective,
    convergence = opt$convergence
  )
}

results <- do.call(
  rbind,
  lapply(c(100, 1000, 5000, 10000), run_one)
)

write.csv(
  results,
  file = file.path(out_dir, "tmb_random_intercept_compare.csv"),
  row.names = FALSE
)

print(results)
