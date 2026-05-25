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

elapsed_ms <- function(expr) {
  start <- proc.time()[["elapsed"]]
  value <- force(expr)
  end <- proc.time()[["elapsed"]]
  list(value = value, ms = 1000.0 * (end - start))
}

run_one <- function(n_obs) {
  data <- list(n_obs = as.integer(n_obs))

  parameters <- list(
    mu = 0.0,
    log_sigma = log(0.5),
    u = 0.0
  )

  setup <- elapsed_ms({
    TMB::MakeADFun(
      data = data,
      parameters = parameters,
      random = "u",
      DLL = "random_intercept_tmb",
      silent = TRUE
    )
  })

  obj <- setup$value
  objective_reps <- 1000L

  objective_time <- elapsed_ms({
    sink <- 0.0
    for (i in seq_len(objective_reps)) {
      sink <- sink + obj$fn(obj$par)
    }
    sink
  })

  gradient_time <- elapsed_ms({
    sink <- 0.0
    for (i in seq_len(objective_reps)) {
      sink <- sink + sum(obj$gr(obj$par))
    }
    sink
  })

  opt_time <- elapsed_ms({
    stats::nlminb(
      start = obj$par,
      objective = obj$fn,
      gradient = obj$gr
    )
  })

  opt <- opt_time$value

  data.frame(
    engine = "tmb",
    n_obs = n_obs,
    n_random = 1L,
    setup_ms = setup$ms,
    objective_eval_ms = objective_time$ms / objective_reps,
    gradient_eval_ms = gradient_time$ms / objective_reps,
    objective_reps = objective_reps,
    optimization_ms = opt_time$ms,
    total_wall_ms = setup$ms + objective_time$ms + gradient_time$ms + opt_time$ms,
    objective = opt$objective,
    success = as.integer(opt$convergence == 0L)
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
