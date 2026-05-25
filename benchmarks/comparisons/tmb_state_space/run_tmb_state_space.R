#!/usr/bin/env Rscript

if (!requireNamespace("TMB", quietly = TRUE)) {
  stop("TMB is not installed. Install TMB before running this benchmark.")
}

library(TMB)
library(Matrix)

out_dir <- "benchmarks/comparisons/tmb_state_space/comparison_outputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

template <- "benchmarks/comparisons/tmb_state_space/state_space_tmb.cpp"

TMB::compile(
  template,
  flags = "-O3 -w -DTMB_EIGEN_DISABLE_WARNINGS"
)
dyn.load(TMB::dynlib(sub("\\.cpp$", "", template)))

elapsed_ms <- function(expr) {
  start <- proc.time()[["elapsed"]]
  value <- force(expr)
  end <- proc.time()[["elapsed"]]
  list(value = value, ms = 1000.0 * (end - start))
}

run_one <- function(n_state) {
  data <- list(n_state = as.integer(n_state))

  parameters <- list(
    mu = 1.0,
    log_sigma_obs = log(0.4),
    log_sigma_rw = log(0.3),
    x = rep(0.0, n_state)
  )

  setup <- elapsed_ms({
    TMB::MakeADFun(
      data = data,
      parameters = parameters,
      random = "x",
      DLL = "state_space_tmb",
      silent = TRUE
    )
  })

  obj <- setup$value

  objective_reps <- if (n_state <= 50) 1000L else 200L

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

  fn_gr_time <- elapsed_ms({
    sink <- 0.0
    for (i in seq_len(objective_reps)) {
      sink <- sink + obj$fn(obj$par) + sum(obj$gr(obj$par))
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

  hessian_nnz <- NA_integer_
  hessian_density <- NA_real_
  factor_nnz <- NA_integer_
  fill_ratio <- NA_real_

  structure_time <- elapsed_ms({
    H <- NULL

    if (!is.null(obj$env$spHess)) {
      H <- tryCatch(
        obj$env$spHess(random = TRUE),
        error = function(e) NULL
      )
    }

    if (is.null(H) && !is.null(obj$env$hessian)) {
      H <- tryCatch(
        obj$env$hessian(random = TRUE),
        error = function(e) NULL
      )
    }

    if (!is.null(H)) {
      hessian_nnz <<- length(H@x)
      hessian_density <<- hessian_nnz / (n_state * n_state)

      chol_H <- tryCatch(
        Matrix::Cholesky(H, LDL = TRUE, perm = TRUE),
        error = function(e) NULL
      )

      if (!is.null(chol_H)) {
        L <- tryCatch(
          Matrix::expand(chol_H)$L,
          error = function(e) NULL
        )

        if (!is.null(L)) {
          factor_nnz <<- length(L@x)

          if (!is.na(hessian_nnz) && hessian_nnz > 0) {
            fill_ratio <<- factor_nnz / hessian_nnz
          }
        }
      }
    }
  })

  data.frame(
    engine = "tmb",
    model = "state_space",
    n_state = n_state,
    n_random = n_state,
    setup_ms = setup$ms,
    objective_eval_ms = objective_time$ms / objective_reps,
    gradient_eval_ms = gradient_time$ms / objective_reps,
    fn_gr_eval_ms = fn_gr_time$ms / objective_reps,
    objective_reps = objective_reps,
    optimization_ms = opt_time$ms,
    structure_ms = structure_time$ms,
    total_wall_ms = setup$ms + objective_time$ms + gradient_time$ms + opt_time$ms + structure_time$ms,
    hessian_nnz = hessian_nnz,
    hessian_density = hessian_density,
    factor_nnz = factor_nnz,
    fill_ratio = fill_ratio,
    objective = opt$objective,
    success = as.integer(opt$convergence == 0L)
  )
}

results <- do.call(
  rbind,
  lapply(c(25, 50, 100, 250), run_one)
)

write.csv(
  results,
  file = file.path(out_dir, "tmb_state_space_compare.csv"),
  row.names = FALSE
)

print(results)
