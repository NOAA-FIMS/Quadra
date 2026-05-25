#!/usr/bin/env Rscript

out <- "benchmarks/analysis/quadra_benchmark_plots.pdf"

safe_read <- function(path) {
  if (!file.exists(path)) return(NULL)
  read.csv(path)
}

plot_metric <- function(df, x, y, title, xlab, ylab) {
  if (is.null(df) || !(x %in% names(df)) || !(y %in% names(df))) {
    plot.new()
    title(main = title)
    text(0.5, 0.5, "No data available")
    return()
  }

  yy <- suppressWarnings(as.numeric(df[[y]]))
  xx <- suppressWarnings(as.numeric(df[[x]]))
  ok <- is.finite(xx) & is.finite(yy)

  if (!any(ok)) {
    plot.new()
    title(main = title)
    text(0.5, 0.5, "No finite data available")
    return()
  }

  plot(
    xx[ok],
    yy[ok],
    type = "b",
    pch = 1,
    lty = 1,
    xlab = xlab,
    ylab = ylab,
    main = title
  )
}

random <- safe_read("benchmarks/normalized/random_intercept_normalized.csv")
state <- safe_read("benchmarks/normalized/state_space_normalized.csv")
exact <- safe_read("benchmarks/exact_laplace_gradient/state_space_exact_gradient_benchmark.csv")
reuse <- safe_read("benchmarks/exact_laplace_gradient/exact_gradient_reuse_benchmark.csv")
factor <- safe_read("benchmarks/exact_laplace_gradient/factorization_reuse_benchmark.csv")

pdf(out, width = 11, height = 8.5)

plot_metric(
  random,
  "n_obs",
  "total_wall_ms",
  "Random Intercept Total Wall Time",
  "Number of observations",
  "Total wall time (ms)"
)

plot_metric(
  random,
  "n_obs",
  "objective_eval_ms",
  "Random Intercept Objective Evaluation",
  "Number of observations",
  "Objective eval time (ms)"
)

plot_metric(
  state,
  "n_state",
  "total_wall_ms",
  "State-Space Total Wall Time",
  "Number of latent states",
  "Total wall time (ms)"
)

plot_metric(
  state,
  "n_state",
  "hessian_nnz",
  "State-Space Hessian Nonzeros",
  "Number of latent states",
  "nnz(H_uu)"
)

plot_metric(
  state,
  "n_state",
  "fill_ratio",
  "State-Space Fill Ratio",
  "Number of latent states",
  "factor nnz / Hessian nnz"
)

plot_metric(
  exact,
  "n_state",
  "total_ms",
  "Exact Laplace Gradient Scaling",
  "Number of latent states",
  "Total exact-gradient time (ms)"
)

plot_metric(
  exact,
  "hessian_nnz",
  "reverse_pass_ms",
  "Exact Gradient Reverse Pass vs Hessian Structure",
  "nnz(H_uu)",
  "Reverse pass time (ms)"
)

plot_metric(
  reuse,
  "iteration",
  "total_ms",
  "Exact Gradient Reuse Timing",
  "Iteration",
  "Total exact-gradient time (ms)"
)

plot_metric(
  factor,
  "n_state",
  "reuse_ratio",
  "Sparse Factorization Reuse Ratio",
  "Number of latent states",
  "reuse_ms / fresh_ms"
)

dev.off()

cat("Wrote", out, "\n")
