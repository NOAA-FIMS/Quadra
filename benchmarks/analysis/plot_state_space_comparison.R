#!/usr/bin/env Rscript

csv_path <- "benchmarks/normalized/state_space_normalized.csv"

if (!file.exists(csv_path)) {
  stop("Normalized state-space benchmark CSV not found.")
}

df <- read.csv(csv_path)

if (nrow(df) == 0) {
  stop("Normalized state-space benchmark CSV is empty.")
}

plot_metric <- function(metric, ylab, title, output) {
  if (!(metric %in% names(df))) {
    stop(paste("Missing metric:", metric))
  }

  sub <- df[!is.na(df[[metric]]) & df[[metric]] != "", ]

  if (nrow(sub) == 0) {
    warning(paste("No data for metric:", metric))

    png(filename = output, width = 1000, height = 700)
    plot.new()
    title(main = title)
    text(0.5, 0.5, paste("No data available for", metric))
    dev.off()

    return(invisible(FALSE))
  }

  sub[[metric]] <- as.numeric(sub[[metric]])

  png(filename = output, width = 1000, height = 700)

  plot(
    sub$n_state,
    sub[[metric]],
    type = "n",
    xlab = "Number of latent states",
    ylab = ylab,
    main = title
  )

  engines <- unique(sub$engine)

  pch_values <- c(1, 2, 3, 4, 5, 6)
  lty_values <- c(1, 2, 3, 4, 5, 6)

  for (i in seq_along(engines)) {
    eng <- engines[i]
    idx <- sub$engine == eng
    ord <- order(sub$n_state[idx])

    lines(
      sub$n_state[idx][ord],
      sub[[metric]][idx][ord],
      type = "b",
      pch = pch_values[i],
      lty = lty_values[i]
    )
  }

  legend(
    "topleft",
    legend = engines,
    lty = lty_values[seq_along(engines)],
    pch = pch_values[seq_along(engines)]
  )

  dev.off()

  invisible(TRUE)
}

plot_metric(
  "total_wall_ms",
  "Total wall time (ms)",
  "State-Space Total Wall-Time Scaling",
  "benchmarks/analysis/state_space_total_wall_comparison.png"
)

plot_metric(
  "workspace_ms",
  "Workspace/setup time (ms)",
  "State-Space Workspace/Setup Scaling",
  "benchmarks/analysis/state_space_workspace_setup_comparison.png"
)

plot_metric(
  "factorization_ms",
  "Factorization time (ms)",
  "State-Space Factorization Scaling",
  "benchmarks/analysis/state_space_factorization_comparison.png"
)

plot_metric(
  "hessian_nnz",
  "nnz(H_uu)",
  "State-Space Hessian Nonzeros",
  "benchmarks/analysis/state_space_hessian_nnz_comparison.png"
)

plot_metric(
  "factor_nnz",
  "Factor nnz",
  "State-Space Factor Nonzeros",
  "benchmarks/analysis/state_space_factor_nnz_comparison.png"
)

plot_metric(
  "fill_ratio",
  "factor nnz / Hessian nnz",
  "State-Space Fill Ratio",
  "benchmarks/analysis/state_space_fill_ratio_comparison.png"
)

cat("Wrote state-space comparative plots.\n")
