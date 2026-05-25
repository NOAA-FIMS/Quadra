#!/usr/bin/env Rscript

csv_path <- "benchmarks/normalized/random_intercept_normalized.csv"

if (!file.exists(csv_path)) {
  stop("Normalized random-intercept benchmark CSV not found.")
}

df <- read.csv(csv_path)

if (nrow(df) == 0) {
  stop("Normalized random-intercept benchmark CSV is empty.")
}

plot_metric <- function(metric, ylab, title, output) {
  if (!(metric %in% names(df))) {
    stop(paste("Missing metric:", metric))
  }

  sub <- df[!is.na(df[[metric]]) & df[[metric]] != "", ]

  if (nrow(sub) == 0) {
    warning(paste("No data for metric:", metric))
    return(invisible(FALSE))
  }

  sub[[metric]] <- as.numeric(sub[[metric]])

  png(
    filename = output,
    width = 1000,
    height = 700
  )

  plot(
    sub$n_obs,
    sub[[metric]],
    type = "n",
    xlab = "Number of observations",
    ylab = ylab,
    main = title
  )

  engines <- unique(sub$engine)

  for (eng in engines) {
    idx <- sub$engine == eng
    ord <- order(sub$n_obs[idx])

    lines(
      sub$n_obs[idx][ord],
      sub[[metric]][idx][ord],
      type = "b"
    )
  }

  legend(
    "topleft",
    legend = engines,
    lty = 1,
    pch = 1
  )

  dev.off()

  invisible(TRUE)
}

plot_metric(
  "objective_eval_ms",
  "Objective evaluation time (ms)",
  "Random Intercept Objective Evaluation Scaling",
  "benchmarks/analysis/random_intercept_objective_eval_comparison.png"
)

plot_metric(
  "gradient_eval_ms",
  "Gradient evaluation time (ms)",
  "Random Intercept Gradient Evaluation Scaling",
  "benchmarks/analysis/random_intercept_gradient_eval_comparison.png"
)

plot_metric(
  "workspace_ms",
  "Workspace/setup time (ms)",
  "Random Intercept Workspace/Setup Scaling",
  "benchmarks/analysis/random_intercept_workspace_setup_comparison.png"
)

plot_metric(
  "total_wall_ms",
  "Total wall time (ms)",
  "Random Intercept Total Wall-Time Scaling",
  "benchmarks/analysis/random_intercept_total_wall_comparison.png"
)

cat("Wrote random-intercept comparative plots.\n")
