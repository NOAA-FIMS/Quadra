#!/usr/bin/env Rscript

csv_path <- "benchmarks/normalized/random_intercept_normalized.csv"

if (!file.exists(csv_path)) {
  stop("Normalized benchmark CSV not found.")
}

df <- read.csv(csv_path)

if (nrow(df) == 0) {
  stop("Normalized benchmark CSV is empty.")
}

png(
  filename = "benchmarks/analysis/random_intercept_scaling.png",
  width = 1000,
  height = 700
)

plot(
  df$n_obs,
  df$total_wall_ms,
  type = "n",
  xlab = "Number of observations",
  ylab = "Total wall time (ms)",
  main = "Quadra vs TMB Random Intercept Scaling"
)

engines <- unique(df$engine)

pch_values <- c(1, 2, 3, 4, 5, 6)
lty_values <- c(1, 2, 3, 4, 5, 6)

for (i in seq_along(engines)) {
  eng <- engines[i]
  idx <- df$engine == eng
  ord <- order(df$n_obs[idx])

  lines(
    df$n_obs[idx][ord],
    df$total_wall_ms[idx][ord],
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

cat("Wrote scaling plot PNG.\n")
