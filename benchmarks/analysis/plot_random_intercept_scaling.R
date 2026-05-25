#!/usr/bin/env Rscript

csv_path <- "benchmarks/normalized/random_intercept_normalized.csv"

if (!file.exists(csv_path)) {
  stop("Normalized benchmark CSV not found.")
}

df <- read.csv(csv_path)

png(
  filename = "benchmarks/analysis/random_intercept_scaling.png",
  width = 900,
  height = 600
)

plot(
  df$n_obs,
  df$total_wall_ms,
  type = "n",
  xlab = "Number of observations",
  ylab = "Total wall time (ms)",
  main = "Random Intercept Scaling"
)

engines <- unique(df$engine)

for (eng in engines) {
  idx <- df$engine == eng

  lines(
    df$n_obs[idx],
    df$total_wall_ms[idx],
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

cat("Wrote scaling plot.\n")
