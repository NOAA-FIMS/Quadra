#!/usr/bin/env Rscript

csv_path <- "benchmarks/exact_laplace_gradient/exact_gradient_reuse_benchmark.csv"

if (!file.exists(csv_path)) {
  stop("Exact gradient reuse benchmark CSV not found.")
}

df <- read.csv(csv_path)

png(
  filename = "benchmarks/analysis/exact_gradient_reuse_total_ms.png",
  width = 1000,
  height = 700
)

plot(
  df$iteration,
  df$total_ms,
  type = "b",
  pch = 1,
  lty = 1,
  xlab = "Iteration",
  ylab = "Total exact-gradient time (ms)",
  main = "Exact Laplace Gradient Reuse Benchmark"
)

dev.off()

png(
  filename = "benchmarks/analysis/exact_gradient_reuse_components.png",
  width = 1000,
  height = 700
)

ylim <- range(
  c(
    df$objective_ms,
    df$tape_setup_ms,
    df$reverse_pass_ms,
    df$gradient_extract_ms
  ),
  na.rm = TRUE
)

plot(
  df$iteration,
  df$objective_ms,
  type = "b",
  pch = 1,
  lty = 1,
  ylim = ylim,
  xlab = "Iteration",
  ylab = "Time (ms)",
  main = "Exact Laplace Gradient Component Timing"
)

lines(df$iteration, df$tape_setup_ms, type = "b", pch = 2, lty = 2)
lines(df$iteration, df$reverse_pass_ms, type = "b", pch = 3, lty = 3)
lines(df$iteration, df$gradient_extract_ms, type = "b", pch = 4, lty = 4)

legend(
  "topright",
  legend = c("objective", "tape_setup", "reverse_pass", "gradient_extract"),
  pch = c(1, 2, 3, 4),
  lty = c(1, 2, 3, 4)
)

dev.off()

cat("Wrote exact gradient reuse plots.\n")
