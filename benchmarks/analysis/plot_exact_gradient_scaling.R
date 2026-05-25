#!/usr/bin/env Rscript

csv_path <- "benchmarks/exact_laplace_gradient/state_space_exact_gradient_benchmark.csv"

if (!file.exists(csv_path)) {
  stop("Exact gradient scaling CSV not found.")
}

df <- read.csv(csv_path)

required <- c(
  "n_state",
  "objective_ms",
  "tape_setup_ms",
  "reverse_pass_ms",
  "gradient_extract_ms",
  "total_gradient_ms",
  "total_ms"
)

missing <- setdiff(required, names(df))

if (length(missing) > 0) {
  stop(paste("Missing columns:", paste(missing, collapse = ", ")))
}

png(
  filename = "benchmarks/analysis/exact_gradient_total_scaling.png",
  width = 1000,
  height = 700
)

plot(
  df$n_state,
  df$total_ms,
  type = "b",
  pch = 1,
  lty = 1,
  xlab = "Number of latent states",
  ylab = "Total exact-gradient time (ms)",
  main = "Exact Laplace Gradient Scaling"
)

dev.off()

png(
  filename = "benchmarks/analysis/exact_gradient_component_scaling.png",
  width = 1000,
  height = 700
)

component_matrix <- rbind(
  df$objective_ms,
  df$tape_setup_ms,
  df$reverse_pass_ms,
  df$gradient_extract_ms
)

matplot(
  df$n_state,
  t(component_matrix),
  type = "b",
  pch = c(1,2,3,4),
  lty = c(1,2,3,4),
  xlab = "Number of latent states",
  ylab = "Component time (ms)",
  main = "Exact Gradient Component Scaling"
)

legend(
  "topleft",
  legend = c(
    "objective",
    "tape_setup",
    "reverse_pass",
    "gradient_extract"
  ),
  pch = c(1,2,3,4),
  lty = c(1,2,3,4)
)

dev.off()

png(
  filename = "benchmarks/analysis/exact_gradient_stacked_scaling.png",
  width = 1000,
  height = 700
)

stack_matrix <- rbind(
  df$objective_ms,
  df$tape_setup_ms,
  df$reverse_pass_ms,
  df$gradient_extract_ms
)

barplot(
  stack_matrix,
  beside = FALSE,
  names.arg = df$n_state,
  xlab = "Number of latent states",
  ylab = "Time (ms)",
  main = "Stacked Exact Gradient Timing Decomposition"
)

legend(
  "topleft",
  legend = c(
    "objective",
    "tape_setup",
    "reverse_pass",
    "gradient_extract"
  ),
  fill = seq_len(4)
)

dev.off()

png(
  filename = "benchmarks/analysis/exact_gradient_reverse_vs_structure.png",
  width = 1000,
  height = 700
)

plot(
  df$hessian_nnz,
  df$reverse_pass_ms,
  type = "b",
  pch = 1,
  lty = 1,
  xlab = "nnz(H_uu)",
  ylab = "Reverse pass time (ms)",
  main = "Reverse Pass Cost vs Hessian Structure"
)

dev.off()

cat("Wrote exact gradient scaling plots.\n")
