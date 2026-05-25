#!/usr/bin/env Rscript

out <- "benchmarks/analysis/quadra_benchmark_plots.pdf"

plots <- c(
  "benchmarks/analysis/random_intercept_scaling.png",
  "benchmarks/analysis/random_intercept_objective_eval_comparison.png",
  "benchmarks/analysis/random_intercept_gradient_eval_comparison.png",
  "benchmarks/analysis/random_intercept_workspace_setup_comparison.png",
  "benchmarks/analysis/random_intercept_total_wall_comparison.png",
  "benchmarks/analysis/state_space_total_wall_comparison.png",
  "benchmarks/analysis/state_space_workspace_setup_comparison.png",
  "benchmarks/analysis/state_space_factorization_comparison.png",
  "benchmarks/analysis/state_space_hessian_nnz_comparison.png",
  "benchmarks/analysis/state_space_factor_nnz_comparison.png",
  "benchmarks/analysis/state_space_fill_ratio_comparison.png",
  "benchmarks/analysis/exact_gradient_total_scaling.png",
  "benchmarks/analysis/exact_gradient_component_scaling.png",
  "benchmarks/analysis/exact_gradient_stacked_scaling.png",
  "benchmarks/analysis/exact_gradient_reverse_vs_structure.png",
  "benchmarks/analysis/exact_gradient_reuse_total_ms.png",
  "benchmarks/analysis/exact_gradient_reuse_components.png",
  "benchmarks/analysis/factorization_reuse_ratio.png",
  "benchmarks/analysis/factorization_reuse_times.png",
  "benchmarks/analysis/factorization_reuse_vs_fill.png"
)

plots <- plots[file.exists(plots)]

if (length(plots) == 0) {
  stop("No benchmark plots found.")
}

pdf(out, width = 11, height = 8.5)

for (plot_path in plots) {
  img <- png::readPNG(plot_path)

  plot.new()
  rasterImage(img, 0, 0, 1, 1)
  title(main = basename(plot_path), line = -1)
}

dev.off()

cat("Wrote", out, "\n")
