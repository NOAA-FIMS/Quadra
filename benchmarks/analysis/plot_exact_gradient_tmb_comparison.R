#!/usr/bin/env Rscript

csv_path <- "benchmarks/normalized/exact_gradient_state_space_comparison.csv"

if (!file.exists(csv_path)) {
  stop("Exact-gradient comparison CSV not found.")
}

df <- read.csv(csv_path)

if (nrow(df) == 0) {
  stop("Exact-gradient comparison CSV is empty.")
}

png(
  filename = "benchmarks/analysis/exact_gradient_tmb_state_space_comparison.png",
  width = 1000,
  height = 700
)

plot(
  df$n_state,
  rep(NA_real_, nrow(df)),
  type = "n",
  xlab = "Number of latent states",
  ylab = "Time (ms)",
  main = "State-Space Exact Gradient / fn+gr Comparison"
)

legend_entries <- c()
pch_values <- c()
lty_values <- c()

q <- df[df$engine == "quadra" & !is.na(df$exact_gradient_ms), ]

if (nrow(q) > 0) {
  q$exact_gradient_ms <- as.numeric(q$exact_gradient_ms)
  ord <- order(q$n_state)
  lines(q$n_state[ord], q$exact_gradient_ms[ord], type = "b", pch = 1, lty = 1)
  legend_entries <- c(legend_entries, "Quadra exact gradient")
  pch_values <- c(pch_values, 1)
  lty_values <- c(lty_values, 1)
}

t <- df[df$engine == "tmb" & !is.na(df$fn_gr_eval_ms), ]

if (nrow(t) > 0) {
  t$fn_gr_eval_ms <- as.numeric(t$fn_gr_eval_ms)
  ord <- order(t$n_state)
  lines(t$n_state[ord], t$fn_gr_eval_ms[ord], type = "b", pch = 2, lty = 2)
  legend_entries <- c(legend_entries, "TMB fn+gr eval")
  pch_values <- c(pch_values, 2)
  lty_values <- c(lty_values, 2)
}

if (length(legend_entries) > 0) {
  legend(
    "topleft",
    legend = legend_entries,
    pch = pch_values,
    lty = lty_values
  )
}

dev.off()

cat("Wrote exact-gradient TMB comparison plot.\n")
