#!/usr/bin/env Rscript

out <- "benchmarks/analysis/quadra_benchmark_report.pdf"

safe_read <- function(path) {
  if (!file.exists(path)) return(NULL)
  read.csv(path)
}

fmt <- function(x, digits = 3) {
  if (is.null(x) || length(x) == 0 || !is.finite(x)) return("NA")
  format(round(x, digits), nsmall = digits)
}

page_title <- function(title, subtitle = NULL) {
  plot.new()
  text(0.5, 0.72, title, cex = 2.2, font = 2)
  if (!is.null(subtitle)) {
    text(0.5, 0.58, subtitle, cex = 1.1)
  }
  text(
    0.5, 0.36,
    paste(
      "Generated from benchmark CSV outputs.",
      "Interpret results by phase; not all engine timings are equivalent.",
      sep = "\n"
    ),
    cex = 0.95
  )
}

page_text <- function(title, bullets) {
  plot.new()
  text(0.05, 0.92, title, adj = 0, cex = 1.6, font = 2)

  y <- 0.82
  for (b in bullets) {
    text(0.08, y, paste0("• ", b), adj = 0, cex = 0.95)
    y <- y - 0.075
  }
}

plot_by_engine <- function(df, x, y, title, xlab, ylab, log_y = FALSE) {
  if (is.null(df) || !(x %in% names(df)) || !(y %in% names(df))) {
    plot.new(); title(main = title); text(0.5, 0.5, "No data available"); return()
  }

  df[[x]] <- suppressWarnings(as.numeric(df[[x]]))
  df[[y]] <- suppressWarnings(as.numeric(df[[y]]))
  df <- df[is.finite(df[[x]]) & is.finite(df[[y]]), ]

  if (nrow(df) == 0) {
    plot.new(); title(main = title); text(0.5, 0.5, "No finite data available"); return()
  }

  yvals <- if (log_y) log10(df[[y]]) else df[[y]]
  ylab2 <- if (log_y) paste0(ylab, " log10 scale") else ylab

  plot(
    df[[x]], yvals,
    type = "n",
    xlab = xlab,
    ylab = ylab2,
    main = title
  )

  engines <- unique(df$engine)
  pch_values <- c(1, 2, 3, 4, 5, 6)
  lty_values <- c(1, 2, 3, 4, 5, 6)

  for (i in seq_along(engines)) {
    eng <- engines[i]
    sub <- df[df$engine == eng, ]
    ord <- order(sub[[x]])

    yy <- if (log_y) log10(sub[[y]][ord]) else sub[[y]][ord]

    lines(
      sub[[x]][ord],
      yy,
      type = "b",
      pch = pch_values[i],
      lty = lty_values[i]
    )
  }

  legend(
    "topleft",
    legend = engines,
    pch = pch_values[seq_along(engines)],
    lty = lty_values[seq_along(engines)],
    bty = "n"
  )
}

plot_metric <- function(df, x, y, title, xlab, ylab, caption = NULL, log_y = FALSE) {
  par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
  plot_by_engine(df, x, y, title, xlab, ylab, log_y = log_y)
  if (!is.null(caption)) {
    mtext(caption, side = 1, line = 4, cex = 0.75)
  }
}

plot_single <- function(df, x, y, title, xlab, ylab, caption = NULL) {
  if (is.null(df) || !(x %in% names(df)) || !(y %in% names(df))) {
    plot.new(); title(main = title); text(0.5, 0.5, "No data available"); return()
  }

  xx <- suppressWarnings(as.numeric(df[[x]]))
  yy <- suppressWarnings(as.numeric(df[[y]]))
  ok <- is.finite(xx) & is.finite(yy)

  if (!any(ok)) {
    plot.new(); title(main = title); text(0.5, 0.5, "No finite data available"); return()
  }

  plot(
    xx[ok], yy[ok],
    type = "b",
    pch = 1,
    lty = 1,
    xlab = xlab,
    ylab = ylab,
    main = title
  )

  if (!is.null(caption)) {
    mtext(caption, side = 1, line = 4, cex = 0.75)
  }
}

plot_stacked_exact <- function(exact) {
  if (is.null(exact)) {
    plot.new(); title(main = "Exact Gradient Decomposition"); text(0.5, 0.5, "No data available"); return()
  }

  required <- c("n_state", "objective_ms", "tape_setup_ms", "reverse_pass_ms", "gradient_extract_ms")
  if (!all(required %in% names(exact))) {
    plot.new(); title(main = "Exact Gradient Decomposition"); text(0.5, 0.5, "Missing columns"); return()
  }

  mat <- rbind(
    exact$objective_ms,
    exact$tape_setup_ms,
    exact$reverse_pass_ms,
    exact$gradient_extract_ms
  )

  barplot(
    mat,
    beside = FALSE,
    names.arg = exact$n_state,
    xlab = "Number of latent states",
    ylab = "Time (ms)",
    main = "Exact Laplace Gradient Timing Decomposition"
  )

  legend(
    "topleft",
    legend = c("Laplace objective path", "Tape setup", "Reverse pass", "Gradient extraction"),
    fill = seq_len(4),
    bty = "n",
    cex = 0.8
  )

  mtext(
    "This decomposes exact-gradient cost into the Laplace objective path and AD-gradient extraction phases.",
    side = 1,
    line = 4,
    cex = 0.75
  )
}

summary_table_page <- function(random, state, exact, factor) {
  plot.new()
  text(0.05, 0.94, "Key Benchmark Signals", adj = 0, cex = 1.6, font = 2)

  y <- 0.84

  add_line <- function(label, value) {
    text(0.08, y <<- y, label, adj = 0, cex = 0.95, font = 2)
    text(0.50, y, value, adj = 0, cex = 0.95)
    y <<- y - 0.07
  }

  if (!is.null(random) && nrow(random) > 0) {
    q <- random[random$engine == "quadra", ]
    t <- random[random$engine == "tmb", ]

    if (nrow(q) > 0) {
      add_line("Quadra random-intercept max wall ms", fmt(max(as.numeric(q$total_wall_ms), na.rm = TRUE)))
    }
    if (nrow(t) > 0) {
      add_line("TMB random-intercept max wall ms", fmt(max(as.numeric(t$total_wall_ms), na.rm = TRUE)))
    }
  }

  if (!is.null(state) && nrow(state) > 0) {
    q <- state[state$engine == "quadra", ]
    if (nrow(q) > 0) {
      add_line("State-space max Hessian nnz", fmt(max(as.numeric(q$hessian_nnz), na.rm = TRUE), 0))
      add_line("State-space max fill ratio", fmt(max(as.numeric(q$fill_ratio), na.rm = TRUE)))
    }
  }

  if (!is.null(exact) && nrow(exact) > 0) {
    add_line("Exact-gradient max total ms", fmt(max(as.numeric(exact$total_ms), na.rm = TRUE)))
  }

  if (!is.null(factor) && nrow(factor) > 0) {
    add_line("Best factorization reuse ratio", fmt(min(as.numeric(factor$reuse_ratio), na.rm = TRUE)))
  }

  text(
    0.05, 0.18,
    paste(
      "Caution: timing comparisons are phase-specific.",
      "Quadra workspace/exact-gradient timings and TMB setup/fn/gr timings are not always identical workloads.",
      "The strongest conclusions come from scaling shape, decomposition, and structural diagnostics.",
      sep = "\n"
    ),
    adj = 0,
    cex = 0.85
  )
}

random <- safe_read("benchmarks/normalized/random_intercept_normalized.csv")
state <- safe_read("benchmarks/normalized/state_space_normalized.csv")
exact <- safe_read("benchmarks/exact_laplace_gradient/state_space_exact_gradient_benchmark.csv")
reuse <- safe_read("benchmarks/exact_laplace_gradient/exact_gradient_reuse_benchmark.csv")
factor <- safe_read("benchmarks/exact_laplace_gradient/factorization_reuse_benchmark.csv")
exact_cmp <- safe_read("benchmarks/normalized/exact_gradient_state_space_comparison.csv")

pdf(out, width = 11, height = 8.5)

page_title(
  "Quadra Benchmark Report",
  "Sparse mixed-effects inference, exact-gradient decomposition, and reuse diagnostics"
)

page_text(
  "Executive Summary",
  c(
    "This report emphasizes benchmark interpretation, not just raw timings.",
    "Random-intercept benchmarks measure overhead and simple scaling.",
    "State-space benchmarks stress sparse Hessian structure and latent dimension growth.",
    "Exact-gradient benchmarks decompose Laplace-gradient cost into objective and AD phases.",
    "Reuse benchmarks estimate the value of cached symbolic sparse factorization.",
    "TMB comparisons should be interpreted carefully because exposed timing phases are not identical."
  )
)

summary_table_page(random, state, exact, factor)


page_text(
  "Benchmark Fairness and Limitations",
  c(
    "These benchmarks are best interpreted as architectural and scaling diagnostics rather than definitive framework rankings.",
    "Several comparisons are phase-level comparisons rather than mathematically identical end-to-end workloads.",
    "Quadra benchmarks are native C++ executables, while current TMB benchmarks include R/TMB evaluation pathways.",
    "Exact-gradient implementations are not yet mathematically identical because the derivative of 0.5 * logdet(H_uu) is not currently included in Quadra's exact-gradient path.",
    "Sparse ordering strategies, factorization reuse behavior, compiler flags, and backend libraries may differ between engines.",
    "Quadra currently exposes more internal decomposition metrics than TMB through this benchmark harness.",
    "Some reuse-oriented workloads have not yet been implemented equivalently for both frameworks.",
    "The strongest conclusions currently come from scaling shape, decomposition behavior, sparse structure diagnostics, and reuse experiments."
  )
)

page_text(
  "Benchmark Method Notes",
  c(
    "Quadra timings are generated from C++ benchmark drivers.",
    "TMB timings are generated through R/TMB runners when TMB is available.",
    "CI artifacts include normalized CSVs, raw timing logs, plots, and this PDF report.",
    "Several comparisons are phase-level comparisons rather than full end-to-end equivalence tests.",
    "Structural metrics such as Hessian nonzeros and fill ratio help explain runtime behavior."
  )
)

plot_metric(
  random,
  "n_obs",
  "total_wall_ms",
  "Random-Intercept Total Wall Time",
  "Number of observations",
  "Total wall time (ms)",
  "Useful for overhead/scaling checks; not the strongest sparse-structure benchmark.",
  log_y = TRUE
)

plot_metric(
  random,
  "n_obs",
  "objective_eval_ms",
  "Random-Intercept Objective Evaluation",
  "Number of observations",
  "Objective evaluation time (ms)",
  "Compares objective-evaluation scaling where available.",
  log_y = FALSE
)

plot_metric(
  state,
  "n_state",
  "total_wall_ms",
  "State-Space Total Wall Time",
  "Number of latent states",
  "Total wall time (ms)",
  "State-space models better stress sparse random-effect structure.",
  log_y = TRUE
)

plot_metric(
  state,
  "n_state",
  "hessian_nnz",
  "State-Space Hessian Nonzeros",
  "Number of latent states",
  "nnz(H_uu)",
  "For random-walk structure, nonzeros should grow roughly linearly."
)

plot_metric(
  state,
  "n_state",
  "fill_ratio",
  "State-Space Fill Ratio",
  "Number of latent states",
  "factor nnz / Hessian nnz",
  "Fill ratio helps explain whether factorization cost is driven by symbolic structure."
)

plot_metric(
  state,
  "n_state",
  "factor_nnz",
  "State-Space Factor Nonzeros",
  "Number of latent states",
  "factor nnz",
  "Factor nonzeros are often more predictive of solve/factor cost than raw state dimension."
)

plot_single(
  exact,
  "n_state",
  "total_ms",
  "Exact Laplace Gradient Scaling",
  "Number of latent states",
  "Total exact-gradient time (ms)",
  "This is Quadra's combined Laplace objective plus exact envelope-gradient path."
)

plot_stacked_exact(exact)

plot_single(
  exact,
  "hessian_nnz",
  "reverse_pass_ms",
  "Reverse Pass Cost vs Hessian Structure",
  "nnz(H_uu)",
  "Reverse pass time (ms)",
  "Helps distinguish AD traversal cost from sparse linear algebra cost."
)

plot_single(
  reuse,
  "iteration",
  "total_ms",
  "Exact Gradient Reuse Timing",
  "Iteration",
  "Total exact-gradient time (ms)",
  "Repeated nearby evaluations expose warm/cold behavior."
)

plot_single(
  factor,
  "n_state",
  "reuse_ratio",
  "Sparse Factorization Reuse Ratio",
  "Number of latent states",
  "reuse_ms / fresh_ms",
  "Values below 1 indicate benefit from symbolic factorization reuse."
)

plot_single(
  factor,
  "fill_ratio",
  "reuse_ratio",
  "Reuse Efficiency vs Fill Ratio",
  "Fill ratio",
  "reuse_ms / fresh_ms",
  "Connects structural fill-in to reuse effectiveness."
)

plot_metric(
  exact_cmp,
  "n_state",
  "exact_gradient_ms",
  "Quadra Exact Gradient Comparison Data",
  "Number of latent states",
  "Exact-gradient time (ms)",
  "Shows Quadra exact-gradient rows when normalized comparison data are present."
)

page_text(
  "Supported Conclusions",
  c(
    "Quadra now exposes enough timing decomposition to identify dominant computational phases.",
    "State-space benchmarks are more informative than scalar random-intercept benchmarks for sparse inference.",
    "Factorization reuse is measurable and can be tracked as a first-class performance metric.",
    "Exact-gradient costs can now be studied by objective path, tape setup, reverse pass, and extraction.",
    "TMB comparison infrastructure exists, but headline comparisons should wait until phase-equivalent workloads are finalized."
  )
)

dev.off()

cat("Wrote", out, "\n")
