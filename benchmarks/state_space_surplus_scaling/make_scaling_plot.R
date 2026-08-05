#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("usage: make_scaling_plot.R results.csv output.png")
}

csv_path <- args[[1]]
out_path <- args[[2]]

d <- read.csv(csv_path)

png(out_path, width = 1800, height = 800, res = 160)

ylim <- range(c(d$quadra_ms, d$tmb_ms), finite = TRUE)
par(mfrow = c(1, 2))

plot(
  d$n,
  d$quadra_ms,
  type = "b",
  log = "y",
  pch = 16,
  lwd = 2,
  ylim = ylim,
  xlab = "Number of years / latent-state scale",
  ylab = "Milliseconds per fixed-theta Laplace evaluation",
  main = "Runtime"
)

lines(d$n, d$tmb_ms, type = "b", pch = 17, lwd = 2)

legend(
  "topleft",
  legend = c("Quadra persistent tridiagonal", "TMB AD/Laplace"),
  pch = c(16, 17),
  lwd = 2,
  bty = "n"
)

grid()

rss_ylim <- range(c(d$quadra_peak_rss_mib, d$tmb_peak_rss_mib), finite = TRUE)
plot(
  d$n,
  d$quadra_peak_rss_mib,
  type = "b",
  log = "y",
  pch = 16,
  lwd = 2,
  ylim = rss_ylim,
  xlab = "Number of years / latent-state scale",
  ylab = "Peak RSS (MiB)",
  main = "Peak resident memory"
)
lines(d$n, d$tmb_peak_rss_mib, type = "b", pch = 17, lwd = 2)
legend(
  "topleft",
  legend = c("Quadra persistent tridiagonal", "TMB AD/Laplace"),
  pch = c(16, 17),
  lwd = 2,
  bty = "n"
)
grid()

dev.off()

cat("wrote", out_path, "\n")
