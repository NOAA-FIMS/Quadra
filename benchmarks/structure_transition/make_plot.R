#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
d <- read.csv(args[[1]])

png(args[[2]], width = 1800, height = 800, res = 160)
par(mfrow = c(1, 2))

cases <- unique(d$case)
symbols <- seq_along(cases) + 14

plot(NULL, xlim = range(d$n), ylim = range(d$ms), log = "xy",
     xlab = "Latent dimension", ylab = "Milliseconds per update",
     main = "Automatic backend transition")
for (i in seq_along(cases)) {
  x <- d[d$case == cases[[i]], ]
  lines(x$n, x$ms, type = "b", pch = symbols[[i]], lwd = 2)
}
legend("topleft", legend = cases, pch = symbols, lwd = 2, bty = "n", cex = 0.8)
grid()

plot(NULL, xlim = range(d$n), ylim = range(d$peak_rss_mib), log = "xy",
     xlab = "Latent dimension", ylab = "Peak RSS (MiB)",
     main = "Process peak resident memory")
for (i in seq_along(cases)) {
  x <- d[d$case == cases[[i]], ]
  lines(x$n, x$peak_rss_mib, type = "b", pch = symbols[[i]], lwd = 2)
}
legend("topleft", legend = cases, pch = symbols, lwd = 2, bty = "n", cex = 0.8)
grid()

dev.off()
cat("wrote", args[[2]], "\n")
