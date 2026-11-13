# Optional scientific PDF; uncertainty itself requires no R dependency.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Usage: Rscript plot_reference_uncertainty.R INPUT_CSV OUTPUT_PDF")
values <- read.csv(args[1], stringsAsFactors=FALSE)
required <- c("benchmark", "quantity", "estimate", "standard_error", "lower95", "upper95")
if (!all(required %in% names(values))) stop("Invalid reference-point uncertainty table")
dir.create(dirname(args[2]), recursive=TRUE, showWarnings=FALSE)
pdf(args[2], width=8.5, height=11, title="Quadra swordfish reference-point uncertainty")
par(mfrow=c(3,1), mar=c(4,6,2.6,1.4), oma=c(3.2,0,3.4,0), las=1, family="sans")

plot_points <- function(quantity, title, label, log_scale=FALSE) {
  # Each panel uses native units and the model's positive pointwise 95% interval.
  rows <- values[values$quantity == quantity & values$estimate > 0, ]
  rows <- rows[match(intersect(c("unfished", "Btgt20", "SPR20", "MSY"), rows$benchmark), rows$benchmark), ]
  if (!nrow(rows) || any(!is.finite(as.matrix(rows[,c("estimate","lower95","upper95")])))) stop("Invalid interval values")
  positions <- rev(seq_len(nrow(rows)))
  limits <- range(c(rows$lower95, rows$upper95))
  if (log_scale) limits <- limits*c(.92,1.08) else limits <- limits+c(-1,1)*diff(limits)*.08
  plot(NA, xlim=limits, ylim=c(.5,nrow(rows)+.5), log=if(log_scale) "x" else "",
       yaxt="n", ylab="", xlab=label, main=title, bty="n", cex.main=1.1)
  abline(h=positions, col="#E6EAED", lty=1)
  segments(rows$lower95, positions, rows$upper95, positions, col="#147C89", lwd=2.5)
  segments(rows$lower95, positions-.07, rows$lower95, positions+.07, col="#147C89", lwd=1.5)
  segments(rows$upper95, positions-.07, rows$upper95, positions+.07, col="#147C89", lwd=1.5)
  points(rows$estimate, positions, pch=19, col="#173D5E", cex=1.1)
  axis(2, at=positions, labels=rows$benchmark, tick=FALSE)
}

plot_points("SSB_tonnes", "Equilibrium female spawning biomass", "Tonnes (logarithmic scale)", TRUE)
plot_points("dead_catch_tonnes", "Equilibrium dead catch", "Tonnes per year")
plot_points("annual_F_SS3_units4_ages1_10", "Equilibrium annual fishing mortality", "Annual F: SS3 units 4, pooled ages 1-10")
mtext("North Pacific swordfish: reference-point uncertainty", outer=TRUE, side=3, line=1.4, font=2, cex=1.2)
mtext("Points: estimates. Bars: pointwise 95% log-scale Wald intervals.", outer=TRUE, side=1, line=.5, cex=.85)
mtext("Full correlated local delta method; conditional on model/fixed settings. No process uncertainty.", outer=TRUE, side=1, line=1.7, cex=.75)
dev.off()
