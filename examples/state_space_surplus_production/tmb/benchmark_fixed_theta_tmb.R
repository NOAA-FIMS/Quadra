#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[[1]]) else 20L

cat("TMB fixed-theta state-space surplus benchmark\n")
cat("============================================\n")

if (!requireNamespace("TMB", quietly = TRUE)) {
  cat("TMB is not installed. Skipping.\n")
  quit(status = 0)
}

library(TMB)

template <- file.path("examples", "state_space_surplus_production", "tmb", "state_space_surplus_tmb.cpp")
dynlib_name <- "state_space_surplus_tmb"

catch_observed <- c(
  80, 88, 95, 105, 115, 125, 130, 128,
  120, 110, 100, 90, 82, 78, 75
)

index_observed <- c(
  1.55, 1.50, 1.43, 1.34, 1.22, 1.10, 0.98, 0.88,
  0.81, 0.78, 0.77, 0.80, 0.84, 0.88, 0.92
)

parameters <- list(
  log_r = log(0.5),
  log_K = log(700.0),
  log_q = log(0.0024),
  log_sigma_process = log(0.15),
  log_sigma_index = log(0.10),
  logit_B0_frac = log(0.90 / 0.10),
  u = rep(0, length(catch_observed) - 1)
)

data <- list(catch_observed = catch_observed, index_observed = index_observed)

if (!file.exists(TMB::dynlib(file.path("examples", "state_space_surplus_production", "tmb", dynlib_name)))) {
  cat("Compiling TMB template...\n")
  TMB::compile(template, flags = "-O2")
}

dyn.load(TMB::dynlib(file.path("examples", "state_space_surplus_production", "tmb", dynlib_name)))

obj <- TMB::MakeADFun(
  data = data,
  parameters = parameters,
  random = "u",
  DLL = dynlib_name,
  silent = TRUE
)

# Warmup.
last <- obj$fn()

gc()
t0 <- proc.time()
for (i in seq_len(reps)) {
  last <- obj$fn()
}
t1 <- proc.time()

elapsed <- as.numeric((t1 - t0)[["elapsed"]])
avg_ms <- elapsed * 1000 / reps

cat("reps =", reps, "\n")
cat("objective =", format(last, digits = 12), "\n")
cat("total_ms =", format(elapsed * 1000, digits = 12), "\n")
cat("avg_ms =", format(avg_ms, digits = 12), "\n")
cat("expected Quadra dense FD objective ~= -10.642184\n")
