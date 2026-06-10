#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[[1]]) else 10L
lengths <- if (length(args) >= 2) as.integer(strsplit(args[[2]], ",")[[1]]) else c(25L, 50L, 100L, 250L)

cat("TMB scaled fixed-theta state-space surplus benchmark\n")
cat("===================================================\n")
cat("reps per n =", reps, "\n\n")

if (!requireNamespace("TMB", quietly = TRUE)) {
  cat("TMB is not installed. Skipping.\n")
  quit(status = 0)
}

library(TMB)

template <- file.path("examples", "state_space_surplus_production", "tmb", "state_space_surplus_tmb.cpp")
dynlib_name <- "state_space_surplus_tmb"

if (!file.exists(TMB::dynlib(file.path("examples", "state_space_surplus_production", "tmb", dynlib_name)))) {
  TMB::compile(template, flags = "-O2")
}

dyn.load(TMB::dynlib(file.path("examples", "state_space_surplus_production", "tmb", dynlib_name)))

make_scaled_data <- function(n) {
  r <- 0.5
  K <- 700
  q <- 0.0024
  B <- 0.90 * K
  catch_observed <- numeric(n)
  index_observed <- numeric(n)

  for (t0 in seq_len(n)) {
    t <- t0 - 1
    seasonal <- sin(2 * pi * t / 17)
    trend <- 1 + 0.10 * sin(2 * pi * t / 53)
    C <- 88 * trend + 18 * seasonal
    catch_observed[t0] <- max(40, C)

    obs_error <- 0.05 * sin(2 * pi * t / 11) +
      0.025 * cos(2 * pi * t / 7)

    index_observed[t0] <- q * B * exp(obs_error)

    if (t0 < n) {
      production <- r * B * (1 - B / K)
      B <- max(B + production - catch_observed[t0], 1e-9)
    }
  }

  list(catch_observed = catch_observed, index_observed = index_observed)
}

cat(sprintf("%8s%14s%14s\n", "n", "objective", "avg_ms"))

for (n in lengths) {
  data <- make_scaled_data(n)

  parameters <- list(
    log_r = log(0.5),
    log_K = log(700.0),
    log_q = log(0.0024),
    log_sigma_process = log(0.15),
    log_sigma_index = log(0.10),
    logit_B0_frac = log(0.90 / 0.10),
    u = rep(0, n - 1)
  )

  obj <- TMB::MakeADFun(
    data = data,
    parameters = parameters,
    random = "u",
    DLL = dynlib_name,
    silent = TRUE
  )

  last <- obj$fn()
  gc()

  t0 <- proc.time()
  for (i in seq_len(reps)) last <- obj$fn()
  t1 <- proc.time()

  avg_ms <- as.numeric((t1 - t0)[["elapsed"]]) * 1000 / reps

  cat(sprintf("%8d%14.6f%14.6f\n", n, last, avg_ms))
}
