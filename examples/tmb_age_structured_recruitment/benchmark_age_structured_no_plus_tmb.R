#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[[1]]) else 10L
lengths <- if (length(args) >= 2) as.integer(strsplit(args[[2]], ",")[[1]]) else c(25L, 50L, 100L, 250L, 500L, 1000L)
n_ages <- if (length(args) >= 3) as.integer(args[[3]]) else 10L

cat("TMB no-plus age-structured recruitment benchmark\n")
cat("================================================\n")
cat("reps per n =", reps, ", ages =", n_ages, "\n\n")

if (!requireNamespace("TMB", quietly = TRUE)) {
  cat("TMB is not installed. Skipping.\n")
  quit(status = 0)
}

library(TMB)

template <- file.path("examples", "tmb_age_structured_recruitment", "age_structured_recruitment_no_plus_tmb.cpp")
dynlib_name <- "age_structured_recruitment_no_plus_tmb"

TMB::compile(template, flags = "-O2")
dyn.load(TMB::dynlib(file.path("examples", "tmb_age_structured_recruitment", dynlib_name)))

logistic <- function(x) 1 / (1 + exp(-x))

make_index <- function(n_years, n_ages) {
  R0 <- 1000
  M <- 0.20
  q <- 0.001

  N <- R0 * exp(-M * 0:(n_ages - 1))
  sel <- logistic(((1:n_ages) - 4) / 0.8)

  index <- numeric(n_years)

  for (y0 in seq_len(n_years)) {
    y <- y0 - 1
    vulnerable <- sum(sel * N)

    obs_error <- 0.04 * sin(2 * pi * y / 13) +
      0.02 * cos(2 * pi * y / 7)

    index[y0] <- q * vulnerable * exp(obs_error)

    recruitment_dev <- 0.15 * sin(2 * pi * y / 17)

    N_next <- numeric(n_ages)
    N_next[1] <- R0 * exp(recruitment_dev)
    for (a in 2:n_ages) {
      N_next[a] <- N[a - 1] * exp(-M)
    }

    # No plus group.

    N <- N_next
  }

  index
}

cat(sprintf("%8s%14s%14s\n", "n", "objective", "avg_ms"))

for (n in lengths) {
  data <- list(index_obs = make_index(n, n_ages), n_ages = n_ages)

  parameters <- list(
    log_R0 = log(1000),
    log_M = log(0.20),
    log_q = log(0.001),
    log_sigma_R = log(0.35),
    log_sigma_index = log(0.10),
    x = rep(0, n)
  )

  obj <- TMB::MakeADFun(
    data = data,
    parameters = parameters,
    random = "x",
    DLL = dynlib_name,
    silent = TRUE
  )

  last <- obj$fn()
  gc()

  t0 <- proc.time()
  for (i in seq_len(reps)) {
    last <- obj$fn()
  }
  t1 <- proc.time()

  avg_ms <- as.numeric((t1 - t0)[["elapsed"]]) * 1000 / reps

  cat(sprintf("%8d%14.6f%14.6f\n", n, last, avg_ms))
}
