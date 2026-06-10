#!/usr/bin/env Rscript

cat("TMB state-space surplus production comparison\n")
cat("=============================================\n\n")

if (!requireNamespace("TMB", quietly = TRUE)) {
  cat("TMB is not installed. Skipping TMB comparison.\n")
  cat("Install in R with install.packages('TMB') if needed.\n")
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

data <- list(
  catch_observed = catch_observed,
  index_observed = index_observed
)

cat("Compiling TMB template...\n")
TMB::compile(template, flags = "-O2")
dyn.load(TMB::dynlib(file.path("examples", "state_space_surplus_production", "tmb", dynlib_name)))

cat("Building objective with u as random effects...\n")
obj <- TMB::MakeADFun(
  data = data,
  parameters = parameters,
  random = "u",
  DLL = dynlib_name,
  silent = TRUE
)

cat("Evaluating initial marginal objective...\n")
t0 <- proc.time()
initial_obj <- obj$fn()
initial_grad <- obj$gr()
t1 <- proc.time()

cat("Optimizing fixed effects with nlminb...\n")
opt <- nlminb(
  start = obj$par,
  objective = obj$fn,
  gradient = obj$gr,
  control = list(iter.max = 200, eval.max = 500)
)

t2 <- proc.time()

rep <- TMB::sdreport(obj)
u_hat <- obj$env$last.par.best[names(obj$env$last.par.best) == "u"]

cat("\nTMB results\n")
cat("-----------\n")
cat("initial marginal objective:", format(initial_obj, digits = 12), "\n")
cat("initial gradient norm:", format(sqrt(sum(initial_grad^2)), digits = 12), "\n")
cat("optimized marginal objective:", format(opt$objective, digits = 12), "\n")
cat("nlminb convergence:", opt$convergence, "\n")
cat("nlminb message:", opt$message, "\n")
cat("initial eval elapsed seconds:", as.numeric((t1 - t0)[["elapsed"]]), "\n")
cat("optimization elapsed seconds:", as.numeric((t2 - t1)[["elapsed"]]), "\n")

cat("\nFixed effects estimates\n")
print(opt$par)

cat("\nRandom effects summary\n")
cat("n:", length(u_hat), "\n")
cat("mean:", mean(u_hat), "\n")
cat("sd:", sd(u_hat), "\n")
cat("min:", min(u_hat), "\n")
cat("max:", max(u_hat), "\n")
cat("norm:", sqrt(sum(u_hat^2)), "\n")

cat("\nFirst random effects\n")
print(round(u_hat, 6))

cat("\nADREPORT summary\n")
print(summary(rep, "report"))

cat("\nQuadra dense Laplace fixed-theta reference from C++ run should be approximately:\n")
cat("joint(theta, u_hat) = -33.863092\n")
cat("laplace marginal nll = -10.642184\n")
cat("Note: TMB initial marginal objective should be close to the fixed-theta Laplace value.\n")
