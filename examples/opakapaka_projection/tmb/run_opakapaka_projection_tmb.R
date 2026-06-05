suppressPackageStartupMessages({
  library(TMB)
})

cat("Synthetic and public-data-safe. Not an official assessment.\n\n")

set.seed(123)

n_years <- 30L
r <- 0.34
K <- 950.0
q_true <- 0.00112
sigma_process <- 0.10
sigma_index <- 0.08
B0 <- 0.82 * K

catch_obs <- rep(35.0, n_years)
B <- numeric(n_years)
B[1] <- B0

for (y in seq_len(n_years - 1L)) {
  pred <- B[y] + r * B[y] * (1.0 - B[y] / K) - catch_obs[y]
  B[y + 1L] <- max(1.0, pred) * exp(rnorm(1L, 0.0, sigma_process * 0.25))
}

index_obs <- q_true * B * exp(rnorm(n_years, 0.0, sigma_index))

cpp <- "examples/opakapaka_projection/tmb/opakapaka_projection_tmb.cpp"
dyn <- sub("\\.cpp$", "", basename(cpp))

compile(cpp, flags = "-O2 -DNDEBUG")
dyn.load(dynlib(file.path("examples/opakapaka_projection/tmb", dyn)))

data <- list(
  index_obs = index_obs,
  catch_obs = catch_obs,
  n_years = n_years
)

parameters <- list(
  log_q = log(0.0010),
  log_B = rep(log(B0), n_years)
)

obj <- MakeADFun(
  data = data,
  parameters = parameters,
  random = "log_B",
  DLL = dyn,
  silent = TRUE
)

lower <- obj$par
upper <- obj$par
lower[] <- -Inf
upper[] <- Inf
lower["log_q"] <- log(1.0e-5)
upper["log_q"] <- log(1.0e-1)

safe_fn <- function(x) {
  val <- tryCatch(obj$fn(x), error = function(e) NA_real_)
  if (!is.finite(val)) return(.Machine$double.xmax / 1e100)
  val
}

safe_gr <- function(x) {
  g <- tryCatch(obj$gr(x), error = function(e) rep(NA_real_, length(x)))
  if (length(g) != length(x) || any(!is.finite(g))) {
    return(rep(0.0, length(x)))
  }
  g
}

t0 <- proc.time()[["elapsed"]]
fit <- nlminb(
  obj$par, safe_fn, safe_gr,
  lower = lower,
  upper = upper,
  control = list(eval.max = 400, iter.max = 400, rel.tol = 1e-10, x.tol = 1e-10)
)
elapsed_ms <- (proc.time()[["elapsed"]] - t0) * 1000.0

obj$fn(fit$par)
obj$gr(fit$par)

cat("TMB synthetic opakapaka-style fit + projection benchmark\n")
cat("=======================================================\n\n")
cat("Fit summary\n")
cat("-----------\n")
cat(sprintf("objective             %.6f\n", fit$objective))
cat(sprintf("convergence           %s\n", fit$convergence))
cat(sprintf("message               %s\n", fit$message))
cat(sprintf("runtime_ms            %.6f\n", elapsed_ms))
cat(sprintf("log_q_hat             %.9f\n", fit$par[["log_q"]]))
cat(sprintf("q_hat                 %.9f\n", exp(fit$par[["log_q"]])))

last_random <- obj$env$last.par.best[-seq_along(fit$par)]
last_B <- exp(tail(last_random, 1L))

proj_years <- 20L
scenarios <- data.frame(
  scenario = c("zero_catch", "status_quo", "high_catch"),
  catch = c(0.0, 35.0, 45.0)
)

projection <- list()
for (i in seq_len(nrow(scenarios))) {
  Bp <- numeric(proj_years + 1L)
  Bp[1] <- last_B
  for (y in seq_len(proj_years)) {
    Bp[y + 1L] <- max(1.0, Bp[y] + r * Bp[y] * (1.0 - Bp[y] / K) - scenarios$catch[i])
  }
  projection[[i]] <- data.frame(
    scenario = scenarios$scenario[i],
    year = 0:proj_years,
    biomass = Bp
  )
}

projection <- do.call(rbind, projection)

outdir <- "examples/opakapaka_projection/outputs"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  data.frame(
    objective = fit$objective,
    convergence = fit$convergence,
    runtime_ms = elapsed_ms,
    log_q_hat = fit$par[["log_q"]],
    q_hat = exp(fit$par[["log_q"]])
  ),
  file.path(outdir, "tmb_synthetic_fit_summary.csv"),
  row.names = FALSE
)

write.csv(
  projection,
  file.path(outdir, "tmb_synthetic_projection_scenarios.csv"),
  row.names = FALSE
)

cat("\nWrote outputs:\n")
cat("  examples/opakapaka_projection/outputs/tmb_synthetic_fit_summary.csv\n")
cat("  examples/opakapaka_projection/outputs/tmb_synthetic_projection_scenarios.csv\n")
