# suppressPackageStartupMessages({
library(TMB)
# })

cat("Synthetic and public-data-safe. Not an official assessment.\n\n")

# Shared synthetic/public-data-safe dataset used by the Quadra example.
# This keeps the TMB and Quadra objective comparisons apples-to-apples.
data_csv <- read.csv("examples/NMFS/pifsc_opakapaka/data/synthetic_opakapaka_projection_data.csv")

data_csv$index <- as.numeric(data_csv$index)
data_csv$catch_mt <- as.numeric(data_csv$catch_mt)

fit_rows <- is.finite(data_csv$index) & is.finite(data_csv$catch_mt)
if ("phase" %in% names(data_csv)) {
  fit_rows <- fit_rows & data_csv$phase == "history"
}
data_fit <- data_csv[fit_rows, , drop = FALSE]

index_obs <- data_fit$index
catch_obs <- data_fit$catch_mt
n_years <- length(index_obs)

cat("Loaded shared CSV fit rows:", n_years, "

")

sigma_process <- 0.10
sigma_index <- 0.08
sigma_initial <- 0.15


cpp <- "examples/NMFS/pifsc_opakapaka/tmb/opakapaka_projection_tmb.cpp"
dyn <- sub("\\.cpp$", "", basename(cpp))

compile(cpp, flags = "-O2 -DNDEBUG -std=gnu++17")
dyn.load(dynlib(file.path("examples/NMFS/pifsc_opakapaka/tmb", dyn)))

data <- list(
  index_obs = index_obs,
  catch_obs = catch_obs,
  n_years = n_years
)
log_q_init <- log(0.001)
log_B_init <- log(779.0 - 1.425 * (seq_len(n_years) - 1))

parameters <- list(
  log_q = log_q_init,
  log_r = log(0.34),
  log_K = log(950.0),
  log_B = log_B_init


)

obj <- MakeADFun(
  data = data,
  parameters = parameters,
  random = "log_B",
  DLL = "opakapaka_projection_tmb",
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
  if (!is.finite(val)) {
    return(.Machine$double.xmax / 1e100)
  }
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

# Objective convention diagnostics.
# obj$env$last.par.best contains fixed + random estimates on the joint scale.
joint_at_mode <- tryCatch(
  {
    full_par <- obj$env$last.par.best
    if (is.null(full_par)) full_par <- obj$env$last.par
    as.numeric(obj$env$f(full_par))
  },
  error = function(e) NA_real_
)

laplace_adjustment_reported <- fit$objective - joint_at_mode

cat(sprintf("joint_at_mode_no_constants %.6f\n", joint_at_mode))
cat(sprintf("reported_minus_joint       %.6f\n", laplace_adjustment_reported))

cat(sprintf("convergence           %s\n", fit$convergence))
cat(sprintf("message               %s\n", fit$message))
cat(sprintf("runtime_ms            %.6f\n", elapsed_ms))
cat(sprintf("log_q_hat             %.9f\n", fit$par[["log_q"]]))
cat(sprintf("q_hat                 %.9f\n", exp(fit$par[["log_q"]])))
cat(sprintf("r_hat                 %.9f\n", exp(fit$par[["log_r"]])))
cat(sprintf("K_hat                 %.9f\n", exp(fit$par[["log_K"]])))

r <- exp(fit$par[["log_r"]])
K <- exp(fit$par[["log_K"]])

last_random <- obj$env$last.par.best[obj$env$random]
last_B <- exp(tail(last_random, 1L))

outdir <- "examples/NMFS/pifsc_opakapaka/outputs"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  data.frame(
    index = seq_along(last_random) - 1,
    log_B = as.numeric(last_random),
    B = exp(as.numeric(last_random))
  ),
  file.path(outdir, "tmb_fitted_states.csv"),
  row.names = FALSE
)

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


write.csv(
  data.frame(
    objective = fit$objective,
    convergence = fit$convergence,
    runtime_ms = elapsed_ms,
    log_q_hat = fit$par[["log_q"]],
    q_hat = exp(fit$par[["log_q"]]),
    log_r_hat = fit$par[["log_r"]],
    r_hat = r,
    log_K_hat = fit$par[["log_K"]],
    K_hat = K
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
cat("  examples/NMFS/pifsc_opakapaka/outputs/tmb_synthetic_fit_summary.csv\n")
cat("  examples/NMFS/pifsc_opakapaka/outputs/tmb_synthetic_projection_scenarios.csv\n")
