#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(TMB))

data_candidates <- c(
  "examples/NMFS/sefsc_red_snapper/data/red_snapper_synthetic_observations.csv",
  "examples/NMFS/sefsc_red_snapper/data/synthetic_red_snapper_observations.csv",
  "examples/NMFS/sefsc_red_snapper/data/red_snapper_observations.csv"
)

data_path <- data_candidates[file.exists(data_candidates)][1]
if (is.na(data_path)) {
  stop("Could not find SEFSC red snapper observation CSV")
}

obs <- read.csv(data_path)

catch_col <- grep("catch", names(obs), value = TRUE)[1]
index_col <- grep("index", names(obs), value = TRUE)[1]
age_cols <- grep("^age_|^age[0-9]+|comp", names(obs), value = TRUE)
age_cols <- age_cols[sapply(obs[age_cols], is.numeric)]

if (is.na(catch_col) || is.na(index_col) || length(age_cols) == 0) {
  stop("Could not infer catch/index/age-composition columns from data CSV")
}

catch_obs <- as.numeric(obs[[catch_col]])
index_obs <- as.numeric(obs[[index_col]])
age_comp_obs <- as.matrix(obs[, age_cols, drop = FALSE])
age_comp_obs <- age_comp_obs / rowSums(age_comp_obs)

cpp <- "examples/NMFS/sefsc_red_snapper/tmb/red_snapper_tmb.cpp"
TMB::compile(cpp)
dyn.load(TMB::dynlib(sub("\\.cpp$", "", cpp)))

parameters <- list(
  log_r0 = log(1200.0),
  log_fbar = log(0.025),
  log_q = log(0.00005),
  logit_sel_a50 = 0.0,
  log_sel_slope = log(1.2),
  log_rec_dev = rep(0.0, length(catch_obs))
)

obj <- MakeADFun(
  data = list(catch_obs = catch_obs, index_obs = index_obs, age_comp_obs = age_comp_obs),
  parameters = parameters,
  random = "log_rec_dev",
  DLL = "red_snapper_tmb",
  silent = TRUE
)

fit <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1000, iter.max = 1000))
pl <- obj$env$parList()

summary_path <- "examples/NMFS/sefsc_red_snapper/outputs/tmb_fit_summary.csv"
out <- data.frame(
  field = c("objective", "convergence", "message", "log_r0", "r0",
            "log_fbar", "fbar", "log_q", "q", "logit_sel_a50",
            "sel_a50", "log_sel_slope", "sel_slope", "random_effects"),
  value = c(fit$objective, fit$convergence, fit$message,
            pl$log_r0, exp(pl$log_r0),
            pl$log_fbar, exp(pl$log_fbar),
            pl$log_q, exp(pl$log_q),
            pl$logit_sel_a50,
            1.0 + 9.0 / (1.0 + exp(-pl$logit_sel_a50)),
            pl$log_sel_slope, exp(pl$log_sel_slope),
            length(pl$log_rec_dev))
)
write.csv(out, summary_path, row.names = FALSE, quote = FALSE)

rec_path <- "examples/NMFS/sefsc_red_snapper/outputs/tmb_recruitment_deviations.csv"
write.csv(data.frame(year = seq_along(pl$log_rec_dev),
                     log_rec_dev = as.numeric(pl$log_rec_dev),
                     rec_multiplier = exp(as.numeric(pl$log_rec_dev))),
          rec_path, row.names = FALSE, quote = FALSE)

cat("wrote:", summary_path, "\n")
cat("wrote:", rec_path, "\n")
