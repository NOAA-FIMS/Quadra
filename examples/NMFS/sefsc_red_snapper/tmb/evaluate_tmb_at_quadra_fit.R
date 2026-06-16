library(TMB)

obs <- read.csv("examples/NMFS/sefsc_red_snapper/data/synthetic_red_snapper_observations.csv")
catch_obs <- obs$catch_mt
index_obs <- obs$index
age_cols <- grep("^age[0-9]+$", names(obs), value = TRUE)
age_comp_obs <- as.matrix(obs[, age_cols, drop = FALSE])
age_comp_obs <- age_comp_obs / rowSums(age_comp_obs)

cpp <- "examples/NMFS/sefsc_red_snapper/tmb/red_snapper_tmb.cpp"
dyn <- sub("\\.cpp$", "", basename(cpp))
if (!file.exists(file.path("examples/NMFS/sefsc_red_snapper/tmb", paste0(dyn, .Platform$dynlib.ext)))) {
  TMB::compile(cpp)
}
dyn.load(dynlib(file.path("examples/NMFS/sefsc_red_snapper/tmb", dyn)))

qsum <- read.csv("examples/NMFS/sefsc_red_snapper/outputs/quadra_fit_summary.csv")
qval <- setNames(qsum$value, qsum$field)

qrec <- read.csv("examples/NMFS/sefsc_red_snapper/outputs/recruitment_deviations.csv")

parameters <- list(
  log_r0 = as.numeric(qval["log_r0"]),
  log_fbar = as.numeric(qval["log_fbar"]),
  log_q = as.numeric(qval["log_q"]),
  logit_sel_a50 = as.numeric(qval["logit_sel_a50"]),
  log_sel_slope = as.numeric(qval["log_sel_slope"]),
  log_rec_dev = qrec$log_rec_dev
)

obj <- MakeADFun(
  data = list(
    catch_obs = catch_obs,
    index_obs = index_obs,
    age_comp_obs = age_comp_obs
  ),
  parameters = parameters,
  random = "log_rec_dev",
  DLL = "red_snapper_tmb",
  silent = TRUE
)

cat("TMB Laplace objective at Quadra fit:", obj$fn(), "\n")
cat("TMB gradient at Quadra fit:\n")
print(obj$gr())

H <- obj$env$spHess(random = TRUE)
H <- as.matrix(H)
write.csv(H,
  "examples/NMFS/sefsc_red_snapper/outputs/tmb_Huu_at_quadra_fit.csv",
  row.names = FALSE
)
cat(
  "TMB Huu logdet:",
  as.numeric(determinant(H, logarithm = TRUE)$modulus),
  "\n"
)

rep <- obj$report()
write.csv(
  data.frame(
    component = c("fixed_prior_nll", "rec_prior_nll", "index_nll", "catch_nll", "age_comp_nll"),
    value = c(rep$fixed_prior_nll, rep$rec_prior_nll, rep$index_nll, rep$catch_nll, rep$age_comp_nll)
  ),
  "examples/NMFS/sefsc_red_snapper/outputs/tmb_components_at_quadra_fit.csv",
  row.names = FALSE
)
print(read.csv("examples/NMFS/sefsc_red_snapper/outputs/tmb_components_at_quadra_fit.csv"))

cat("TMB random gradient at Quadra u:\n")
cat("length last.par:", length(obj$env$last.par), "\n")
cat("length last.par.best:", length(obj$env$last.par.best), "\n")
cat("length random:", length(obj$env$random), "\n")
cat("random indices:\n")
print(obj$env$random)

u_from_last <- obj$env$last.par[obj$env$random]
u_from_best <- obj$env$last.par.best[obj$env$random]

cat("max |last random - Quadra u|:\n")
print(max(abs(u_from_last - qrec$log_rec_dev)))

cat("max |best random - Quadra u|:\n")
print(max(abs(u_from_best - qrec$log_rec_dev)))

g <- obj$gr()
cat("TMB grad norm at Quadra fit:", sqrt(sum(g * g)), "\n")

obj_joint <- MakeADFun(
  data = list(
    catch_obs = catch_obs,
    index_obs = index_obs,
    age_comp_obs = age_comp_obs
  ),
  parameters = parameters,
  DLL = "red_snapper_tmb",
  silent = TRUE
)

joint_gr <- obj_joint$gr()[1:5]
laplace_gr <- obj$gr()
implied_logdet_gr <- laplace_gr - joint_gr

cat("TMB joint fixed gradient at Quadra theta/u:\n")
print(joint_gr)

cat("TMB Laplace gradient at Quadra fit:\n")
print(laplace_gr)

cat("TMB implied logdet contribution:\n")
print(implied_logdet_gr)


cat("\nTMB profiled logdet FD at Quadra fit:\n")

fixed_names <- c("log_r0", "log_fbar", "log_q", "logit_sel_a50", "log_sel_slope")
theta0 <- as.numeric(qval[fixed_names])
names(theta0) <- fixed_names
u0 <- qrec$log_rec_dev

make_obj_for_theta <- function(theta_vec, u_start = u0) {
  pars <- list(
    log_r0 = as.numeric(theta_vec["log_r0"]),
    log_fbar = as.numeric(theta_vec["log_fbar"]),
    log_q = as.numeric(theta_vec["log_q"]),
    logit_sel_a50 = as.numeric(theta_vec["logit_sel_a50"]),
    log_sel_slope = as.numeric(theta_vec["log_sel_slope"]),
    log_rec_dev = u_start
  )

  MakeADFun(
    data = list(
      catch_obs = catch_obs,
      index_obs = index_obs,
      age_comp_obs = age_comp_obs
    ),
    parameters = pars,
    random = "log_rec_dev",
    DLL = "red_snapper_tmb",
    silent = TRUE
  )
}

get_profiled_u_and_logdet <- function(theta_vec, u_start = u0) {
  o <- make_obj_for_theta(theta_vec, u_start)

  # Force evaluation so TMB performs the inner random-effect optimization.
  invisible(o$fn())

  # In this TMB version, profiled random modes are stored in last.par[random].
  u_prof <- o$env$last.par[o$env$random]

  H <- as.matrix(o$env$spHess(random = TRUE))
  logdet <- as.numeric(determinant(H, logarithm = TRUE)$modulus)

  list(u = u_prof, logdet = logdet, obj = o)
}

eps <- 1e-5
profiled_logdet_fd <- numeric(length(theta0))
profiled_u_fd_norm <- numeric(length(theta0))
names(profiled_logdet_fd) <- fixed_names
names(profiled_u_fd_norm) <- fixed_names

for (j in seq_along(theta0)) {
  th_plus <- theta0
  th_minus <- theta0
  th_plus[j] <- th_plus[j] + eps
  th_minus[j] <- th_minus[j] - eps

  plus <- get_profiled_u_and_logdet(th_plus, u0)
  minus <- get_profiled_u_and_logdet(th_minus, u0)

  profiled_logdet_fd[j] <- 0.5 * (plus$logdet - minus$logdet) / (2 * eps)

  # This is du*/dtheta_j from true profiling, useful for comparison later.
  u_fd <- (plus$u - minus$u) / (2 * eps)
  profiled_u_fd_norm[j] <- sqrt(sum(u_fd * u_fd))
}

cat("0.5 * profiled logdet FD gradient:\n")
print(profiled_logdet_fd)

cat("profiled u FD column norms:\n")
print(profiled_u_fd_norm)

cat("TMB implied logdet contribution from obj$gr - joint_gr:\n")
print(implied_logdet_gr)

cat("difference: profiled FD - implied TMB logdet contribution:\n")
print(profiled_logdet_fd - implied_logdet_gr)


cat("\nTMB manual-random-optimized profiled logdet FD:\n")

fixed_names <- c("log_r0", "log_fbar", "log_q", "logit_sel_a50", "log_sel_slope")
theta0 <- as.numeric(qval[fixed_names])
names(theta0) <- fixed_names
u0 <- qrec$log_rec_dev

make_joint_obj <- function(theta_vec, u_vec) {
  pars <- list(
    log_r0 = as.numeric(theta_vec["log_r0"]),
    log_fbar = as.numeric(theta_vec["log_fbar"]),
    log_q = as.numeric(theta_vec["log_q"]),
    logit_sel_a50 = as.numeric(theta_vec["logit_sel_a50"]),
    log_sel_slope = as.numeric(theta_vec["log_sel_slope"]),
    log_rec_dev = u_vec
  )

  MakeADFun(
    data = list(
      catch_obs = catch_obs,
      index_obs = index_obs,
      age_comp_obs = age_comp_obs
    ),
    parameters = pars,
    DLL = "red_snapper_tmb",
    silent = TRUE
  )
}

profile_u_for_theta <- function(theta_vec, u_start) {
  joint <- make_joint_obj(theta_vec, u_start)

  full_start <- c(theta_vec, u_start)
  ntheta <- length(theta_vec)
  nu <- length(u_start)

  fn_u <- function(u_vec) {
    par <- c(theta_vec, u_vec)
    joint$fn(par)
  }

  gr_u <- function(u_vec) {
    par <- c(theta_vec, u_vec)
    as.numeric(joint$gr(par)[(ntheta + 1):(ntheta + nu)])
  }

  opt <- nlminb(
    start = u_start,
    objective = fn_u,
    gradient = gr_u,
    control = list(
      eval.max = 1000,
      iter.max = 1000,
      rel.tol = 1e-12,
      x.tol = 1e-12
    )
  )

  list(
    u = opt$par,
    objective = opt$objective,
    convergence = opt$convergence,
    message = opt$message,
    grad_norm = sqrt(sum(gr_u(opt$par)^2)),
    joint = joint
  )
}

logdet_at_theta_u <- function(theta_vec, u_vec) {
  # Use random-enabled TMB object only as a convenient Huu provider at fixed theta/u.
  o <- make_obj_for_theta(theta_vec, u_vec)
  invisible(o$fn())
  H <- as.matrix(o$env$spHess(random = TRUE))
  as.numeric(determinant(H, logarithm = TRUE)$modulus)
}

eps <- 1e-5
manual_profiled_logdet_fd <- numeric(length(theta0))
manual_u_fd_norm <- numeric(length(theta0))
manual_u_opt_grad_norm_plus <- numeric(length(theta0))
manual_u_opt_grad_norm_minus <- numeric(length(theta0))
manual_u_opt_conv_plus <- integer(length(theta0))
manual_u_opt_conv_minus <- integer(length(theta0))

names(manual_profiled_logdet_fd) <- fixed_names
names(manual_u_fd_norm) <- fixed_names
names(manual_u_opt_grad_norm_plus) <- fixed_names
names(manual_u_opt_grad_norm_minus) <- fixed_names
names(manual_u_opt_conv_plus) <- fixed_names
names(manual_u_opt_conv_minus) <- fixed_names

for (j in seq_along(theta0)) {
  th_plus <- theta0
  th_minus <- theta0
  th_plus[j] <- th_plus[j] + eps
  th_minus[j] <- th_minus[j] - eps

  plus <- profile_u_for_theta(th_plus, u0)
  minus <- profile_u_for_theta(th_minus, u0)

  ld_plus <- logdet_at_theta_u(th_plus, plus$u)
  ld_minus <- logdet_at_theta_u(th_minus, minus$u)

  manual_profiled_logdet_fd[j] <- 0.5 * (ld_plus - ld_minus) / (2 * eps)
  manual_u_fd <- (plus$u - minus$u) / (2 * eps)

  manual_u_fd_norm[j] <- sqrt(sum(manual_u_fd * manual_u_fd))
  manual_u_opt_grad_norm_plus[j] <- plus$grad_norm
  manual_u_opt_grad_norm_minus[j] <- minus$grad_norm
  manual_u_opt_conv_plus[j] <- plus$convergence
  manual_u_opt_conv_minus[j] <- minus$convergence
}

cat("0.5 * manually profiled logdet FD gradient:\n")
print(manual_profiled_logdet_fd)

cat("manual profiled u FD column norms:\n")
print(manual_u_fd_norm)

cat("random optimizer convergence plus/minus:\n")
print(manual_u_opt_conv_plus)
print(manual_u_opt_conv_minus)

cat("random optimizer gradient norms plus:\n")
print(manual_u_opt_grad_norm_plus)

cat("random optimizer gradient norms minus:\n")
print(manual_u_opt_grad_norm_minus)

cat("TMB implied logdet contribution from obj$gr - joint_gr:\n")
print(implied_logdet_gr)

cat("difference: manual profiled FD - implied TMB logdet contribution:\n")
print(manual_profiled_logdet_fd - implied_logdet_gr)
