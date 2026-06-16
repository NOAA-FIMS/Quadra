#!/usr/bin/env bash
set -euo pipefail

FILE="examples/NMFS/sefsc_red_snapper/tmb/evaluate_tmb_at_quadra_fit.R"

if [[ ! -f "$FILE" ]]; then
  echo "ERROR: $FILE not found. Run from Quadra repo root."
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
BACKUP="${FILE}.before_manual_random_profiled_logdet_fd.${STAMP}"
cp "$FILE" "$BACKUP"
echo "Backed up $FILE to:"
echo "  $BACKUP"

python3 - <<'PY'
from pathlib import Path

path = Path("examples/NMFS/sefsc_red_snapper/tmb/evaluate_tmb_at_quadra_fit.R")
text = path.read_text()

if "TMB manual-random-optimized profiled logdet FD" in text:
    print("Manual profiled logdet FD diagnostic already installed.")
    raise SystemExit(0)

append = r'''

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
'''

text = text + append
path.write_text(text)
print("Installed manual-random-optimized profiled logdet FD diagnostic.")
PY

echo
echo "Run:"
echo "Rscript $FILE"
