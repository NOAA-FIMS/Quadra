#!/usr/bin/env bash
set -euo pipefail

FILE="examples/NMFS/sefsc_red_snapper/tmb/evaluate_tmb_at_quadra_fit.R"

if [[ ! -f "$FILE" ]]; then
  echo "ERROR: $FILE not found. Run from Quadra repo root."
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
BACKUP="${FILE}.before_profiled_logdet_fd.${STAMP}"
cp "$FILE" "$BACKUP"
echo "Backed up $FILE to:"
echo "  $BACKUP"

python3 - <<'PY'
from pathlib import Path

path = Path("examples/NMFS/sefsc_red_snapper/tmb/evaluate_tmb_at_quadra_fit.R")
text = path.read_text()

if "TMB profiled logdet FD at Quadra fit" in text:
    print("Profiled logdet FD diagnostic already installed.")
    raise SystemExit(0)

append = r'''

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
'''

text = text + append
path.write_text(text)
print("Installed TMB profiled logdet FD diagnostic.")
PY

echo
echo "Run:"
echo "Rscript $FILE"
