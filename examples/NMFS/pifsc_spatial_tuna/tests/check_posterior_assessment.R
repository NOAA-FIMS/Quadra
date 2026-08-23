args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1L) args[[1L]] else
  "build/assessment_outputs/data"
report_dir <- if (length(args) >= 2L) args[[2L]] else
  "build/assessment_outputs/report"

fail <- function(...) stop(..., call. = FALSE)
require_file <- function(name) {
  path <- file.path(data_dir, name)
  if (!file.exists(path) || file.info(path)$size <= 0L) {
    fail("posterior assessment check: missing or empty ", path)
  }
  path
}
require_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing) > 0L) {
    fail(label, " missing columns: ", paste(missing, collapse = ", "))
  }
}
metric <- function(x, name) {
  value <- x$value[x$metric == name]
  if (length(value) != 1L) fail("missing summary metric: ", name)
  value[[1L]]
}
numeric_metric <- function(x, name) {
  value <- suppressWarnings(as.numeric(metric(x, name)))
  if (!is.finite(value)) fail("non-finite summary metric: ", name)
  value
}

summary <- read.csv(require_file("sampler_summary.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)
reconstruction <- read.csv(
  require_file("posterior_reconstruction_summary.csv"),
  stringsAsFactors = FALSE, check.names = FALSE)
random_draws <- read.csv(require_file("posterior_random_effect_draws.csv"),
                         check.names = FALSE)
references <- read.csv(require_file("posterior_reference_points.csv"),
                       check.names = FALSE)
projections <- read.csv(require_file("posterior_projection_draws.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)

if (metric(summary, "health") != "PASS") fail("sampler health is not PASS")
draw_count <- as.integer(numeric_metric(summary, "draws"))
if (numeric_metric(summary, "max_rhat") >
    numeric_metric(summary, "rhat_threshold")) fail("R-hat threshold failed")
if (numeric_metric(summary, "min_bulk_ess") <
    numeric_metric(summary, "bulk_ess_threshold")) fail("bulk ESS threshold failed")
if (numeric_metric(summary, "min_tail_ess") <
    numeric_metric(summary, "tail_ess_threshold")) fail("tail ESS threshold failed")
if (numeric_metric(summary, "divergences") != 0) fail("divergences detected")

if (numeric_metric(reconstruction, "retained_draws") != draw_count ||
    numeric_metric(reconstruction, "valid_reconstructions") != draw_count ||
    numeric_metric(reconstruction, "failed_reconstructions") != 0) {
  fail("not every retained draw has a valid latent reconstruction")
}
require_columns(random_draws, c("chain", "iteration", "valid"),
                "random-effect draws")
if (nrow(random_draws) != draw_count || any(random_draws$valid != 1L)) {
  fail("random-effect draw count or validity does not match sampler output")
}
random_values <- as.matrix(random_draws[, -(1:3), drop = FALSE])
if (ncol(random_values) == 0L || any(!is.finite(random_values))) {
  fail("random-effect reconstruction contains missing or non-finite values")
}

management_draws <- as.integer(numeric_metric(reconstruction, "management_draws"))
require_columns(
  references,
  c("chain", "iteration", "valid", "B0", "B_MSY", "MSY",
    "F_MSY_multiplier", "B_terminal_over_B_MSY",
    "F_status_quo_over_F_MSY"),
  "posterior reference points")
if (nrow(references) != management_draws || any(references$valid != 1L)) {
  fail("reference-point count or validity does not match management draws")
}
reference_values <- as.matrix(references[, -(1:3), drop = FALSE])
if (any(!is.finite(reference_values)) ||
    any(references$B0 <= 0) || any(references$B_MSY <= 0) ||
    any(references$MSY <= 0) || any(references$F_MSY_multiplier < 0) ||
    any(references$B_terminal_over_B_MSY <= 0) ||
    any(references$F_status_quo_over_F_MSY < 0)) {
  fail("posterior reference points violate finite/positive constraints")
}

require_columns(
  projections,
  c("chain", "iteration", "scenario", "projection_year",
    "fishing_multiplier", "spawning_biomass", "depletion",
    "retained_yield", "discard_yield", "total_yield"),
  "posterior projections")
expected_scenarios <- c("F_MSY", "half_status_quo", "no_fishing", "status_quo")
if (!identical(sort(unique(projections$scenario)), expected_scenarios)) {
  fail("posterior projections do not contain the four required scenarios")
}
expected_projection_rows <- management_draws * length(expected_scenarios) * 10L
if (nrow(projections) != expected_projection_rows) {
  fail("posterior projection row count is ", nrow(projections),
       "; expected ", expected_projection_rows)
}
projection_values <- projections[, c(
  "projection_year", "fishing_multiplier", "spawning_biomass", "depletion",
  "retained_yield", "discard_yield", "total_yield"), drop = FALSE]
if (any(!is.finite(as.matrix(projection_values))) ||
    any(projections$fishing_multiplier < 0) ||
    any(projections$spawning_biomass <= 0) || any(projections$depletion <= 0) ||
    any(projections$retained_yield < 0) || any(projections$discard_yield < 0) ||
    any(projections$total_yield < 0)) {
  fail("posterior projections violate finite/non-negative constraints")
}

if (dir.exists(report_dir)) {
  data_path <- normalizePath(data_dir, mustWork = TRUE)
  report_path <- normalizePath(report_dir, mustWork = TRUE)
  if (identical(data_path, report_path)) fail("data and report directories coincide")
  report_csv <- list.files(report_dir, pattern = "\\.csv$", recursive = TRUE)
  if (length(report_csv) > 0L) {
    fail("CSV files found in report directory: ", paste(report_csv, collapse = ", "))
  }
}

cat(sprintf(
  paste0("posterior assessment check: PASS (%d draws, %d reconstructions, ",
         "%d management draws, %d projections)\n"),
  draw_count, nrow(random_draws), management_draws, nrow(projections)))
