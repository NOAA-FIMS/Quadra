args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1L) args[[1L]] else "data/opal_raw"
output_file <- if (length(args) >= 2L) args[[2L]] else "build/opal_data_audit.csv"

load(file.path(input_dir, "wcpo_bet_data.rda"))
load(file.path(input_dir, "wcpo_bet_lf.rda"))
load(file.path(input_dir, "wcpo_bet_wf.rda"))
load(file.path(input_dir, "wcpo_bet_parameters.rda"))

d <- wcpo_bet_data
rows <- data.frame(
  item = character(), value = character(), status = character(),
  stringsAsFactors = FALSE
)
add <- function(item, value, status = "ok") {
  rows <<- rbind(rows, data.frame(
    item = item, value = as.character(value), status = status,
    stringsAsFactors = FALSE
  ))
}

add("quarterly_timesteps", d$n_year)
add("calendar_years", length(unique(d$cpue_data$year)))
add("fleets", d$n_fishery)
add("ages", d$n_age)
add("length_bins", d$n_len)
add("catch_positive_cells", sum(d$catch_obs_ysf > 0, na.rm = TRUE))
add("catch_missing_cells", sum(is.na(d$catch_obs_ysf)))
add("catch_units", paste(sort(unique(d$catch_units_f)), collapse = ";"),
    "conversion_required")
add("cpue_rows", nrow(d$cpue_data))
add("cpue_log_se_range", paste(range(d$cpue_data$se), collapse = ";"))
add("length_frequency_rows", nrow(wcpo_bet_lf))
add("length_frequency_instances",
    nrow(unique(wcpo_bet_lf[c("year", "month", "fishery", "week")])))
add("weight_frequency_rows", nrow(wcpo_bet_wf))
add("weight_frequency_instances",
    nrow(unique(wcpo_bet_wf[c("year", "month", "fishery", "week")])))
add("all_source_missing_values",
    sum(is.na(d$catch_obs_ysf)) + sum(is.na(d$cpue_data)) +
      sum(is.na(wcpo_bet_lf)) + sum(is.na(wcpo_bet_wf)))
add("explicit_regions", 0, "single_region_assumption_required")
add("effort_observations", 0, "resolved_by_catch_conditioning")
add("discard_observations", 0, "discard_likelihood_must_be_disabled")
add("age_compositions", 0, "age_length_conversion_required")

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.csv(rows, output_file, row.names = FALSE, quote = TRUE)
cat("Opal WCPO bigeye audit written to", output_file, "\n")
print(rows, row.names = FALSE)
