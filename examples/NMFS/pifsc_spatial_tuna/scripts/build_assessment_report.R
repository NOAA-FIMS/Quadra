args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1L) args[[1L]] else "build/assessment_outputs/data"
report_dir <- if (length(args) >= 2L) args[[2L]] else "build/assessment_outputs/report"
report_path <- file.path(report_dir, "assessment_report.md")
figure_dir <- file.path(report_dir, "figures")

required_packages <- c("jsonlite", "ggplot2", "svglite")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0L) {
  stop("assessment report requires R packages: ",
       paste(missing_packages, collapse = ", "))
}

required_files <- c(
  "acceptance_summary.json",
  "residual_summary.csv",
  "likelihood_decomposition.csv",
  "parameter_diagnostics.csv"
)
missing_files <- required_files[!file.exists(file.path(data_dir, required_files))]
if (length(missing_files) > 0L) {
  stop("assessment report missing required outputs: ",
       paste(missing_files, collapse = ", "))
}

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
summary <- jsonlite::fromJSON(file.path(data_dir, "acceptance_summary.json"),
                              simplifyVector = TRUE)
residuals <- read.csv(file.path(data_dir, "residual_summary.csv"),
                      check.names = FALSE)
likelihood <- read.csv(file.path(data_dir, "likelihood_decomposition.csv"),
                       check.names = FALSE)
parameters <- read.csv(file.path(data_dir, "parameter_diagnostics.csv"),
                       check.names = FALSE)
biomass_path <- file.path(data_dir, "biomass_trajectory.csv")
biomass <- if (file.exists(biomass_path)) read.csv(biomass_path) else NULL
reference_path <- file.path(data_dir, "reference_points.csv")
reference_metadata_path <- file.path(data_dir, "reference_point_metadata.csv")
projection_path <- file.path(data_dir, "projection_summary.csv")
reference <- if (file.exists(reference_path)) read.csv(reference_path) else NULL
reference_metadata <- if (file.exists(reference_metadata_path)) {
  read.csv(reference_metadata_path, check.names = FALSE)
} else NULL
projections <- if (file.exists(projection_path)) read.csv(projection_path) else NULL
spatial_path <- file.path(data_dir, "spatial_animation.csv")
spatial_section <- if (file.exists(spatial_path)) {
  paste0(
    "[Open the interactive spatial fishery pulse map]",
    "(spatial_fishery_pulse.html). It animates seasonal regional biomass, ",
    "movement, fleet-specific retained catch and discards, and the projected ",
    "management-scenario race."
  )
} else {
  "Spatial animation data are unavailable for this run."
}
nuts_path <- file.path(data_dir, "sampler_summary.csv")
if (!file.exists(nuts_path)) {
  nuts_path <- file.path(data_dir, "nuts_summary.csv")
}
nuts <- if (file.exists(nuts_path)) read.csv(nuts_path, stringsAsFactors = FALSE) else NULL
sampler_identity_path <- file.path(data_dir, "sampler_identity.csv")
sampler_identity <- if (file.exists(sampler_identity_path)) {
  read.csv(sampler_identity_path, stringsAsFactors = FALSE)
} else NULL
sampler_name <- if (!is.null(sampler_identity)) {
  method <- sampler_identity$value[sampler_identity$metric == "method"]
  if (length(method) == 0L) "sampler" else method[[1L]]
} else "AD-NUTS"
posterior_reference_path <- file.path(data_dir, "posterior_reference_points.csv")
posterior_projection_path <- file.path(data_dir, "posterior_projection_draws.csv")
posterior_reference <- if (file.exists(posterior_reference_path)) {
  read.csv(posterior_reference_path, check.names = FALSE)
} else NULL
posterior_projections <- if (file.exists(posterior_projection_path)) {
  read.csv(posterior_projection_path, check.names = FALSE)
} else NULL

assert_columns <- function(data, columns, label) {
  absent <- setdiff(columns, names(data))
  if (length(absent) > 0L) {
    stop(label, " missing columns: ", paste(absent, collapse = ", "))
  }
}
assert_columns(residuals, c("component", "fleet", "sdnr"), "residual summary")
assert_columns(likelihood, c("component", "nll", "absolute_share"),
               "likelihood decomposition")
assert_columns(parameters, c("name", "value", "displacement", "random_effect"),
               "parameter diagnostics")

theme_assessment <- ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.title.position = "plot",
    legend.position = "bottom"
  )

save_svg <- function(filename, plot, width = 7.2, height = 4.4) {
  ggplot2::ggsave(
    filename = file.path(figure_dir, filename), plot = plot,
    device = svglite::svglite, width = width, height = height,
    units = "in", bg = "white"
  )
}

residuals$fleet <- factor(residuals$fleet)
p_residual <- ggplot2::ggplot(
  residuals,
  ggplot2::aes(x = component, y = sdnr, fill = fleet)
) +
  ggplot2::geom_hline(yintercept = 1, colour = "#333333", linewidth = 0.5) +
  ggplot2::geom_hline(
    yintercept = c(0.85, 1.15), colour = "#b2182b",
    linewidth = 0.45, linetype = "dashed"
  ) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.72) +
  ggplot2::coord_cartesian(ylim = c(0, max(1.25, residuals$sdnr, na.rm = TRUE))) +
  ggplot2::labs(
    title = "Residual dispersion by data component",
    x = NULL, y = "SDNR", fill = "Fleet"
  ) + theme_assessment
save_svg("residual_sdnr.svg", p_residual)

p_likelihood <- ggplot2::ggplot(
  likelihood,
  ggplot2::aes(x = reorder(component, absolute_share), y = absolute_share)
) +
  ggplot2::geom_col(fill = "#2166ac", width = 0.72) +
  ggplot2::coord_flip() +
  ggplot2::scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100 * x)) +
  ggplot2::labs(
    title = "Absolute contribution to the objective",
    x = NULL, y = "Absolute share"
  ) + theme_assessment
save_svg("likelihood_components.svg", p_likelihood)

fixed_parameters <- parameters[parameters$random_effect == 0, , drop = FALSE]
fixed_parameters <- fixed_parameters[order(abs(fixed_parameters$displacement),
                                           decreasing = TRUE), , drop = FALSE]
fixed_parameters <- head(fixed_parameters, 15L)
p_parameters <- ggplot2::ggplot(
  fixed_parameters,
  ggplot2::aes(x = reorder(name, displacement), y = displacement,
               fill = displacement > 0)
) +
  ggplot2::geom_hline(yintercept = 0, colour = "#333333", linewidth = 0.4) +
  ggplot2::geom_col(width = 0.72, show.legend = FALSE) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = c("#b2182b", "#2166ac")) +
  ggplot2::labs(
    title = "Largest fixed-parameter movements from initial values",
    x = NULL, y = "Estimate minus initial value"
  ) + theme_assessment
save_svg("parameter_displacement.svg", p_parameters, height = 5.8)

biomass_section <- "Historical biomass trajectory is unavailable."
if (!is.null(biomass)) {
  assert_columns(biomass, c("year", "spawning_biomass", "depletion"),
                 "biomass trajectory")
  p_biomass <- ggplot2::ggplot(
    biomass, ggplot2::aes(x = year, y = depletion)
  ) +
    ggplot2::geom_hline(yintercept = 1, colour = "#555555", linewidth = 0.45,
                        linetype = "dashed") +
    ggplot2::geom_line(colour = "#2166ac", linewidth = 0.9) +
    ggplot2::geom_point(colour = "#2166ac", size = 1.7) +
    ggplot2::labs(title = "Spawning-biomass trajectory",
                  x = "Model year", y = "SSB / SSB0") + theme_assessment
  save_svg("spawning_biomass.svg", p_biomass)
  biomass_section <- paste0(
    "![Spawning-biomass trajectory](figures/spawning_biomass.svg)\n\n",
    "This is depletion relative to the model's unfished spawning biomass, ",
    "not B/BMSY."
  )
}

fmt <- function(x, digits = 4L) {
  if (is.null(x) || length(x) == 0L) return("NA")
  formatted <- formatC(x, digits = digits, format = "fg", flag = "#")
  formatted[!is.finite(x)] <- "NA"
  formatted
}

fit <- summary$fit
converged <- isTRUE(fit$converged)

reference_section <- paste(
  "Reference-point calculations are unavailable.",
  "BMSY, FMSY, MSY, B/BMSY, and F/FMSY are not inferred from depletion."
)
projection_section <- paste(
  "Projection results are unavailable.",
  "No management projection is reported from an unconverged or invalid fit."
)
uncertainty_section <- paste(
  "Posterior sampling was not requested for this run.",
  "The Laplace fit remains the primary assessment estimate."
)
posterior_projection_section <- NULL
if (!is.null(nuts)) {
  assert_columns(nuts, c("metric", "value"), "sampler summary")
  nuts_value <- function(metric) {
    value <- nuts$value[nuts$metric == metric]
    if (length(value) == 0L) "NA" else value[[1L]]
  }
  nuts_health <- nuts_value("health")
  uncertainty_section <- paste0(
    if (identical(nuts_health, "PASS")) "" else
      paste0("> [!CAUTION]\n> ", sampler_name,
             " health checks failed; posterior summaries must not be used.\n\n"),
    "| Sampler diagnostic | Value |\n|---|---:|\n",
    "| Health | **", nuts_health, "** |\n",
    "| Chains | ", nuts_value("chains"), " |\n",
    "| Total draws | ", nuts_value("draws"), " |\n",
    "| Maximum R-hat | ", nuts_value("max_rhat"), " |\n",
    "| Minimum bulk ESS | ", nuts_value("min_bulk_ess"), " |\n",
    "| Minimum tail ESS | ", nuts_value("min_tail_ess"), " |\n",
    "| Divergences | ", nuts_value("divergences"), " |"
  )
  if (identical(nuts_health, "PASS") && converged &&
      !is.null(posterior_reference)) {
    valid_reference <- posterior_reference[posterior_reference$valid == 1, ,
                                           drop = FALSE]
    management_metrics <- c("B0", "B_MSY", "MSY", "F_MSY_multiplier",
                            "B_terminal_over_B_MSY",
                            "F_status_quo_over_F_MSY")
    if (nrow(valid_reference) > 1L &&
        all(management_metrics %in% names(valid_reference))) {
      interval_rows <- vapply(management_metrics, function(metric) {
        values <- valid_reference[[metric]]
        interval <- stats::quantile(values, c(0.05, 0.5, 0.95), na.rm = TRUE,
                                    names = FALSE)
        paste0("| ", metric, " | ", fmt(interval[[2L]]), " | ",
               fmt(interval[[1L]]), " | ", fmt(interval[[3L]]), " |")
      }, character(1))
      uncertainty_section <- paste0(
        uncertainty_section, "\n\n### Posterior management quantities\n\n",
        "Intervals combine parameter uncertainty and fitted recruitment effects.\n\n",
        "| Quantity | Median | 5% | 95% |\n|---|---:|---:|---:|\n",
        paste(interval_rows, collapse = "\n"), "\n\n",
        "Valid management draws: ", nrow(valid_reference), " of ",
        nrow(posterior_reference), "."
      )
      p_kobe <- ggplot2::ggplot(
        valid_reference,
        ggplot2::aes(x = F_status_quo_over_F_MSY,
                     y = B_terminal_over_B_MSY)
      ) +
        ggplot2::annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = Inf,
                          fill = "#1a9850", alpha = 0.16) +
        ggplot2::annotate("rect", xmin = 1, xmax = Inf, ymin = 1, ymax = Inf,
                          fill = "#fee08b", alpha = 0.18) +
        ggplot2::annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                          fill = "#fee08b", alpha = 0.18) +
        ggplot2::annotate("rect", xmin = 1, xmax = Inf, ymin = 0, ymax = 1,
                          fill = "#d73027", alpha = 0.16) +
        ggplot2::geom_vline(xintercept = 1, linewidth = 0.5) +
        ggplot2::geom_hline(yintercept = 1, linewidth = 0.5) +
        ggplot2::geom_point(colour = "#252525", alpha = 0.5, size = 1.5) +
        ggplot2::labs(title = "Posterior terminal stock status",
                      x = "F / FMSY", y = "B / BMSY") + theme_assessment
      save_svg("posterior_kobe.svg", p_kobe, width = 6.2, height = 5.2)
      uncertainty_section <- paste0(
        uncertainty_section, "\n\n",
        "![Posterior Kobe status](figures/posterior_kobe.svg)"
      )
    }
  }
  if (identical(nuts_health, "PASS") && converged &&
      !is.null(posterior_projections) && nrow(posterior_projections) > 1L) {
    groups <- split(
      posterior_projections,
      interaction(posterior_projections$scenario,
                  posterior_projections$projection_year, drop = TRUE)
    )
    projection_intervals <- do.call(rbind, lapply(groups, function(group) {
      interval <- stats::quantile(group$depletion, c(0.05, 0.5, 0.95),
                                  na.rm = TRUE, names = FALSE)
      data.frame(scenario = group$scenario[[1L]],
                 projection_year = group$projection_year[[1L]],
                 lower = interval[[1L]], median = interval[[2L]],
                 upper = interval[[3L]])
    }))
    p_posterior_projection <- ggplot2::ggplot(
      projection_intervals,
      ggplot2::aes(x = projection_year, y = median, colour = scenario,
                   fill = scenario, group = scenario)
    ) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                           alpha = 0.16, colour = NA) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::labs(title = "Posterior biomass projections",
                    subtitle = "Median and 90% posterior predictive interval",
                    x = "Projection year", y = "SSB / SSB0",
                    colour = "Scenario", fill = "Scenario") + theme_assessment
    save_svg("posterior_biomass_projections.svg", p_posterior_projection)
    posterior_projection_section <- paste0(
      "![Posterior biomass projections]",
      "(figures/posterior_biomass_projections.svg)\n\n",
      "The ribbons propagate joint posterior parameter uncertainty and ",
      "lognormal recruitment variability. Scenario trajectories share each ",
      "draw's recruitment realization for paired comparison."
    )
    if (!is.null(posterior_reference)) {
      valid_reference <- posterior_reference[posterior_reference$valid == 1,
                                             drop = FALSE]
      joined <- merge(
        posterior_projections, valid_reference,
        by = c("chain", "iteration"), suffixes = c("_projection", "_reference")
      )
      if (nrow(joined) > 0L) {
        joined$B_over_BMSY <- joined$spawning_biomass / joined$B_MSY
        joined$F_over_FMSY <- joined$fishing_multiplier /
          joined$F_MSY_multiplier
        decision_years <- joined[joined$projection_year %in% c(5, 10), ,
                                 drop = FALSE]
        decision_groups <- split(
          decision_years,
          interaction(decision_years$scenario, decision_years$projection_year,
                      drop = TRUE)
        )
        decision_rows <- vapply(decision_groups, function(group) {
          paste0(
            "| ", group$scenario[[1L]], " | ",
            group$projection_year[[1L]], " | ",
            fmt(mean(group$B_over_BMSY > 1)), " | ",
            fmt(mean(group$F_over_FMSY < 1)), " | ",
            fmt(mean(group$depletion > 0.2)), " |"
          )
        }, character(1))
        posterior_projection_section <- paste0(
          posterior_projection_section,
          "\n\n### Probabilistic decision table\n\n",
          "| Scenario | Year | P(B > BMSY) | P(F < FMSY) | P(B > 0.2 B0) |\n",
          "|---|---:|---:|---:|---:|\n",
          paste(decision_rows, collapse = "\n")
        )
      }
    }
  }
}
if (!is.null(reference) && !is.null(reference_metadata)) {
  assert_columns(reference, c("metric", "value", "valid"), "reference points")
  assert_columns(reference_metadata,
                 c("valid", "grid_boundary", "fishing_pattern_source", "message"),
                 "reference-point metadata")
  reference_valid <- nrow(reference_metadata) == 1L &&
    reference_metadata$valid[[1L]] == 1L && converged
  if (reference_valid) {
    reference_rows <- paste0(
      "| ", reference$metric, " | ", fmt(reference$value), " |",
      collapse = "\n"
    )
    reference_section <- paste0(
      "Fishing pattern: `", reference_metadata$fishing_pattern_source[[1L]],
      "`. F quantities are multipliers on that pattern.\n\n",
      "| Quantity | Estimate |\n|---|---:|\n", reference_rows
    )
    if (!is.null(projections) && nrow(projections) > 0L) {
      assert_columns(projections,
                     c("scenario", "projection_year", "depletion", "total_yield"),
                     "projection summary")
      p_projection <- ggplot2::ggplot(
        projections,
        ggplot2::aes(x = projection_year, y = depletion,
                     colour = scenario, group = scenario)
      ) +
        ggplot2::geom_line(linewidth = 0.9) +
        ggplot2::geom_point(size = 1.4) +
        ggplot2::labs(title = "Deterministic biomass projections",
                      x = "Projection year", y = "SSB / SSB0",
                      colour = "Scenario") + theme_assessment
      save_svg("biomass_projections.svg", p_projection)
      terminal_projection <- projections[
        !duplicated(projections$scenario, fromLast = TRUE), , drop = FALSE
      ]
      projection_rows <- paste0(
        "| ", terminal_projection$scenario, " | ",
        fmt(terminal_projection$fishing_multiplier), " | ",
        fmt(terminal_projection$depletion), " | ",
        fmt(terminal_projection$total_yield), " |", collapse = "\n"
      )
      projection_section <- paste0(
        "![Deterministic biomass projections](figures/biomass_projections.svg)\n\n",
        "Ten-year deterministic projections use mean Beverton–Holt recruitment; ",
        "they do not yet represent recruitment or parameter uncertainty.\n\n",
        "| Scenario | Fishing multiplier | Year-10 SSB/SSB0 | Year-10 yield |\n",
        "|---|---:|---:|---:|\n", projection_rows
      )
    }
  } else {
    reference_section <- paste0(
      "> [!CAUTION]\n> Reference points withheld: ",
      reference_metadata$message[[1L]], "."
    )
  }
}

if (!is.null(posterior_projection_section))
  projection_section <- posterior_projection_section

status <- if (converged) "PASS" else "FAIL — NOT CONVERGED"
warning_block <- if (converged) {
  "> **Assessment status:** optimizer convergence checks passed."
} else {
  paste0(
    "> [!CAUTION]\n",
    "> **Assessment status: NOT CONVERGED.** Estimates, stock status, reference ",
    "points, and projections must not be used for management advice.\n",
    "> Final gradient norm: `", fmt(fit$gradient_norm), "`."
  )
}

largest <- parameters[order(abs(parameters$displacement), decreasing = TRUE),
                      c("name", "value", "displacement"), drop = FALSE]
largest <- head(largest, 10L)
parameter_rows <- paste0(
  "| `", largest$name, "` | ", fmt(largest$value), " | ",
  fmt(largest$displacement), " |",
  collapse = "\n"
)

residual_rows <- paste0(
  "| ", residuals$component, " | ", residuals$fleet, " | ", residuals$n,
  " | ", fmt(residuals$mean_residual), " | ", fmt(residuals$sdnr), " |",
  collapse = "\n"
)

report <- paste0(
  "# Tuna Stock Assessment Report\n\n",
  warning_block, "\n\n",
  "## Executive summary\n\n",
  "| Quantity | Value |\n|---|---:|\n",
  "| Assessment status | **", status, "** |\n",
  "| Objective (NLL) | ", fmt(fit$nll), " |\n",
  "| Laplace gradient norm | ", fmt(fit$gradient_norm), " |\n",
  "| Total optimizer iterations | ", fit$total_iterations, " |\n",
  "| Terminal depletion (SSB/SSB0) | ", fmt(fit$depletion_terminal), " |\n\n",
  "The terminal depletion value is a model-derived spawning-biomass ratio. ",
  "It is not equivalent to B/BMSY until equilibrium reference points have been ",
  "calculated and validated.\n\n",
  "## Reference points and stock status\n\n",
  reference_section, "\n\n",
  "## Projections\n\n",
  projection_section, "\n\n",
  "## Model fit and diagnostics\n\n",
  "### Population trajectory\n\n", biomass_section, "\n\n",
  "### Spatial fishery pulse\n\n", spatial_section, "\n\n",
  "### Residual diagnostics\n\n",
  "![Residual dispersion](figures/residual_sdnr.svg)\n\n",
  "The dashed lines show the current SDNR screening interval of 0.85–1.15.\n\n",
  "| Component | Fleet | N | Mean residual | SDNR |\n",
  "|---|---:|---:|---:|---:|\n", residual_rows, "\n\n",
  "![Likelihood components](figures/likelihood_components.svg)\n\n",
  "## Parameter diagnostics\n\n",
  "![Parameter displacement](figures/parameter_displacement.svg)\n\n",
  "| Parameter | Estimate | Displacement |\n",
  "|---|---:|---:|\n", parameter_rows, "\n\n",
  "## Uncertainty and sampler diagnostics\n\n",
  uncertainty_section, "\n\n",
  "## Reproducibility\n\n",
  "This report was generated directly from the machine-readable files in this ",
  "directory. Missing or non-finite JSON values are represented as `null`. ",
  "Figures are SVG files linked with relative paths so the report remains ",
  "portable.\n"
)

writeLines(report, report_path, useBytes = TRUE)
cat("Assessment report written to", report_path, "\n")

if (file.exists(spatial_path)) {
  source(file.path("scripts", "build_spatial_pulse_map.R"), local = new.env())
}
