#ifndef QUADRA_TUNA_REFERENCE_POINTS_HPP
#define QUADRA_TUNA_REFERENCE_POINTS_HPP
#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "tuna_spatial_assessment_model.hpp"

namespace quadra
{
    struct TunaFishingPattern
    {
        // Seasonal fishing rate, [fleet, season, region, age]. In an
        // effort-driven model this is instantaneous F; in a catch-conditioned
        // model it is the fitted removal fraction used by the assessment.
        std::vector<double> rate_m;
        std::string source_m;
        bool removal_fraction_m = false;
    };

    struct TunaProjectionPoint
    {
        int projection_year_m = 0;
        double fishing_multiplier_m = 0.0;
        double spawning_biomass_m = std::numeric_limits<double>::quiet_NaN();
        double depletion_m = std::numeric_limits<double>::quiet_NaN();
        double retained_yield_m = 0.0;
        double discard_yield_m = 0.0;
        double total_yield_m = 0.0;
    };

    struct TunaReferencePoints
    {
        bool valid_m = false;
        std::string message_m;
        std::string fishing_pattern_source_m;
        double b0_m = std::numeric_limits<double>::quiet_NaN();
        double b_msy_m = std::numeric_limits<double>::quiet_NaN();
        double msy_m = std::numeric_limits<double>::quiet_NaN();
        double f_msy_multiplier_m = std::numeric_limits<double>::quiet_NaN();
        double terminal_b_over_b_msy_m = std::numeric_limits<double>::quiet_NaN();
        double status_quo_f_over_f_msy_m = std::numeric_limits<double>::quiet_NaN();
        bool grid_boundary_m = false;
    };

    namespace tuna_reference_detail
    {
        struct Dynamics
        {
            double r0 = 0.0;
            double steepness = 0.0;
            double sigma_recruit = 0.0;
            std::vector<double> selectivity;
            std::vector<double> retention;
            std::vector<double> availability;
            std::vector<double> movement;
        };

        inline size_t fsra(const TunaSpatialAssessmentData &data, int f, int s,
                           int r, int a)
        {
            return data.fleet_season_region_age_index(f, s, r, a);
        }

        inline size_t srr(const TunaSpatialAssessmentData &data, int s, int from,
                          int to)
        {
            return data.season_region_region_index(s, from, to);
        }

        inline double report_value(const ModelReportContext &ctx,
                                   const std::string &name)
        {
            for (const auto &value : ctx.reports().values())
            {
                if (value.name_m == name)
                    return value.value_m;
            }
            throw std::runtime_error("missing model report value: " + name);
        }

        inline Dynamics unpack_dynamics(const TunaSpatialAssessmentData &data,
                                        const TunaAssessmentControls &controls,
                                        const std::vector<double> &parameters)
        {
            AdvancedSpatialTunaAssessmentModel model(data, controls);
            if (parameters.size() != model.parameter_set().size())
                throw std::invalid_argument(
                    "unpack_dynamics: parameter vector length mismatch");

            Dynamics out;
            size_t pos = 0;
            out.r0 = std::exp(parameters[pos++]);
            out.steepness = 0.2 + 0.8 / (1.0 + std::exp(-parameters[pos++]));
            out.selectivity.assign(
                static_cast<size_t>(data.n_fleets_m * data.n_ages_m), 0.0);
            out.retention.assign(out.selectivity.size(), 0.0);
            out.availability.assign(
                static_cast<size_t>(data.n_fleets_m * data.n_seasons_m *
                                    data.n_regions_m * data.n_ages_m),
                1.0);
            std::vector<double> log_q(static_cast<size_t>(data.n_fleets_m));

            for (int f = 0; f < data.n_fleets_m; ++f)
            {
                log_q[static_cast<size_t>(f)] = parameters[pos++];
                ++pos; // index q
                const double sel_raw = parameters[pos++];
                const double sel_slope = std::exp(parameters[pos++]);
                const double ret_raw = parameters[pos++];
                const double ret_slope = std::exp(parameters[pos++]);
                pos += 4; // observation scale and composition parameters

                const double sel50 = 1.0 + (data.n_ages_m - 1.0) /
                                               (1.0 + std::exp(-sel_raw));
                const double ret50 = 1.0 + (data.n_ages_m - 1.0) /
                                               (1.0 + std::exp(-ret_raw));
                const auto logistic = [](double age, double a50, double slope)
                {
                    return 1.0 / (1.0 + std::exp(-slope * (age - a50)));
                };
                const double oldest_sel =
                    logistic(static_cast<double>(data.n_ages_m), sel50, sel_slope);
                for (int a = 0; a < data.n_ages_m; ++a)
                {
                    const size_t fa = static_cast<size_t>(f * data.n_ages_m + a);
                    out.selectivity[fa] =
                        logistic(static_cast<double>(a + 1), sel50, sel_slope) /
                        (oldest_sel + 1e-12);
                    out.retention[fa] =
                        logistic(static_cast<double>(a + 1), ret50, ret_slope);
                }

                std::vector<double> scale(
                    static_cast<size_t>(data.n_seasons_m * data.n_regions_m), 1.0);
                if (controls.estimate_availability_scales_m)
                {
                    if (controls.availability_by_fleet_only_m)
                    {
                        const double value = std::exp(parameters[pos++]);
                        std::fill(scale.begin(), scale.end(), value);
                    }
                    else
                    {
                        for (double &value : scale)
                            value = std::exp(parameters[pos++]);
                    }
                }
                for (int s = 0; s < data.n_seasons_m; ++s)
                    for (int r = 0; r < data.n_regions_m; ++r)
                        for (int a = 0; a < data.n_ages_m; ++a)
                        {
                            const double surface = data.availability_surface_m.empty()
                                ? 1.0
                                : data.availability_surface_m[fsra(data, f, s, r, a)];
                            out.availability[fsra(data, f, s, r, a)] =
                                surface * scale[static_cast<size_t>(
                                    s * data.n_regions_m + r)];
                        }
            }

            const int movement_seasons =
                controls.share_movement_across_seasons_m ? 1 : data.n_seasons_m;
            std::vector<double> logits(static_cast<size_t>(
                movement_seasons * data.n_regions_m * (data.n_regions_m - 1)));
            for (double &value : logits)
                value = parameters[pos++];
            out.movement.assign(static_cast<size_t>(
                data.n_seasons_m * data.n_regions_m * data.n_regions_m), 0.0);
            for (int s = 0; s < data.n_seasons_m; ++s)
                for (int from = 0; from < data.n_regions_m; ++from)
                {
                    const int ps = controls.share_movement_across_seasons_m ? 0 : s;
                    double denominator = 1.0;
                    for (int to = 0; to < data.n_regions_m - 1; ++to)
                        denominator += std::exp(logits[static_cast<size_t>(
                            (ps * data.n_regions_m + from) *
                                (data.n_regions_m - 1) + to)]);
                    for (int to = 0; to < data.n_regions_m - 1; ++to)
                        out.movement[srr(data, s, from, to)] =
                            std::exp(logits[static_cast<size_t>(
                                (ps * data.n_regions_m + from) *
                                    (data.n_regions_m - 1) + to)]) /
                            denominator;
                    out.movement[srr(data, s, from, data.n_regions_m - 1)] =
                        1.0 / denominator;
                }
            out.sigma_recruit = std::exp(parameters[pos++]);
            return out;
        }

        inline std::vector<double> terminal_state(
            const TunaSpatialAssessmentData &data,
            const TunaAssessmentControls &controls,
            const std::vector<double> &parameters,
            ModelReportContext *context_out = nullptr)
        {
            TunaAssessmentControls diagnostic_controls = controls;
            diagnostic_controls.report_observation_predictions_m = true;
            AdvancedSpatialTunaAssessmentModel model(data, diagnostic_controls);
            ModelReportContext context;
            model.evaluate(parameters, context);
            std::vector<double> state(static_cast<size_t>(
                data.n_regions_m * data.n_ages_m));
            for (int r = 0; r < data.n_regions_m; ++r)
                for (int a = 0; a < data.n_ages_m; ++a)
                    state[data.region_age_index(r, a)] = report_value(
                        context, "terminal_numbers_region_" +
                                     std::to_string(r + 1) + "_age_" +
                                     std::to_string(a + 1));
            if (context_out)
                *context_out = std::move(context);
            return state;
        }

        inline double spawning_biomass(const TunaSpatialAssessmentData &data,
                                       const std::vector<double> &state)
        {
            double out = 0.0;
            const int y = data.n_years_m - 1;
            for (int r = 0; r < data.n_regions_m; ++r)
                for (int a = 0; a < data.n_ages_m; ++a)
                    out += state[data.region_age_index(r, a)] *
                           data.weight_at_age_m[data.year_age_index(y, a)] *
                           data.maturity_at_age_m[static_cast<size_t>(a)] *
                           std::exp(-data.natural_mortality_at_age_m[
                                        static_cast<size_t>(a)] *
                                    data.spawning_fraction_m);
            return out;
        }

        inline double beverton_holt(double ssb, const Dynamics &d, double b0)
        {
            return 4.0 * d.steepness * d.r0 * ssb /
                   (b0 * (1.0 - d.steepness) +
                        ssb * (5.0 * d.steepness - 1.0) + 1e-12);
        }

        inline std::vector<double> unfished_state(
            const TunaSpatialAssessmentData &data, double r0)
        {
            std::vector<double> state(static_cast<size_t>(
                data.n_regions_m * data.n_ages_m), 0.0);
            for (int r = 0; r < data.n_regions_m; ++r)
            {
                const double recruit = r0 * data.regional_recruit_proportions_m[
                                                static_cast<size_t>(r)];
                state[data.region_age_index(r, 0)] = recruit;
                double survival = 1.0;
                for (int a = 1; a < data.n_ages_m - 1; ++a)
                {
                    survival *= std::exp(-data.natural_mortality_at_age_m[
                        static_cast<size_t>(a - 1)]);
                    state[data.region_age_index(r, a)] = recruit * survival;
                }
                const int plus = data.n_ages_m - 1;
                survival *= std::exp(-data.natural_mortality_at_age_m[
                    static_cast<size_t>(plus - 1)]);
                const double plus_survival = std::exp(
                    -data.natural_mortality_at_age_m[static_cast<size_t>(plus)]);
                state[data.region_age_index(r, plus)] =
                    recruit * survival / (1.0 - plus_survival + 1e-12);
            }
            return state;
        }

        struct YearResult
        {
            std::vector<double> post_season_state;
            double retained = 0.0;
            double discard = 0.0;
        };

        inline YearResult fish_one_year(const TunaSpatialAssessmentData &data,
                                        const Dynamics &d,
                                        const TunaFishingPattern &pattern,
                                        const std::vector<double> &initial,
                                        double multiplier)
        {
            YearResult result;
            result.post_season_state = initial;
            const int weight_year = data.n_years_m - 1;
            for (int s = 0; s < data.n_seasons_m; ++s)
            {
                std::vector<double> survivors(result.post_season_state.size(), 0.0);
                std::vector<double> moved(result.post_season_state.size(), 0.0);
                for (int r = 0; r < data.n_regions_m; ++r)
                    for (int a = 0; a < data.n_ages_m; ++a)
                    {
                        const size_t ra = data.region_age_index(r, a);
                        const double n = result.post_season_state[ra];
                        const double m = data.natural_mortality_at_age_m[
                                             static_cast<size_t>(a)] /
                                         data.n_seasons_m;
                        double total_f = 0.0;
                        for (int f = 0; f < data.n_fleets_m; ++f)
                            total_f += multiplier * pattern.rate_m[
                                fsra(data, f, s, r, a)];
                        const double weight = data.weight_at_age_m[
                            data.year_age_index(weight_year, a)];
                        for (int f = 0; f < data.n_fleets_m; ++f)
                        {
                            const size_t fa = static_cast<size_t>(
                                f * data.n_ages_m + a);
                            const double fleet_rate = multiplier * pattern.rate_m[
                                fsra(data, f, s, r, a)];
                            const double z = m + total_f + 1e-12;
                            const double exploitation =
                                (1.0 - std::exp(-z)) / z;
                            const double capture = pattern.removal_fraction_m
                                ? n * fleet_rate
                                : n * fleet_rate * exploitation;
                            result.retained += capture * d.retention[fa] * weight;
                            result.discard += capture * (1.0 - d.retention[fa]) * weight;
                        }
                        if (pattern.removal_fraction_m)
                        {
                            const double remaining = 1.0 - total_f;
                            const double positive_remaining =
                                std::log1p(std::exp(40.0 * remaining)) / 40.0;
                            survivors[ra] = n * positive_remaining * std::exp(-m);
                        }
                        else
                        {
                            survivors[ra] = n * std::exp(-(m + total_f));
                        }
                    }
                for (int from = 0; from < data.n_regions_m; ++from)
                    for (int to = 0; to < data.n_regions_m; ++to)
                        for (int a = 0; a < data.n_ages_m; ++a)
                            moved[data.region_age_index(to, a)] +=
                                survivors[data.region_age_index(from, a)] *
                                d.movement[srr(data, s, from, to)];
                result.post_season_state.swap(moved);
            }
            return result;
        }

        inline std::vector<double> recruit_and_age(
            const TunaSpatialAssessmentData &data, const Dynamics &d,
            const std::vector<double> &post_season, double b0,
            double recruitment_multiplier = 1.0)
        {
            const double recruitment =
                beverton_holt(spawning_biomass(data, post_season), d, b0) *
                recruitment_multiplier;
            std::vector<double> next(post_season.size(), 0.0);
            for (int r = 0; r < data.n_regions_m; ++r)
            {
                next[data.region_age_index(r, 0)] = recruitment *
                    data.regional_recruit_proportions_m[static_cast<size_t>(r)];
                for (int a = 0; a < data.n_ages_m - 1; ++a)
                    next[data.region_age_index(r, a + 1)] +=
                        post_season[data.region_age_index(r, a)];
                const int plus = data.n_ages_m - 1;
                next[data.region_age_index(r, plus)] +=
                    post_season[data.region_age_index(r, plus)];
            }
            return next;
        }

        struct Equilibrium
        {
            bool converged = false;
            double ssb = 0.0;
            double yield = 0.0;
        };

        inline Equilibrium equilibrium(const TunaSpatialAssessmentData &data,
                                       const Dynamics &d,
                                       const TunaFishingPattern &pattern,
                                       double multiplier, double b0)
        {
            std::vector<double> state = unfished_state(data, d.r0);
            Equilibrium out;
            for (int iteration = 0; iteration < 1000; ++iteration)
            {
                const YearResult year =
                    fish_one_year(data, d, pattern, state, multiplier);
                std::vector<double> next =
                    recruit_and_age(data, d, year.post_season_state, b0);
                double max_relative = 0.0;
                for (size_t i = 0; i < state.size(); ++i)
                    max_relative = std::max(
                        max_relative, std::abs(next[i] - state[i]) /
                                          std::max(1.0, std::abs(state[i])));
                state.swap(next);
                if (iteration >= 20 && max_relative < 1e-10)
                {
                    out.converged = true;
                    const YearResult final_year =
                        fish_one_year(data, d, pattern, state, multiplier);
                    out.ssb = spawning_biomass(data, state);
                    out.yield = final_year.retained + final_year.discard;
                    break;
                }
            }
            return out;
        }
    } // namespace tuna_reference_detail

    inline TunaFishingPattern make_tuna_fishing_pattern(
        const TunaSpatialAssessmentData &data,
        const TunaAssessmentControls &controls,
        const std::vector<double> &parameters,
        int recent_years = 3)
    {
        data.validate();
        if (recent_years <= 0)
            throw std::invalid_argument("recent_years must be positive");
        const auto dynamics = tuna_reference_detail::unpack_dynamics(
            data, controls, parameters);
        TunaFishingPattern pattern;
        pattern.rate_m.assign(static_cast<size_t>(
            data.n_fleets_m * data.n_seasons_m * data.n_regions_m *
            data.n_ages_m), 0.0);

        if (controls.use_catch_conditioning_m)
        {
            ModelReportContext context;
            (void)tuna_reference_detail::terminal_state(
                data, controls, parameters, &context);
            const int y = data.n_years_m - 1;
            for (int f = 0; f < data.n_fleets_m; ++f)
                for (int s = 0; s < data.n_seasons_m; ++s)
                    for (int r = 0; r < data.n_regions_m; ++r)
                    {
                        const double scale = tuna_reference_detail::report_value(
                            context, "diag_conditioned_scale_y" +
                                         std::to_string(y + 1) + "_s" +
                                         std::to_string(s + 1) + "_f" +
                                         std::to_string(f + 1) + "_r" +
                                         std::to_string(r + 1));
                        for (int a = 0; a < data.n_ages_m; ++a)
                        {
                            const size_t fa = static_cast<size_t>(
                                f * data.n_ages_m + a);
                            pattern.rate_m[tuna_reference_detail::fsra(
                                data, f, s, r, a)] =
                                scale * dynamics.selectivity[fa] *
                                dynamics.availability[tuna_reference_detail::fsra(
                                    data, f, s, r, a)];
                        }
                    }
            pattern.source_m = "terminal fitted catch-conditioned F-at-age";
            pattern.removal_fraction_m = true;
        }
        else
        {
            const int n = std::min(recent_years, data.n_years_m);
            const int first = data.n_years_m - n;
            size_t pos = 2;
            std::vector<double> q(static_cast<size_t>(data.n_fleets_m));
            for (int f = 0; f < data.n_fleets_m; ++f)
            {
                q[static_cast<size_t>(f)] = std::exp(parameters[pos]);
                pos += 10 + (controls.estimate_availability_scales_m
                    ? (controls.availability_by_fleet_only_m
                           ? 1
                           : data.n_seasons_m * data.n_regions_m)
                    : 0);
            }
            for (int f = 0; f < data.n_fleets_m; ++f)
                for (int s = 0; s < data.n_seasons_m; ++s)
                    for (int r = 0; r < data.n_regions_m; ++r)
                    {
                        double mean_effort = 0.0;
                        for (int y = first; y < data.n_years_m; ++y)
                            mean_effort += data.effort_m[
                                data.fleet_year_season_region_index(f, y, s, r)];
                        mean_effort /= n;
                        for (int a = 0; a < data.n_ages_m; ++a)
                        {
                            const size_t fa = static_cast<size_t>(
                                f * data.n_ages_m + a);
                            pattern.rate_m[tuna_reference_detail::fsra(
                                data, f, s, r, a)] =
                                q[static_cast<size_t>(f)] * mean_effort *
                                dynamics.selectivity[fa] *
                                dynamics.availability[tuna_reference_detail::fsra(
                                    data, f, s, r, a)];
                        }
                    }
            pattern.source_m = "recent mean effort fitted F-at-age";
            pattern.removal_fraction_m = false;
        }
        for (double value : pattern.rate_m)
            if (!std::isfinite(value) || value < 0.0)
                throw std::runtime_error("fishing pattern contains invalid F-at-age");
        return pattern;
    }

    inline TunaReferencePoints calculate_tuna_reference_points(
        const TunaSpatialAssessmentData &data,
        const TunaAssessmentControls &controls,
        const std::vector<double> &parameters,
        double terminal_ssb,
        double maximum_multiplier = 3.0,
        int grid_intervals = 150)
    {
        TunaReferencePoints out;
        if (!(maximum_multiplier > 0.0) || grid_intervals < 10)
        {
            out.message_m = "invalid reference-point grid settings";
            return out;
        }
        try
        {
            const auto dynamics = tuna_reference_detail::unpack_dynamics(
                data, controls, parameters);
            const TunaFishingPattern pattern = make_tuna_fishing_pattern(
                data, controls, parameters);
            out.fishing_pattern_source_m = pattern.source_m;
            const auto unfished = tuna_reference_detail::equilibrium(
                data, dynamics, pattern, 0.0,
                tuna_reference_detail::spawning_biomass(
                    data, tuna_reference_detail::unfished_state(data, dynamics.r0)));
            if (!unfished.converged || !(unfished.ssb > 0.0))
                throw std::runtime_error("unfished equilibrium failed");
            out.b0_m = unfished.ssb;
            double best_yield = -1.0;
            int best_index = -1;
            double search_maximum = maximum_multiplier;
            for (int expansion = 0; expansion < 5; ++expansion)
            {
                best_yield = -1.0;
                best_index = -1;
                for (int i = 0; i <= grid_intervals; ++i)
                {
                    const double multiplier =
                        search_maximum * i / grid_intervals;
                    const auto equilibrium = tuna_reference_detail::equilibrium(
                        data, dynamics, pattern, multiplier, out.b0_m);
                    if (!equilibrium.converged)
                        continue;
                    if (equilibrium.yield > best_yield)
                    {
                        best_yield = equilibrium.yield;
                        best_index = i;
                        out.msy_m = equilibrium.yield;
                        out.b_msy_m = equilibrium.ssb;
                        out.f_msy_multiplier_m = multiplier;
                    }
                }
                if (best_index >= 0 && best_index < grid_intervals)
                    break;
                // Effort-driven rates are instantaneous F and can safely be
                // searched more broadly. Catch-conditioned rates are removal
                // fractions, so extrapolation beyond the fitted range is not
                // biologically valid.
                if (pattern.removal_fraction_m)
                    break;
                search_maximum *= 2.0;
            }
            if (best_index < 0 || !(out.b_msy_m > 0.0) ||
                !(out.f_msy_multiplier_m > 0.0))
                throw std::runtime_error("no finite MSY equilibrium found");
            out.grid_boundary_m = best_index == grid_intervals;
            out.terminal_b_over_b_msy_m = terminal_ssb / out.b_msy_m;
            out.status_quo_f_over_f_msy_m = 1.0 / out.f_msy_multiplier_m;
            out.valid_m = !out.grid_boundary_m;
            out.message_m = out.grid_boundary_m
                ? "MSY occurred at the upper search boundary"
                : "equilibrium reference points calculated";
        }
        catch (const std::exception &error)
        {
            out.message_m = error.what();
        }
        return out;
    }

    inline std::vector<TunaProjectionPoint> project_tuna_scenario(
        const TunaSpatialAssessmentData &data,
        const TunaAssessmentControls &controls,
        const std::vector<double> &parameters,
        double fishing_multiplier,
        int years,
        const std::vector<double> &recruitment_multipliers = {})
    {
        if (!(fishing_multiplier >= 0.0) || !std::isfinite(fishing_multiplier) ||
            years <= 0)
            throw std::invalid_argument("invalid tuna projection scenario");
        if (!recruitment_multipliers.empty() &&
            recruitment_multipliers.size() != static_cast<size_t>(years))
            throw std::invalid_argument(
                "projection recruitment multipliers must match years");
        const auto dynamics = tuna_reference_detail::unpack_dynamics(
            data, controls, parameters);
        const TunaFishingPattern pattern = make_tuna_fishing_pattern(
            data, controls, parameters);
        std::vector<double> state = tuna_reference_detail::terminal_state(
            data, controls, parameters);
        const double b0 = tuna_reference_detail::equilibrium(
            data, dynamics, pattern, 0.0,
            tuna_reference_detail::spawning_biomass(
                data, tuna_reference_detail::unfished_state(data, dynamics.r0))).ssb;
        // The fitted terminal state is after the last fishing season. Advance
        // recruitment and ageing before reporting projection year 1.
        state = tuna_reference_detail::recruit_and_age(
            data, dynamics, state, b0,
            recruitment_multipliers.empty() ? 1.0
                                            : recruitment_multipliers.front());
        std::vector<TunaProjectionPoint> out;
        out.reserve(static_cast<size_t>(years));
        for (int year = 1; year <= years; ++year)
        {
            TunaProjectionPoint point;
            point.projection_year_m = year;
            point.fishing_multiplier_m = fishing_multiplier;
            point.spawning_biomass_m =
                tuna_reference_detail::spawning_biomass(data, state);
            point.depletion_m = point.spawning_biomass_m / b0;
            const auto fished = tuna_reference_detail::fish_one_year(
                data, dynamics, pattern, state, fishing_multiplier);
            point.retained_yield_m = fished.retained;
            point.discard_yield_m = fished.discard;
            point.total_yield_m = fished.retained + fished.discard;
            out.push_back(point);
            if (year < years)
                state = tuna_reference_detail::recruit_and_age(
                    data, dynamics, fished.post_season_state, b0,
                    recruitment_multipliers.empty()
                        ? 1.0
                        : recruitment_multipliers[static_cast<size_t>(year)]);
        }
        return out;
    }

    inline std::string tuna_reference_points_csv(const TunaReferencePoints &x)
    {
        const auto number = [](double value)
        {
            if (!std::isfinite(value))
                return std::string();
            std::ostringstream text;
            text.precision(17);
            text << value;
            return text.str();
        };
        std::ostringstream out;
        out << "metric,value,valid,description\n"
            << "B0," << number(x.b0_m) << "," << x.valid_m
            << ",equilibrium unfished spawning biomass\n"
            << "B_MSY," << number(x.b_msy_m) << "," << x.valid_m
            << ",equilibrium spawning biomass at MSY\n"
            << "MSY," << number(x.msy_m) << "," << x.valid_m
            << ",equilibrium total catch biomass at MSY\n"
            << "F_MSY_multiplier," << number(x.f_msy_multiplier_m) << "," << x.valid_m
            << ",multiplier on the stated fishing pattern\n"
            << "B_terminal_over_B_MSY," << number(x.terminal_b_over_b_msy_m) << ","
            << x.valid_m << ",terminal spawning biomass status\n"
            << "F_status_quo_over_F_MSY," << number(x.status_quo_f_over_f_msy_m) << ","
            << x.valid_m << ",status-quo fishing pressure status\n";
        return out.str();
    }

    inline std::string tuna_projection_csv(
        const std::string &scenario,
        const std::vector<TunaProjectionPoint> &points)
    {
        std::ostringstream out;
        out << "scenario,projection_year,fishing_multiplier,spawning_biomass,"
               "depletion,retained_yield,discard_yield,total_yield\n";
        for (const auto &point : points)
            out << scenario << "," << point.projection_year_m << ","
                << point.fishing_multiplier_m << "," << point.spawning_biomass_m
                << "," << point.depletion_m << "," << point.retained_yield_m
                << "," << point.discard_yield_m << "," << point.total_yield_m
                << "\n";
        return out.str();
    }
} // namespace quadra

#endif
