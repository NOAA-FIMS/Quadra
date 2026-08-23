#include <cmath>
#include <iostream>
#include <vector>

#include "tuna/tuna_reference_points.hpp"
#include "tuna/tuna_spatial_assessment_model.hpp"

int main()
{
    quadra::TunaSpatialAssessmentData data;
    data.n_years_m = 2;
    data.n_ages_m = 2;
    data.n_fleets_m = 1;
    data.n_regions_m = 1;
    data.n_seasons_m = 1;
    data.natural_mortality_at_age_m = {0.2, 0.2};
    data.maturity_at_age_m = {0.0, 1.0};
    data.weight_at_age_m = {1.0, 2.0, 1.0, 2.0};
    data.movement_matrix_m = {1.0};
    data.regional_recruit_proportions_m = {1.0};
    data.effort_m.assign(2, 0.0);
    data.observed_index_m.assign(2, 0.0);
    data.observed_retained_biomass_m.assign(2, 0.0);
    data.observed_discard_biomass_m.assign(2, 0.0);
    data.observed_catch_numbers_m.assign(4, 0);
    data.observed_total_catch_m = {10.0, 12.0};
    data.catch_units_m = {1};

    quadra::TunaAssessmentControls controls;
    controls.use_catch_conditioning_m = true;
    controls.use_index_likelihood_m = false;
    controls.use_catch_composition_likelihood_m = false;
    controls.use_retained_biomass_likelihood_m = false;
    controls.use_discard_biomass_likelihood_m = false;
    controls.use_priors_m = false;
    controls.estimate_availability_scales_m = false;

    quadra::AdvancedSpatialTunaAssessmentModel model(data, controls);
    std::vector<double> parameters = model.parameter_set().initials();
    quadra::ModelReportContext first_context;
    const double first = model.evaluate(parameters, first_context);
    parameters[2] = 2.0; // log_q must be irrelevant in catch-conditioned mode.
    quadra::ModelReportContext second_context;
    const double second = model.evaluate(parameters, second_context);

    if (!std::isfinite(first) || !std::isfinite(second) ||
        std::abs(first - second) > 1e-10)
    {
        std::cerr << "catch-conditioning smoke test failed: "
                  << first << " versus " << second << "\n";
        return 1;
    }

    quadra::TunaAssessmentControls effort_controls = controls;
    effort_controls.use_catch_conditioning_m = false;
    data.observed_total_catch_m.clear();
    data.catch_units_m.clear();
    data.effort_m.assign(2, 1.0);
    quadra::AdvancedSpatialTunaAssessmentModel effort_model(data, effort_controls);
    std::vector<double> effort_parameters =
        effort_model.parameter_set().initials();
    effort_parameters[2] = std::log(0.2);
    const quadra::TunaReferencePoints reference =
        quadra::calculate_tuna_reference_points(
            data, effort_controls, effort_parameters, 100.0, 10.0, 100);
    if (!reference.valid_m || !std::isfinite(reference.b0_m) ||
        !std::isfinite(reference.msy_m) || !(reference.msy_m > 0.0))
    {
        std::cerr << "reference-point smoke test failed: "
                  << reference.message_m << "\n";
        return 1;
    }
    const auto projection = quadra::project_tuna_scenario(
        data, effort_controls, effort_parameters,
        reference.f_msy_multiplier_m, 3);
    const auto low_recruitment_projection = quadra::project_tuna_scenario(
        data, effort_controls, effort_parameters,
        reference.f_msy_multiplier_m, 3, {0.5, 0.5, 0.5});
    if (projection.size() != 3 ||
        !std::isfinite(projection.back().spawning_biomass_m) ||
        !std::isfinite(projection.back().total_yield_m) ||
        !(low_recruitment_projection.back().spawning_biomass_m <
          projection.back().spawning_biomass_m))
    {
        std::cerr << "projection smoke test failed\n";
        return 1;
    }
    std::cout << "catch-conditioning smoke test: PASS\n";
    std::cout << "reference-point and projection smoke test: PASS\n";
    return 0;
}
