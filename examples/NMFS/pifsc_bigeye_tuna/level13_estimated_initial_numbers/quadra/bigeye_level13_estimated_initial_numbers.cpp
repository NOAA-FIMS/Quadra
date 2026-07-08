#include "../diagnostics/bigeye_functional_analysis_diagnostics.hpp"
#include "../diagnostics/bigeye_safe_fixed_effect_diagnostics.hpp"
#include "../diagnostics/bigeye_longline_slope_geometry_scan.hpp"
#include "../diagnostics/bigeye_recruitment_diagnostics.hpp"
#include "../diagnostics/bigeye_initial_numbers_diagnostics.hpp"
#include "../diagnostics/bigeye_fixed_effect_geometry.hpp"
#include "../diagnostics/bigeye_fixed_effect_wiggle_diagnostics.hpp"
#include "../objective/bigeye_quadra_objective.hpp"
#include "../reports/bigeye_report_suite.hpp"
#include "bigeye_age_structured.hpp"

#include "../../../../../core/optimizer.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::vector<pifsc_bigeye_tuna::FleetObservation>
read_multifleet_observations(const std::string &path)
{
    std::ifstream in(path);
    if (!in)
    {
        throw std::runtime_error("Could not open multifleet observations CSV: " +
                                 path);
    }

    std::string line;
    std::getline(in, line);

    std::vector<pifsc_bigeye_tuna::FleetObservation> out;

    while (std::getline(in, line))
    {
        if (line.empty())
            continue;

        std::stringstream ss(line);
        std::string year_s;
        std::string fleet;
        std::string catch_s;
        std::string index_s;

        std::getline(ss, year_s, ',');
        std::getline(ss, fleet, ',');
        std::getline(ss, catch_s, ',');
        std::getline(ss, index_s, ',');

        pifsc_bigeye_tuna::FleetObservation obs;
        obs.year = std::stoi(year_s);
        obs.fleet = fleet;
        obs.catch_mt = std::stod(catch_s);
        obs.index = std::stod(index_s);

        double comp_sum = 0.0;
        for (int a = 0; a < pifsc_bigeye_tuna::kAges; ++a)
        {
            std::string comp_s;
            if (std::getline(ss, comp_s, ',') && !comp_s.empty())
            {
                obs.age_comp[static_cast<std::size_t>(a)] = std::stod(comp_s);
                comp_sum += obs.age_comp[static_cast<std::size_t>(a)];
            }
        }

        if (comp_sum > 0.0)
        {
            for (double &p_age : obs.age_comp)
                p_age /= comp_sum;
        }
        else
        {
            const double inv_n_age =
                1.0 / static_cast<double>(pifsc_bigeye_tuna::kAges);
            for (double &p_age : obs.age_comp)
                p_age = inv_n_age;
        }

        out.push_back(obs);
    }

    return out;
}

std::vector<pifsc_bigeye_tuna::Observation>
aggregate_fleet_observations_for_reports(
    const std::vector<pifsc_bigeye_tuna::FleetObservation> &fleet_obs)
{
    struct Accumulator
    {
        double catch_mt = 0.0;
        double index_sum = 0.0;
        int index_count = 0;
    };

    std::map<int, Accumulator> by_year;

    for (const auto &obs : fleet_obs)
    {
        auto &acc = by_year[obs.year];
        acc.catch_mt += obs.catch_mt;
        acc.index_sum += obs.index;
        acc.index_count += 1;
    }

    std::vector<pifsc_bigeye_tuna::Observation> out;
    for (const auto &kv : by_year)
    {
        pifsc_bigeye_tuna::Observation obs;
        obs.year = kv.first;
        obs.catch_mt = kv.second.catch_mt;
        obs.index =
            kv.second.index_count > 0
                ? kv.second.index_sum / static_cast<double>(kv.second.index_count)
                : 0.0;

        const double inv_n_age =
            1.0 / static_cast<double>(pifsc_bigeye_tuna::kAges);
        for (double &p_age : obs.age_comp)
            p_age = inv_n_age;

        out.push_back(obs);
    }

    return out;
}

} // namespace

int main()
{
    const std::string input_path =
        "examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/data/"
        "synthetic_bigeye_level13_estimated_initial_numbers_observations.csv";
    const auto report_paths =
        pifsc_bigeye_tuna::default_bigeye_report_paths();
    const auto fleet_observations = read_multifleet_observations(input_path);
    const auto observations = aggregate_fleet_observations_for_reports(fleet_observations);

    pifsc_bigeye_tuna::BigeyeQuadraObjective objective(fleet_observations);

    quadra::ParameterVector params;
    params.add({"log_r0", std::log(1200.0), quadra::ParameterTransform::Identity,
                false});
    params.add({"log_fbar", std::log(0.025), quadra::ParameterTransform::Identity,
                false});
    params.add({"log_q_purse_seine", std::log(0.00005),
                quadra::ParameterTransform::Identity, false});
    params.add({"logit_sel_a50_longline", 0.0,
                quadra::ParameterTransform::Identity, false});
    params.add({"log_sel_slope_longline", std::log(1.2),
                quadra::ParameterTransform::Identity, false});


    for (int a = 0; a < pifsc_bigeye_tuna::kAges; ++a) {
      params.add({"init_log_number_dev_age_" + std::to_string(a + 1), 0.0,
                  quadra::ParameterTransform::Identity, false});
    }

    for (std::size_t t = 0; t < objective.n_years(); ++t)
    {
        params.add({"log_rec_dev_" + std::to_string(t + 1), 0.0,
                    quadra::ParameterTransform::Identity, true});
    }

    quadra::LaplaceOptions opts;

    auto fit = quadra::optimize_lbfgs(objective, params, opts);

    pifsc_bigeye_tuna::write_bigeye_report_suite(report_paths, observations,
                                                      objective, params, fit);

    {
      std::vector<double> full_parameter_vector;
      full_parameter_vector.reserve(fit.par.size() + fit.u_hat.size());
      for (const auto v : fit.par) {
        full_parameter_vector.push_back(v);
      }
      for (const auto v : fit.u_hat) {
        full_parameter_vector.push_back(v);
      }

      const double direct_joint_objective = objective(full_parameter_vector);
      const double direct_minus_fit = direct_joint_objective - fit.value;

      std::ofstream out(
          "examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/outputs/"
          "bigeye_level13_objective_consistency_check.csv");
      if (!out) {
        throw std::runtime_error("Could not open objective consistency check CSV");
      }

      out << std::setprecision(15);
      out << "metric,value,note\n";
      out << "fit_value," << fit.value
          << ",reported profiled Laplace objective from optimizer\n";
      out << "direct_joint_objective," << direct_joint_objective
          << ",objective evaluated directly at fit.par plus fit.u_hat\n";
      out << "direct_minus_fit_value," << direct_minus_fit
          << ",direct joint objective minus fit.value; expected to differ because fit.value includes Laplace terms\n";
      out << "fit_par_size," << fit.par.size()
          << ",number of fixed effects in returned fit\n";
      out << "fit_u_hat_size," << fit.u_hat.size()
          << ",number of random effects in returned fit\n";
      out << "full_parameter_vector_size," << full_parameter_vector.size()
          << ",fixed plus random effects passed to objective\n";

      std::ofstream txt(
          "examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/outputs/"
          "bigeye_level13_objective_consistency_check.txt");
      if (!txt) {
        throw std::runtime_error("Could not open objective consistency check text");
      }

      txt << std::setprecision(15);
      txt << "Objective Consistency Check\n";
      txt << "===========================\n\n";
      txt << "fit_value:                  " << fit.value << "\n";
      txt << "direct_joint_objective:     " << direct_joint_objective << "\n";
      txt << "direct_minus_fit_value:     " << direct_minus_fit << "\n";
      txt << "fit_par_size:               " << fit.par.size() << "\n";
      txt << "fit_u_hat_size:             " << fit.u_hat.size() << "\n";
      txt << "full_parameter_vector_size: " << full_parameter_vector.size() << "\n\n";
      txt << "Interpretation\n";
      txt << "--------------\n";
      txt << "The direct joint objective is evaluated at the exact returned fixed effects and random effects.\n";
      txt << "The optimizer fit value is the profiled Laplace objective, so it should differ by Laplace terms.\n";
      txt << "A large mismatch with the objective-components report indicates stale parameter ordering or assumptions.\n";
    }

    pifsc_bigeye_tuna::write_level13_initial_numbers_diagnostics(
        "examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/outputs/"
        "bigeye_level13_initial_numbers_diagnostics.txt",
        "examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/outputs/"
        "bigeye_level13_initial_numbers_diagnostics.csv",
        fit);
    std::cout << "wrote:      examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/outputs/bigeye_level13_initial_numbers_diagnostics.txt\n";
    std::cout << "wrote:      examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/outputs/bigeye_level13_initial_numbers_diagnostics.csv\n";

    pifsc_bigeye_tuna::write_recruitment_diagnostics(
        "examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/outputs/"
        "bigeye_level13_recruitment_diagnostics.txt",
        "examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/outputs/"
        "bigeye_level13_recruitment_diagnostics.csv",
        objective, fit);

    pifsc_bigeye_tuna::write_longline_slope_geometry_scan(
        "examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/outputs/"
        "bigeye_level13_longline_slope_geometry_scan.txt",
        "examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/outputs/"
        "bigeye_level13_longline_slope_geometry_scan.csv",
        objective, params, fit, opts);

    pifsc_bigeye_tuna::write_safe_fixed_effect_wiggle_diagnostics(
        "examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/outputs/"
        "bigeye_level13_safe_fixed_effect_wiggle_diagnostics.txt",
        "examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/outputs/"
        "bigeye_level13_safe_fixed_effect_wiggle_diagnostics.csv",
        objective, params, fit, opts);

  std::cout << "wrote:      "
            << "examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/outputs/"
               "bigeye_level13_safe_fixed_effect_wiggle_diagnostics.txt\n";
  std::cout << "wrote:      "
            << "examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/outputs/"
               "bigeye_level13_safe_fixed_effect_wiggle_diagnostics.csv\n";

    std::cout << "wrote:      examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/outputs/bigeye_level13_objective_consistency_check.txt\n";
    std::cout << "wrote:      examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/outputs/bigeye_level13_objective_consistency_check.csv\n";

    return 0;
}
