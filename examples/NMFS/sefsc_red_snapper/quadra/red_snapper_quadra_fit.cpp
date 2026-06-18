#include "../diagnostics/red_snapper_functional_analysis_diagnostics.hpp"
#include "../objective/red_snapper_quadra_objective.hpp"
#include "../reports/red_snapper_report_suite.hpp"
#include "red_snapper_age_structured.hpp"

#include "../../../../core/optimizer.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int main()
{
    const std::string input_path = "examples/NMFS/sefsc_red_snapper/data/"
                                   "synthetic_red_snapper_observations.csv";
    const auto report_paths =
        sefsc_red_snapper::default_red_snapper_report_paths();
    const auto observations = sefsc_red_snapper::read_observations(input_path);

    sefsc_red_snapper::RedSnapperQuadraObjective objective(observations);

    quadra::ParameterVector params;
    params.add({"log_r0", std::log(1200.0), quadra::ParameterTransform::Identity,
                false});
    params.add({"log_fbar", std::log(0.025), quadra::ParameterTransform::Identity,
                false});
    params.add({"log_q", std::log(0.00005), quadra::ParameterTransform::Identity,
                false});
    params.add(
        {"logit_sel_a50", 0.0, quadra::ParameterTransform::Identity, false});
    params.add({"log_sel_slope", std::log(1.2),
                quadra::ParameterTransform::Identity, false});

    for (std::size_t t = 0; t < observations.size(); ++t)
    {
        params.add({"log_rec_dev_" + std::to_string(t + 1), 0.0,
                    quadra::ParameterTransform::Identity, true});
    }

    quadra::LaplaceOptions opts;

    auto fit = quadra::optimize_lbfgs(objective, params, opts);

    sefsc_red_snapper::write_red_snapper_report_suite(report_paths, observations,
                                                      objective, params, fit);
    sefsc_red_snapper::write_red_snapper_functional_analysis_report(
        "examples/NMFS/sefsc_red_snapper/outputs/"
        "red_snapper_functional_analysis_report.txt",
        "examples/NMFS/sefsc_red_snapper/outputs/"
        "red_snapper_functional_analysis_report.csv",
        objective, params, fit);
    std::cout
        << "SEFSC red-snapper-style Quadra Laplace recruitment-deviation fit\n";
    std::cout << "objective:  " << fit.value << "\n";
    std::cout << "grad_norm:  " << fit.grad_norm << "\n";
    std::cout << "converged:  " << (fit.converged ? "yes" : "no") << "\n";
    std::cout << "message:    " << fit.message << "\n";
    std::cout << "wrote:      " << report_paths.summary << "\n";
    std::cout << "wrote:      " << report_paths.trajectory << "\n";
    std::cout << "wrote:      " << report_paths.residual_diagnostics << "\n";
    std::cout << "wrote:      " << report_paths.selectivity << "\n";
    std::cout << "wrote:      " << report_paths.recruitment_deviations << "\n";
    std::cout << "wrote:      " << report_paths.objective_components << "\n";
    std::cout << "wrote:      " << report_paths.laplace_structure_text << "\n";
    std::cout << "wrote:      " << report_paths.laplace_structure_csv << "\n";
    std::cout << "wrote:      "
              << "examples/NMFS/sefsc_red_snapper/outputs/"
                 "red_snapper_functional_analysis_report.txt\n";

    std::cout << "wrote:      "
              << "examples/NMFS/sefsc_red_snapper/outputs/"
                 "red_snapper_functional_analysis_report.csv\n";
    std::cout << "wrote:      " << report_paths.reference_points << "\n";
    return 0;
}
