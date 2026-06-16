#include "../data/pollock_data.hpp"
#include "../data/pollock_io.hpp"
#include "reports/pollock_reports.hpp"
#include "reports/pollock_fit_summary.hpp"
#include "drivers/pollock_driver_output.hpp"
#include "diagnostics/pollock_utilities.hpp"
#include "diagnostics/pollock_fixed_effect_diagnostics.hpp"
#include "diagnostics/pollock_huu_diagnostics.hpp"
#include "diagnostics/pollock_huu_output_diagnostics.hpp"
#include "diagnostics/pollock_fixed_hessian_diagnostics.hpp"
#include "diagnostics/pollock_functional_analysis_diagnostics.hpp"
#include "model/pollock_parameters.hpp"
#include "model/pollock_constants.hpp"
#include "model/pollock_model.hpp"
#include "../../../../core/optimizer.hpp"

#include <iostream>
#include <stdexcept>

int main()
{
  try
  {
    std::cout << "Synthetic AFSC walleye-pollock-style assessment example\n";
    std::cout << "=======================================================\n\n";
    std::cout << "Synthetic and public-data-safe. Not an official assessment.\n";
    std::cout << "Assessment-scale diagnostic: tolerance is relaxed for synthetic profiling/identifiability checks.\n";
    std::cout << "Recruitment deviations use a fixed AR(1) prior: rho=0.60, sigma=0.15.\n";
#ifdef WALLEYE_POLLOCK_RANDOM_RECRUITMENT_COUNT
    std::cout << "Random recruitment enabled for first "
              << WALLEYE_POLLOCK_RANDOM_RECRUITMENT_COUNT
              << " year(s).\n\n";
#else
    std::cout << "Level 1: fixed-effect index fit with observed-catch removals; random recruitment disabled.\n\n";
#endif

    auto obs = read_obs("examples/NMFS/afsc_walleye_pollock/data/synthetic_walleye_pollock_observations.csv");
    std::cout << "Loaded synthetic rows: " << obs.size() << "\n\n";

    PollockModel model(obs);
    auto params = pollock::make_params(obs.size());
    auto opts = quadra::default_laplace_options();

    auto fit = quadra::optimize_lbfgs(model, params, opts);

    write_summary("examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_fit_summary.csv", fit);
    write_fixed_parameter_estimates(
        "examples/NMFS/afsc_walleye_pollock/outputs/"
        "walleye_pollock_fixed_parameter_estimates.csv",
        fit);
    write_fixed_gradient_diagnostics(
        "examples/NMFS/afsc_walleye_pollock/outputs/"
        "walleye_pollock_fixed_gradient_diagnostics.csv",
        fit);

#ifdef WALLEYE_POLLOCK_FIXED_HESSIAN_DIAGNOSTICS
    {
      quadra::LaplaceOptions hess_opts = quadra::default_laplace_options();
      write_fixed_hessian_diagnostics(
          "examples/NMFS/afsc_walleye_pollock/outputs/"
          "walleye_pollock_fixed_hessian_diagnostics.csv",
          "examples/NMFS/afsc_walleye_pollock/outputs/"
          "walleye_pollock_fixed_hessian_matrix.csv",
          model, params, fit, hess_opts);
    }
#endif

#ifdef WALLEYE_POLLOCK_HUU_DIAGNOSTICS
    pollock_write_huu_diagnostics(
        "examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_huu_diagnostics.csv",
        model, params, fit);
#endif

#ifdef WALLEYE_POLLOCK_HUU_MATRIX_DUMP
    pollock_write_huu_matrix(
        "examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_huu_matrix.csv",
        model, params, fit);
    pollock_write_huu_sparsity(
        "examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_huu_sparsity.csv",
        model, params, fit);
#endif

#ifdef WALLEYE_POLLOCK_HUU_PATTERN_COMPARE
    pollock_write_huu_pattern_compare(
        "examples/NMFS/afsc_walleye_pollock/outputs/"
        "walleye_pollock_huu_pattern_compare.csv",
        model, params, fit);
#endif

#ifdef WALLEYE_POLLOCK_HUU_BAND_SUMMARY
    pollock_write_huu_band_summary(
        "examples/NMFS/afsc_walleye_pollock/outputs/"
        "walleye_pollock_huu_band_summary.csv",
        model, params, fit);
#endif

#ifdef WALLEYE_POLLOCK_HUU_BANDLIMIT_DIAGNOSTIC
    pollock_write_huu_bandlimit_diagnostic(
        "examples/NMFS/afsc_walleye_pollock/outputs/"
        "walleye_pollock_huu_bandlimit_diagnostic.csv",
        model, params, fit);
#endif

#ifdef WALLEYE_POLLOCK_HUU_THRESHOLD_DIAGNOSTIC
    pollock_write_huu_threshold_diagnostic(
        "examples/NMFS/afsc_walleye_pollock/outputs/"
        "walleye_pollock_huu_threshold_diagnostic.csv",
        model, params, fit);
#endif

#ifdef WALLEYE_POLLOCK_LAPLACE_STRUCTURE_REPORT
    pollock_write_laplace_structure_report(
        "examples/NMFS/afsc_walleye_pollock/outputs/"
        "walleye_pollock_laplace_structure_report.txt",
        model, params, fit);
#endif

#ifdef WALLEYE_POLLOCK_FUNCTIONAL_ANALYSIS_REPORT
    pollock_write_functional_analysis_report(
        "examples/NMFS/afsc_walleye_pollock/outputs/"
        "walleye_pollock_functional_analysis_report.txt",
        "examples/NMFS/afsc_walleye_pollock/outputs/"
        "walleye_pollock_functional_analysis_report.csv",
        model, params, fit);
#endif

#ifdef WALLEYE_POLLOCK_MARKDOWN_REPORT
    pollock_example::write_pollock_markdown_report(
        "examples/NMFS/afsc_walleye_pollock/outputs/"
        "walleye_pollock_analysis.md",
        "examples/NMFS/afsc_walleye_pollock/outputs/"
        "walleye_pollock_functional_analysis_report.csv",
        "examples/NMFS/afsc_walleye_pollock/outputs/"
        "walleye_pollock_laplace_structure_report.txt");
#endif

    pollock_example::write_recruitment_deviations(
        "examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_recruitment_deviations.csv",
        fit);

    pollock_example::print_fit_and_structure_diagnostics(fit);
    pollock_example::print_output_manifest();

    return fit.converged ? 0 : 2;
  }
  catch (const std::exception &e)
  {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 1;
  }
}
