#include "../architecture/assessment/assessment_cycle.hpp"
#include "../architecture/assessment/generated_assessment_cycle.hpp"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

double max_abs_diff = 0.0;
std::string first_failure;

bool check_close(const std::string &name, double handwritten, double generated,
                 double tol = 1.0e-10) {
  const double diff = std::abs(handwritten - generated);

  if (diff > max_abs_diff) {
    max_abs_diff = diff;
  }

  if (diff > tol && first_failure.empty()) {
    std::ostringstream os;
    os << std::setprecision(17)
       << "FAIL: " << name << "\n"
       << "  handwritten = " << handwritten << "\n"
       << "  generated   = " << generated << "\n"
       << "  abs diff    = " << diff << "\n"
       << "  tolerance   = " << tol << "\n";
    first_failure = os.str();
  }

  return diff <= tol;
}

template <typename T>
bool check_array(const std::string &name,
                 const std::array<T, bigeye_v2::kAges> &handwritten,
                 const std::array<T, bigeye_v2::kAges> &generated,
                 double tol = 1.0e-10) {
  bool ok = true;

  for (int age = 0; age < bigeye_v2::kAges; ++age) {
    ok = check_close(name + "[age " + std::to_string(age + 1) + "]",
                     handwritten[age], generated[age], tol) && ok;
  }

  return ok;
}

template <typename T>
bool check_vector(const std::string &name,
                  const std::vector<T> &handwritten,
                  const std::vector<T> &generated,
                  double tol = 1.0e-10) {
  bool ok = true;

  if (handwritten.size() != generated.size()) {
    if (first_failure.empty()) {
      std::ostringstream os;
      os << "FAIL: " << name << " size mismatch\n"
         << "  handwritten size = " << handwritten.size() << "\n"
         << "  generated size   = " << generated.size() << "\n";
      first_failure = os.str();
    }
    return false;
  }

  for (std::size_t i = 0; i < handwritten.size(); ++i) {
    ok = check_close(name + "[" + std::to_string(i) + "]",
                     handwritten[i], generated[i], tol) && ok;
  }

  return ok;
}

template <typename T>
bool check_vector_array(
    const std::string &name,
    const std::vector<std::array<T, bigeye_v2::kAges>> &handwritten,
    const std::vector<std::array<T, bigeye_v2::kAges>> &generated,
    double tol = 1.0e-10) {
  bool ok = true;

  if (handwritten.size() != generated.size()) {
    if (first_failure.empty()) {
      std::ostringstream os;
      os << "FAIL: " << name << " size mismatch\n"
         << "  handwritten size = " << handwritten.size() << "\n"
         << "  generated size   = " << generated.size() << "\n";
      first_failure = os.str();
    }
    return false;
  }

  for (std::size_t y = 0; y < handwritten.size(); ++y) {
    ok = check_array(name + "[year " + std::to_string(y) + "]",
                     handwritten[y], generated[y], tol) && ok;
  }

  return ok;
}

void print_status(const std::string &name, bool ok) {
  std::cout << std::left << std::setw(28) << name
            << (ok ? "PASS" : "FAIL") << "\n";
}

} // namespace

int main() {
  using namespace bigeye_v2;

  BigeyeModelData<double> data;
  data.n_years = 3;

  AssessmentParameters<double> parameters;
  parameters.life.log_m_young_offset = std::log(0.75);
  parameters.life.log_m_old_offset = std::log(0.65);

  parameters.populations.resize(1);
  parameters.populations[0].r0 = 1000.0;

  parameters.fleets.resize(1);
  parameters.fleets[0].sel_a50 = 5.0;
  parameters.fleets[0].sel_slope = 1.0;
  parameters.fleets[0].fbar = 0.2;
  parameters.fleets[0].q_index = 0.01;
  parameters.fleets[0].catch_sigma = 0.1;
  parameters.fleets[0].index_sigma = 0.2;

  AssessmentState<double> handwritten;
  AssessmentState<double> generated;

  handwritten.populations.resize(1);
  handwritten.fleets.resize(1);
  generated.populations.resize(1);
  generated.fleets.resize(1);

  handwritten.populations[0].numbers_at_age.assign(
      static_cast<std::size_t>(data.n_years), std::array<double, kAges>{});

  generated.populations[0].numbers_at_age.assign(
      static_cast<std::size_t>(data.n_years), std::array<double, kAges>{});

  for (int a = 0; a < kAges; ++a) {
    handwritten.populations[0].numbers_at_age[0][a] = 1000.0;
    generated.populations[0].numbers_at_age[0][a] = 1000.0;
  }

  AssessmentCycle{}(data, parameters, handwritten);
  GeneratedAssessmentCycleFromIR{}(data, parameters, generated);

  std::cout << "== Generated cycle equivalence ==\n";

  bool all_ok = true;

  bool life_ok = true;
  life_ok = check_array("LifeHistoryState.m_at_age",
                        handwritten.life.m_at_age,
                        generated.life.m_at_age) && life_ok;
  life_ok = check_array("LifeHistoryState.weight_at_age",
                        handwritten.life.weight_at_age,
                        generated.life.weight_at_age) && life_ok;
  life_ok = check_array("LifeHistoryState.maturity_at_age",
                        handwritten.life.maturity_at_age,
                        generated.life.maturity_at_age) && life_ok;
  print_status("LifeHistoryState", life_ok);
  all_ok = life_ok && all_ok;

  bool population_ok = true;
  population_ok = check_vector("PopulationState.recruits_by_year",
                               handwritten.populations[0].recruits_by_year,
                               generated.populations[0].recruits_by_year) &&
                  population_ok;
  population_ok = check_vector_array("PopulationState.numbers_at_age",
                                     handwritten.populations[0].numbers_at_age,
                                     generated.populations[0].numbers_at_age) &&
                  population_ok;
  population_ok = check_vector_array("PopulationState.survivors_at_age",
                                     handwritten.populations[0].survivors_at_age,
                                     generated.populations[0].survivors_at_age) &&
                  population_ok;
  population_ok = check_vector(
                      "PopulationState.spawning_biomass_by_year",
                      handwritten.populations[0].spawning_biomass_by_year,
                      generated.populations[0].spawning_biomass_by_year) &&
                  population_ok;
  print_status("PopulationState", population_ok);
  all_ok = population_ok && all_ok;

  bool fleet_ok = true;
  fleet_ok = check_array("FleetState.selectivity_at_age",
                         handwritten.fleets[0].selectivity_at_age,
                         generated.fleets[0].selectivity_at_age) &&
             fleet_ok;
  fleet_ok = check_array("FleetState.f_at_age",
                         handwritten.fleets[0].f_at_age,
                         generated.fleets[0].f_at_age) &&
             fleet_ok;
  fleet_ok = check_array("FleetState.z_at_age",
                         handwritten.fleets[0].z_at_age,
                         generated.fleets[0].z_at_age) &&
             fleet_ok;
  fleet_ok = check_vector_array("FleetState.catch_numbers_at_age",
                                handwritten.fleets[0].catch_numbers_at_age,
                                generated.fleets[0].catch_numbers_at_age) &&
             fleet_ok;
  fleet_ok = check_vector_array("FleetState.catch_biomass_at_age",
                                handwritten.fleets[0].catch_biomass_at_age,
                                generated.fleets[0].catch_biomass_at_age) &&
             fleet_ok;
  fleet_ok = check_vector("FleetState.total_catch_biomass_by_year",
                          handwritten.fleets[0].total_catch_biomass_by_year,
                          generated.fleets[0].total_catch_biomass_by_year) &&
             fleet_ok;
  fleet_ok = check_vector("FleetState.predicted_index_by_year",
                          handwritten.fleets[0].predicted_index_by_year,
                          generated.fleets[0].predicted_index_by_year) &&
             fleet_ok;
  fleet_ok = check_vector_array(
                 "FleetState.predicted_catch_age_proportion",
                 handwritten.fleets[0].predicted_catch_age_proportion,
                 generated.fleets[0].predicted_catch_age_proportion) &&
             fleet_ok;
  print_status("FleetState", fleet_ok);
  all_ok = fleet_ok && all_ok;

  bool likelihood_ok = true;
  likelihood_ok = check_close("LikelihoodState.catch_nll",
                              handwritten.likelihood.catch_nll,
                              generated.likelihood.catch_nll) &&
                  likelihood_ok;
  likelihood_ok = check_close("LikelihoodState.index_nll",
                              handwritten.likelihood.index_nll,
                              generated.likelihood.index_nll) &&
                  likelihood_ok;
  likelihood_ok = check_close("LikelihoodState.agecomp_nll",
                              handwritten.likelihood.agecomp_nll,
                              generated.likelihood.agecomp_nll) &&
                  likelihood_ok;
  likelihood_ok = check_close("LikelihoodState.total_nll",
                              handwritten.likelihood.total_nll,
                              generated.likelihood.total_nll) &&
                  likelihood_ok;
  print_status("LikelihoodState", likelihood_ok);
  all_ok = likelihood_ok && all_ok;

  std::cout << std::setprecision(17)
            << "Maximum absolute difference: " << max_abs_diff << "\n";

  if (!all_ok) {
    std::cerr << "\nFirst mismatch\n";
    std::cerr << "--------------\n";
    std::cerr << first_failure;
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 generated cycle equivalence regression\n";
  return EXIT_SUCCESS;
}
