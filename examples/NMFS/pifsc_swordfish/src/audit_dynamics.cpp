/**
 * @file audit_dynamics.cpp
 * @details
 * Conditional replay audit: supplied report fishing/recruitment isolates seasonal arithmetic
 * from estimation. This is intentionally not an independent fit; audit_population and the
 * assessment executable provide the independent original-catch path.
 */
#include "swordfish/biology.hpp"
#include "swordfish_biology_input.hpp"
#include <algorithm>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

using namespace swordfish;
/// Accumulate finite benchmark comparison counts and worst absolute/tolerance-scaled errors for this audit.
struct Check {
  size_t count = 0;
  double maximum_absolute_error = 0, maximum_tolerance_fraction = 0;
  /// Compare one finite model/reference pair using absolute plus relative tolerance, updating the audit summary.
  void compare(double value, double expected, double atol, double rtol) {
    if (!std::isfinite(value) || !std::isfinite(expected)) throw std::runtime_error("Nonfinite audit value");
    const double error = std::abs(value-expected);
    maximum_absolute_error = std::max(maximum_absolute_error,error);
    maximum_tolerance_fraction = std::max(maximum_tolerance_fraction,error/(atol+rtol*std::abs(expected)));
    ++count;
  }
};
template <class T> void read(std::istream& in, T& x) {
  /// Read one typed comparison-fixture token and throw if the fixture is truncated or malformed.
  if (!(in >> x)) throw std::runtime_error("Truncated or invalid replay fixture");
}
/// Run audit dynamics; return zero on success and nonzero on a failed check or invalid input.
int main(int argc, char** argv) {
  try {
    if (argc != 4) throw std::runtime_error("Usage: audit_dynamics replay.txt audit.csv fitted_parameters.txt");
    const auto data = imported_swordfish_biology_data();
    auto fitted=initial_values<double>(data);
    std::ifstream parameters(argv[3]);std::string parameter_magic;size_t np;
    if (!(parameters>>parameter_magic>>np) || parameter_magic!="SWORDFISH_FITTED_PARAMETERS_V1" || np!=fitted.size())
      throw std::runtime_error("Fitted biology parameter header");
    std::map<std::string,bool> seen;
    for (size_t i=0;i<np;++i) {
      std::string name;double value;
      if (!(parameters>>name>>value) || !fitted.count(name) || !std::isfinite(value) || seen[name])
        throw std::runtime_error("Fitted biology parameter mapping");
      fitted.at(name)=value;seen[name]=true;
    }
    std::string trailing;if(parameters>>trailing)throw std::runtime_error("Unread fitted parameters");
    std::ifstream input(argv[1]);
    if (!input) throw std::runtime_error("Cannot open replay fixture");
    std::string magic; read(input,magic);
    if (magic != "SWORDFISH_CONDITIONAL_REPLAY_V2") throw std::runtime_error("Unsupported fixture");
    int first,last; read(input,first); read(input,last);
    if (first != data.dimensions.first_year || last != data.dimensions.last_year)
      throw std::runtime_error("Fixture year range differs from model inputs");
    std::array<double,nage> mature_weight{};
    for (auto& x : mature_weight) read(input,x);
    std::array<double,nage> calculated_mature_weight{};
    for(int a=0;a<nage;++a)calculated_mature_weight[a]=mature_weight_at_age(data,fitted,a);
    State<double> state{};
    std::map<std::string,Check> checks;
    std::vector<std::vector<double>> independent_curves;
    std::ofstream timing(std::filesystem::path(argv[2]).parent_path()/"spawning_timing_comparison.csv");
    if (!timing) throw std::runtime_error("Cannot open timing diagnostic output");
    timing << "year,ss3_report_ssb,ss33020_conditional_ssb,calendar_conditional_ssb\n" << std::setprecision(12);
    for (int y = first; y <= last; ++y) for (int s = 1; s <= 4; ++s) {
      int fy,fs; read(input,fy); read(input,fs);
      if (fy != y || fs != s) throw std::runtime_error("Replay order");
      State<double> begin{}, mid{}, weights{};
      std::vector<State<double>> fishing(nfishing), observed_catch(nfishing);
      for (int sex = 0; sex < nsex; ++sex) {
        for (int a = 0; a < nage; ++a) read(input,begin[cell(sex,a)]);
        for (int a = 0; a < nage; ++a) read(input,mid[cell(sex,a)]);
        for (int a = 0; a < nage; ++a) read(input,weights[cell(sex,a)]);
        for (int f = 0; f < nfishing; ++f) {
          for (int a = 0; a < nage; ++a) read(input,fishing[f][cell(sex,a)]);
          for (int a = 0; a < nage; ++a) read(input,observed_catch[f][cell(sex,a)]);
        }
      }
      double recruits,bio,summary_bio,ssb;
      read(input,recruits); read(input,bio); read(input,summary_bio); read(input,ssb);
      if(s==1)independent_curves=length_selectivity(data,fitted,y);
      for(int f=0;f<nfishing;++f) {
        double apical_f;read(input,apical_f);
        for(int sex=0;sex<nsex;++sex)for(int age=0;age<nage;++age) {
          const auto selected=selected_at_age(data,fitted,independent_curves[f],sex,age,(s-1)*0.25+0.125);
          const double rate=apical_f*selected.numbers;
          checks["independent_f_at_age"].compare(rate,fishing[f][cell(sex,age)],0.000001,0.00005);
          fishing[f][cell(sex,age)]=rate;
        }
      }
      if (y == first && s == 1) state = begin; // only initial condition from SS3
      if (s == 3) settle(state,recruits,data.biology.female_fraction);
      for (int j = 0; j < nsex*nage; ++j)
        checks["begin_numbers"].compare(state[j],begin[j],0.001,0.00002);
      const double duration = data.dimensions.season_months[s-1]/12;
      const auto result = advance_season(state,data.biology.natural_mortality,fishing,duration);
      const auto midpoint = advance_season(state,data.biology.natural_mortality,fishing,duration/2);
      double calculated_bio = 0, calculated_summary = 0;
      for (int j = 0; j < nsex*nage; ++j) {
        checks["midseason_numbers"].compare(midpoint.survivors[j],mid[j],0.001,0.00002);
        double balance = result.survivors[j] + result.natural_deaths[j];
        for (int f = 0; f < nfishing; ++f) {
          balance += result.catch_numbers[f][j];
          checks["catch_numbers"].compare(result.catch_numbers[f][j],observed_catch[f][j],0.0001,0.00003);
        }
        checks["mass_balance"].compare(balance,state[j],1e-10,1e-12);
        const double biological_weight=weight_at_age(data,fitted,j/nage,j%nage,(s-1)*0.25);
        checks["independent_body_weight"].compare(biological_weight,weights[j],0.001,0.00002);
        calculated_bio += state[j]*biological_weight;
        if (j % nage >= 1) calculated_summary += state[j]*biological_weight;
      }
      checks["total_biomass"].compare(calculated_bio,bio,0.7,0.00002);
      checks["age1plus_biomass"].compare(calculated_summary,summary_bio,0.7,0.00002);
      if (ssb >= 0) {
        const double elapsed = (data.dimensions.spawn_month-1)/12 - (s-1)*duration;
        const double compatible = spawning_biomass_ss33020(state,data.biology.natural_mortality,
                       fishing,calculated_mature_weight,elapsed,duration);
        const double calendar = spawning_biomass(state,data.biology.natural_mortality,
                       fishing,calculated_mature_weight,elapsed);
        checks["spawning_biomass_ss33020"].compare(compatible,ssb,0.5,0.00002);
        timing << y << ',' << ssb << ',' << compatible << ',' << calendar << '\n';
      }
      state = s == 4 ? age_on_january_first(result.survivors) : result.survivors;
    }
    for(int a=0;a<nage;++a)
      checks["independent_mature_weight"].compare(calculated_mature_weight[a],mature_weight[a],0.001,0.00003);
    // Check ordinary age mean lengths using the fitted-report growth parameters.
    // The control/report fixed male L1 discrepancy is deliberately not hidden.
    int nbiology; read(input,nbiology);
    struct LengthRow { int sex,season,age; double begin,mid; };
    std::vector<LengthRow> lengths(nbiology);
    for (auto& row : lengths) {
      read(input,row.sex); read(input,row.season); read(input,row.age);
      read(input,row.begin); read(input,row.mid);
    }
    int ngrowth; read(input,ngrowth);
    if (ngrowth != 2) throw std::runtime_error("Growth sex count");
    std::array<Growth,2> runtime_growth;
    for (int i = 0; i < ngrowth; ++i) {
      int sex; read(input,sex);
      if (sex != i+1) throw std::runtime_error("Growth sex order");
      auto& g = runtime_growth[i];
      for (auto field : {&g.age1,&g.age2,&g.length1,&g.length2,&g.k,&g.cv_young,&g.cv_old,
                         &g.weight_coefficient,&g.weight_power}) read(input,*field);
    }
    std::string extra;
    if (input >> extra) throw std::runtime_error("Unread fixture values");
    for (const auto& row : lengths) {
      if (row.age == nage-1) continue; // aggregate plus-group growth is not yet implemented
      const double birth_fraction = (data.biology.settlement_month-1)/12;
      const double age = row.age+(row.season-1)*0.25-birth_fraction;
      const auto& g = runtime_growth[row.sex-1];
      const double lo = data.population_length_bins.front();
      checks["mean_length_begin"].compare(mean_length(age,g,lo),row.begin,0.001,0.00001);
      checks["mean_length_mid"].compare(mean_length(age+0.125,g,lo),row.mid,0.001,0.00001);
    }
    std::ofstream out(argv[2]);
    if (!out) throw std::runtime_error("Cannot open audit output");
    out << "check,comparisons,max_absolute_error,max_tolerance_fraction,passed\n" << std::setprecision(12);
    bool passed = true;
    for (const auto& [name,check] : checks) {
      const bool ok = check.count > 0 && check.maximum_tolerance_fraction <= 1;
      out << name << ',' << check.count << ',' << check.maximum_absolute_error << ','
          << check.maximum_tolerance_fraction << ',' << (ok ? "true" : "false") << '\n';
      std::cout << name << ": " << (ok ? "PASS" : "FAIL") << " (" << check.count
                << ", max error " << check.maximum_absolute_error << ")\n";
      passed = passed && ok;
    }
    if (!out || !timing) throw std::runtime_error("Could not write audit output");
    return passed ? 0 : 1;
  } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 2; }
}
