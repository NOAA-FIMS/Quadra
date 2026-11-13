/**
 * @file scenarios.hpp
 * @brief Validated scenario definitions shared by point and uncertainty forecasts.
 *
 * Exactly five supported rules describe the supplied 1975..2021 case and 2022..2031 horizon.
 * Historical SPR is evaluated with the historical year selectivity blocks. The uncertainty path
 * recalibrates targets and fleet proportions under each parameter perturbation; freezing fitted
 * numerical rates would define a different policy experiment.
 */
#pragma once
#include "projection.hpp"
#include "csv.hpp"
#include <fstream>
#include <sstream>
namespace swordfish {
/// Parsed scenario identifier, target rule, value and reference-year window for the supported ten-year forecast.
/// @details id is consecutive from one, name is a short unique identifier, and target selects one supported rule. value is a fraction for fractional rules and zero for derived rules; first/last are inclusive historical windows or zero when unused.
struct ProjectionScenario {
  int id,first,last;
  std::string name,target;
  double value;
};
/// Read and validate exactly five unquoted scenario records, rejecting malformed numeric tokens, duplicate names and unsupported targets/reference years.
/// @param path Input file path; records must match the documented version/schema.
/// @return The typed result described above, with the module units and indexing conventions.
inline std::vector<ProjectionScenario> read_projection_scenarios(const std::string& path) {
  std::ifstream in(path);std::string line;
  if(!csv_line(in,line)||line!="scenario,name,target,value,reference_first_year,reference_last_year")throw std::runtime_error("Projection config header");
  std::vector<ProjectionScenario> out;
  while(csv_line(in,line)) {
    const auto c=csv_fields(line);if(c.size()!=6)throw std::runtime_error("Projection config columns");
    ProjectionScenario row{csv_integer(c[0]),csv_integer(c[4]),csv_integer(c[5]),c[1],c[2],csv_number(c[3])};
    if(row.id!=int(out.size()+1)||row.name.empty()||!std::isfinite(row.value))throw std::runtime_error("Projection config ID, name or value");
    if(row.name.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_")!=std::string::npos||row.name.size()>80)throw std::invalid_argument("Scenario name must be a short identifier");
    for(const auto& previous:out)if(previous.name==row.name)throw std::invalid_argument("Duplicate scenario name");
    const bool fraction=row.target=="equilibrium_SPR"||row.target=="equilibrium_SSB_fraction_of_dynamic_mean";
    const bool historical=row.target=="equilibrium_SSB_fraction_of_dynamic_mean"||row.target=="mean_historical_SPR"||row.target=="mean_historical_annual_F";
    if(!fraction&&row.target!="mean_historical_SPR"&&row.target!="mean_historical_annual_F"&&row.target!="maximum_equilibrium_yield")throw std::invalid_argument("Unsupported scenario target");
    if((fraction&&!(row.value>0&&row.value<1))||(!fraction&&row.value!=0))throw std::invalid_argument("Invalid target value");
    if((historical&&(row.first<1975||row.last>2021||row.first>row.last))||(!historical&&(row.first!=0||row.last!=0)))throw std::invalid_argument("Invalid scenario reference years");
    out.push_back(row);
  }
  if(out.size()!=5)throw std::runtime_error("Expected five projection scenarios");return out;
}
/// Reject unsupported assessment dimensions and changes among 2019-2021 selectivity curves before using a common terminal curve.
/// @param data Validated original-input schema for the supported swordfish assessment.
/// @param p Native canonical-name parameter map; all required fixed and estimated values must be present.
inline void check_projection_selectivity(const InputData& data,const ParameterValues<double>& p) {
  if(data.dimensions.first_year!=1975||data.dimensions.last_year!=2021||data.dimensions.seasons!=4||data.dimensions.subseasons!=2||data.dimensions.plus_age!=15||data.dimensions.areas!=1)throw std::invalid_argument("Projection implementation supports the 1975-2021 swordfish assessment only");
  const auto curves=length_selectivity(data,p,2021);
  for(int y=2019;y<=2020;++y)if(curves!=length_selectivity(data,p,y))throw std::runtime_error("Projection selectivity requires explicit averaging for changing terminal curves");
}
/// Scenario definition paired with its recalculated numerical target and equilibrium fishing solution.
/// @details target_value has the units implied by definition.target: SSB tonnes, raw SPR or annual F. The MSY rule does not use a numerical input target. equilibrium contains the resulting fishing multiplier and stationary summaries.
struct ResolvedScenario {ProjectionScenario definition;double target_value;BenchmarkResult equilibrium;};
/// Recalculate dynamic-baseline, historical-SPR, SPR, MSY and status-quo targets from a fitted history and solve their fishing multipliers.
/// @param data Validated original-input schema for the supported swordfish assessment.
/// @param p Native canonical-name parameter map; all required fixed and estimated values must be present.
/// @param dev Year-to-log-recruitment-deviation map; the caller enforces the main-period constraint.
/// @param history Historical trajectory fitted to the original catch series.
/// @param unfished Counterfactual history with fishing removed and recruitment feedback recalculated.
/// @param definitions Five validated scenario rules; each rule is recalibrated for the supplied parameter values.
/// @return The typed result described above, with the module units and indexing conventions.
inline std::vector<ResolvedScenario> resolve_projection_scenarios(const InputData& data,const ParameterValues<double>& p,
    const Deviations<double>& dev,const PopulationRun<double>& history,const PopulationRun<double>& unfished,
    const std::vector<ProjectionScenario>& definitions) {
  Benchmarks b(data,p,dev,2019,2021);const double b0=b(0).ssb;
  std::vector<ResolvedScenario> out;
  for(const auto& def:definitions) {
    double value=def.value;BenchmarkResult solution{};
    if(def.target=="equilibrium_SSB_fraction_of_dynamic_mean"||def.target=="mean_historical_SPR"||def.target=="mean_historical_annual_F")
      if(def.first<data.dimensions.first_year||def.last>data.dimensions.last_year||def.first>def.last)throw std::invalid_argument("Projection reference years");
    if(def.target=="equilibrium_SSB_fraction_of_dynamic_mean") {
      double mean=0;for(int y=def.first;y<=def.last;++y)mean+=population_row(unfished,data,y,2).ssb/(def.last-def.first+1);
      value*=mean;solution=b.target(value/b0,true);
    } else if(def.target=="mean_historical_SPR") {
      value=0;for(int y=def.first;y<=def.last;++y)value+=historical_spr(data,p,history,y)/(def.last-def.first+1);
      solution=b.target(value,false);
    } else if(def.target=="equilibrium_SPR")solution=b.target(value,false);
    else if(def.target=="maximum_equilibrium_yield")solution=b.msy();
    else if(def.target=="mean_historical_annual_F") {
      value=0;for(int y=def.first;y<=def.last;++y)value+=annual_f_units4(data,population_row(history,data,y,1).beginning,population_row(history,data,y,4).removals.survivors)/(def.last-def.first+1);
      double lo=0,hi=1;while(b(hi).annual_f<value){hi*=2;if(hi>100)throw std::runtime_error("F target not bracketed");}
      for(int i=0;i<70;++i){const double mid=(lo+hi)/2;if(b(mid).annual_f<value)lo=mid;else hi=mid;}solution=b((lo+hi)/2);
    } else throw std::runtime_error("Unsupported projection target");
    out.push_back({def,value,solution});
  }return out;
}
/// Resolved scenario and its complete deterministic annual forecast trajectory.
/// @details The resolved rule and its ten annual rows belong to one parameter evaluation. Reusing only its fitted multiplier during uncertainty propagation would freeze the policy instead of recalibrating it.
struct ScenarioProjection {ResolvedScenario scenario;std::vector<ProjectionRow> years;};
/// Rebuild historical state, dynamic baseline, fishing mixture, targets and all five forecasts for one native parameter/deviation pair.
/// @param data Validated original-input schema for the supported swordfish assessment.
/// @param p Native canonical-name parameter map; all required fixed and estimated values must be present.
/// @param dev Year-to-log-recruitment-deviation map; the caller enforces the main-period constraint.
/// @param definitions Five validated scenario rules; each rule is recalibrated for the supplied parameter values.
/// @return The typed result described above, with the module units and indexing conventions.
inline std::vector<ScenarioProjection> evaluate_projection_scenarios(const InputData& data,const ParameterValues<double>& p,
    const Deviations<double>& dev,const std::vector<ProjectionScenario>& definitions) {
  check_projection_selectivity(data,p);
  const auto history=run_population(data,p,dev),unfished=dynamic_unfished(data,p,dev);
  const auto pattern=mean_fishing_pattern(data,history,2019,2021);
  const auto january=age_on_january_first(population_row(history,data,2021,4).removals.survivors);
  std::vector<ScenarioProjection> out;
  for(const auto& row:resolve_projection_scenarios(data,p,dev,history,unfished,definitions))
    out.push_back({row,project(data,p,january,pattern,row.equilibrium.scale,10)});
  return out;
}
}
