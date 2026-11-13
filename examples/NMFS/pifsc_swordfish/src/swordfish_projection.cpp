/**
 * @file swordfish_projection.cpp
 * @details
 * Deterministic scenario and counterfactual reporting. Five shared scenario definitions use
 * fitted history, not the published bootstrap medians as calibration targets. Alternate readings
 * of conflicting source descriptions are explicitly reported as sensitivities.
 */
#include "swordfish/scenarios.hpp"
#include "swordfish/reference_io.hpp"
#include "swordfish_input.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
using namespace swordfish;
/// Run swordfish projection; return zero on success and nonzero on a failed check or invalid input.
int main(){try{
 auto data=imported_swordfish_data();const auto native=read_fitted_parameters(data,"build/quadra_fitted_parameters.csv","build/biology_fitted_parameters.txt","build/population_fitted_deviations.txt");
 const auto& p=native.first;const auto& dev=native.second;const auto history=run_population(data,p,dev),unfished=dynamic_unfished(data,p,dev);
 const std::filesystem::path folder="build/projections";std::filesystem::create_directories(folder);
 std::ofstream dynamic(folder/"dynamic_unfished.csv"),numbers(folder/"dynamic_unfished_numbers.csv"),spr_history(folder/"historical_SPR.csv");
 spr_history<<"year,raw_SPR\n"<<std::setprecision(16);
 dynamic<<"year,spawning_biomass_tonnes,recruitment_thousands,fished_spawning_biomass_tonnes,depletion\n"<<std::setprecision(16);
 numbers<<"year,sex,age,numbers_SS3_annual_thousands\n"<<std::setprecision(16);
 for(int y=data.dimensions.first_year;y<=data.dimensions.last_year;++y){const double b=population_row(unfished,data,y,2).ssb,f=population_row(history,data,y,2).ssb;
  spr_history<<y<<','<<historical_spr(data,p,history,y)<<'\n';
  dynamic<<y<<','<<b<<','<<population_row(unfished,data,y,3).recruits<<','<<f<<','<<f/b<<'\n';
  for(int sex=0;sex<2;++sex)for(int age=0;age<nage;++age)numbers<<y<<','<<sex+1<<','<<age<<','<<(age==0?population_row(unfished,data,y,3).recruits*(sex==0?p.at("FracFemale_GP_1"):1-p.at("FracFemale_GP_1")):population_row(unfished,data,y,1).beginning[cell(sex,age)])<<'\n';
 }
 check_projection_selectivity(data,p);
 const auto pattern=mean_fishing_pattern(data,history,2019,2021);Benchmarks b(data,p,dev,2019,2021);
 const auto january=age_on_january_first(population_row(history,data,2021,4).removals.survivors);
 std::ofstream output(folder/"deterministic_projections.csv"),specs(folder/"resolved_scenarios.csv"),pattern_output(folder/"fishing_pattern.csv");
 output<<"scenario,name,year,spawning_biomass_tonnes,recruitment_thousands,biomass_january_tonnes,dead_catch_tonnes,annual_F_SS3_units4_ages1_10\n"<<std::setprecision(16);
 specs<<"scenario,name,target,target_value,F_multiplier,equilibrium_SPR,equilibrium_SSB_tonnes,equilibrium_annual_F,equilibrium_yield_tonnes\n"<<std::setprecision(16);
 pattern_output<<"season,fleet,mean_apical_F_2019_2021\n"<<std::setprecision(16);for(int s=0;s<4;++s)for(int f=0;f<nfishing;++f)pattern_output<<s+1<<','<<f+1<<','<<pattern[s][f]<<'\n';
 const auto definitions=read_projection_scenarios("data/projection_scenarios.csv");
 const auto scenarios=resolve_projection_scenarios(data,p,dev,history,unfished,definitions);
 for(const auto& resolved:scenarios){const auto& def=resolved.definition;const auto& solution=resolved.equilibrium;
   specs<<def.id<<','<<def.name<<','<<def.target<<','<<resolved.target_value<<','<<solution.scale<<','<<solution.spr<<','<<solution.ssb<<','<<solution.annual_f<<','<<solution.yield<<'\n';
   for(const auto& row:project(data,p,january,pattern,solution.scale,10))output<<def.id<<','<<def.name<<','<<row.year<<','<<row.ssb<<','<<row.recruitment<<','<<row.biomass<<','<<row.yield<<','<<row.annual_f<<'\n';
 }
 // Retain alternate readings of the conflicting source descriptions as explicit sensitivities.
 std::ofstream sensitivities(folder/"projection_sensitivities.csv"),sensitivity_specs(folder/"resolved_sensitivities.csv");
 sensitivities<<"sensitivity,baseline_scenario,year,spawning_biomass_tonnes,dead_catch_tonnes,annual_F_SS3_units4_ages1_10\n"<<std::setprecision(16);
 sensitivity_specs<<"sensitivity,baseline_scenario,F_multiplier,equilibrium_SPR,equilibrium_SSB_tonnes,equilibrium_annual_F\n"<<std::setprecision(16);
 double alternate_spr=0;
 for(int y=1998;y<=2000;++y)alternate_spr+=historical_spr(data,p,history,y)/3;
 struct Sensitivity {std::string name;int baseline;BenchmarkResult result;};
 for(const auto& item:std::vector<Sensitivity>{{"Equilibrium20_instead_of_dynamic20",1,b.target(.2,true)},{"Reference1998_2000_table_label",2,b.target(alternate_spr,false)},{"StatusQuo_direct_mean_apical_F",5,b(1)}}){
   sensitivity_specs<<item.name<<','<<item.baseline<<','<<item.result.scale<<','<<item.result.spr<<','<<item.result.ssb<<','<<item.result.annual_f<<'\n';
   for(const auto& row:project(data,p,january,pattern,item.result.scale,10))sensitivities<<item.name<<','<<item.baseline<<','<<row.year<<','<<row.ssb<<','<<row.yield<<','<<row.annual_f<<'\n';
 }
 if(!sensitivities||!sensitivity_specs)throw std::runtime_error("Sensitivity output failure");
 if(scenarios.size()!=5||!dynamic||!numbers||!output||!specs||!pattern_output||!spr_history)throw std::runtime_error("Projection output/configuration failure");
 std::cout<<"Wrote historical dynamic unfished outputs and five 2022-2031 deterministic scenarios\n";
 return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
