#include "swordfish/projection.hpp"
#include "swordfish/reference_io.hpp"
#include "swordfish_input.hpp"
#include <iostream>
using namespace swordfish;
/// Test helper for test projection: near. See the typed signature and assertions in the linked test source.
void near(double actual,double expected,const std::string& label){if(std::abs(actual-expected)>1e-8*std::max(1.,std::abs(expected)))throw std::runtime_error(label);}
/// Run test projection; return zero on success and nonzero on a failed check or invalid input.
int main(){try{auto data=imported_swordfish_data();auto native=read_reference_parameters(data,"build/biology_fitted_parameters.txt","build/population_fitted_deviations.txt");
 const auto& p=native.first;auto dev=native.second;const auto biology=terminal_biology(data,p);FishingPattern zero{};const double r0=std::exp(p.at("SR_LN(R0)"));const auto virgin=equilibrium(data,p,biology,zero,r0);
 for(const auto& row:project(data,p,virgin.beginning[0],zero,1,20)){near(row.ssb,virgin.ssb,"Unfished stationary SSB");near(row.recruitment,r0,"Unfished stationary recruitment");near(row.yield,0,"Unfished catch");near(row.annual_f,0,"Unfished F");}
 // Removing all fishing and all recruitment variability must preserve virgin equilibrium.
 auto constant_data=data;constant_data.recruitment_advanced[8]=0;for(auto& row:dev)row.second=0;const auto constant=dynamic_unfished(constant_data,p,dev);
 for(int y=1975;y<=2021;++y){near(population_row(constant,constant_data,y,2).ssb,virgin.ssb,"Stationary dynamic B0");near(population_row(constant,constant_data,y,3).recruits,r0,"Stationary dynamic recruitment");}
 const auto history=run_population(data,p,native.second);auto pattern=mean_fishing_pattern(data,history,2019,2021);Benchmarks benchmarks(data,p,native.second,2019,2021);const double scale=1.5;const auto target=benchmarks(scale);auto rates=pattern;for(auto& season:rates)for(auto& f:season)f*=scale;
 const auto equilibrium_state=equilibrium(data,p,biology,rates,target.recruitment);
 for(const auto& row:project(data,p,equilibrium_state.beginning[0],pattern,scale,40)){near(row.ssb,target.ssb,"Fished stationary SSB");near(row.recruitment,target.recruitment,"Fished stationary recruitment");near(row.yield,target.yield,"Fished stationary yield");near(row.annual_f,target.annual_f,"Fished stationary annual F");}
 bool rejected=false;try{project(data,p,virgin.beginning[0],pattern,-1,10);}catch(const std::invalid_argument&){rejected=true;}if(!rejected)throw std::runtime_error("Negative F accepted");
 std::cout<<"Dynamic unfished and projection equilibrium, recruitment, yield and F checks pass\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
