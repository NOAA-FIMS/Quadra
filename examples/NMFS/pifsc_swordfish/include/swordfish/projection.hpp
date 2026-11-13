/**
 * @file projection.hpp
 * @brief Historical no-fishing counterfactual and deterministic forecasts.
 *
 * Forecast states start in January after aging final historical survivors. Future recruitment
 * follows expected Beverton-Holt recruitment without future deviations or log-bias adjustment.
 * Holding apical rates constant need not hold reported pooled annual F constant, because the age
 * composition changes. These trajectories contain no future process uncertainty.
 */
#pragma once
#include "benchmarks.hpp"
namespace swordfish {
// Rerun the historical recruitment process with fishing removed, including initial F.
/// Rerun historical dynamics with all catch and initial fishing removed, retaining recruitment deviations, bias adjustment and counterfactual spawning feedback.
/// @param data Validated original-input schema for the supported swordfish assessment.
/// @param p Native canonical-name parameter map; all required fixed and estimated values must be present.
/// @param deviations Historical fitted log recruitment deviations, including early age-structure effects.
/// @return The typed result described above, with the module units and indexing conventions.
inline PopulationRun<double> dynamic_unfished(const InputData& data,const ParameterValues<double>& p,
                                              const Deviations<double>& deviations) {
  auto no_fishing=data;for(auto& row:no_fishing.catches)row.value=0;
  auto values=p;for(auto& item:values)if(item.first.rfind("InitF_seas_",0)==0)item.second=0;
  return run_population(no_fishing,values,deviations);
}
/// One projected year's SSB, expected recruitment, January biomass, dead yield and realized pooled annual F.
/// @details year is calendar year; ssb/biomass/yield are tonnes and recruitment is thousands of fish. annual_f is recalculated using transient ages 1..10 abundance. These are deterministic expectations, not individual stochastic outcomes.
struct ProjectionRow {int year;double ssb=0,recruitment=0,biomass=0,yield=0,annual_f=0;};
/// Four seasons by 19 fleets of annual apical fishing rates.
/// @details Rows are the four seasons and columns are the 19 fishing fleets. Entries are annual apical mortality rates even though their application duration is one quarter.
using FishingPattern=std::array<std::array<double,nfishing>,4>;
/// Average fitted fleet/season apical rates over checked historical reference years.
/// @param data Validated original-input schema for the supported swordfish assessment.
/// @param history Historical trajectory fitted to the original catch series.
/// @param first First historical reference year, inclusive.
/// @param last Last historical reference year, inclusive.
/// @return The typed result described above, with the module units and indexing conventions.
inline FishingPattern mean_fishing_pattern(const InputData& data,const PopulationRun<double>& history,int first,int last) {
  if(first<data.dimensions.first_year||last>data.dimensions.last_year||first>last)throw std::invalid_argument("Fishing pattern years");
  FishingPattern out{};for(int y=first;y<=last;++y)for(int s=0;s<4;++s)for(int f=0;f<nfishing;++f)
    out[s][f]+=population_row(history,data,y,s+1).hybrid.apical_f[f]/(last-first+1);
  return out;
}
/// Calculate raw per-recruit SPR using that historical year's fishing rates and selectivity blocks.
/// @param data Validated original-input schema for the supported swordfish assessment.
/// @param p Native canonical-name parameter map; all required fixed and estimated values must be present.
/// @param history Historical trajectory fitted to the original catch series.
/// @param year Calendar year used to select observations, deviations or active control blocks.
/// @return Raw SPR, not the report display quantity 1-SPR.
inline double historical_spr(const InputData& data,const ParameterValues<double>& p,const PopulationRun<double>& history,int year) {
  const auto curves=length_selectivity(data,p,year);std::array<SeasonalBiology<double>,4> biology;
  for(int s=0;s<4;++s)biology[s]=seasonal_biology(data,p,curves,s+1);
  FishingPattern zero{};const auto rates=mean_fishing_pattern(data,history,year,year);
  return equilibrium(data,p,biology,rates,1.).ssb/equilibrium(data,p,biology,zero,1.).ssb;
}
/// Calculate pooled numbers-weighted annual F over ages 1-10 from January and end-year abundance, matching SS3 units 4.
/// @param data Validated original-input schema for the supported swordfish assessment.
/// @param beginning January abundance in thousands of fish before that year mortality.
/// @param end December survivors before January aging, on the same scale as beginning.
/// @return Pooled annual fishing mortality log(natural-only survivors / actual survivors) over ages 1..10.
inline double annual_f_units4(const InputData& data,const State<double>& beginning,const State<double>& end) {
  double natural=0,actual=0;for(int sex=0;sex<2;++sex)for(int age=1;age<=10;++age){const int j=cell(sex,age);natural+=beginning[j]*std::exp(-data.biology.natural_mortality[j]);actual+=end[j];}
  if(!(actual>0&&natural>0))throw std::invalid_argument("Annual F requires positive ages 1-10 abundance");
  return std::log(natural/actual);
}
/// Construct all four seasonal biological factors under terminal-year selectivity.
/// @param data Validated original-input schema for the supported swordfish assessment.
/// @param p Native canonical-name parameter map; all required fixed and estimated values must be present.
/// @return The typed result described above, with the module units and indexing conventions.
inline std::array<SeasonalBiology<double>,4> terminal_biology(const InputData& data,const ParameterValues<double>& p) {
  const auto curves=length_selectivity(data,p,data.dimensions.last_year);std::array<SeasonalBiology<double>,4> out;
  for(int s=0;s<4;++s)out[s]=seasonal_biology(data,p,curves,s+1);return out;
}
/// Advance a January state through deterministic future seasons at constant apical fishing rates, with BH expected recruitment and no future process deviations.
/// @param data Validated original-input schema for the supported swordfish assessment.
/// @param p Native canonical-name parameter map; all required fixed and estimated values must be present.
/// @param january Starting January abundance in thousands of fish; supplied by value so the caller state is preserved.
/// @param pattern Fixed fleet/season relative apical fishing-rate pattern for the projection.
/// @param scale Nonnegative common multiplier on the reference apical-rate pattern.
/// @param years Number of projected calendar years; supported helper range is 1..200.
/// @return The typed result described above, with the module units and indexing conventions.
/// @details The input January state is copied. No future recruitment errors are simulated. Annual F is recomputed from transient age weights, so it can vary even when pattern and scale are fixed.
inline std::vector<ProjectionRow> project(const InputData& data,const ParameterValues<double>& p,
                    State<double> january,const FishingPattern& pattern,double scale,int years) {
  if(!(scale>=0&&std::isfinite(scale))||years<1||years>200)throw std::invalid_argument("Projection scale or duration");
  const auto biology=terminal_biology(data,p);auto rates=pattern;for(auto& season:rates)for(auto& f:season)f*=scale;
  FishingPattern zero{};const double r0=std::exp(p.at("SR_LN(R0)"));
  const double b0=equilibrium(data,p,biology,zero,r0).ssb;
  std::vector<ProjectionRow> out;
  for(int y=data.dimensions.last_year+1;y<=data.dimensions.last_year+years;++y){
    ProjectionRow row;row.year=y;const auto beginning=january;
    for(int j=0;j<32;++j)row.biomass+=january[j]*biology[0].weight_begin[j];
    auto state=january;
    for(int s=0;s<4;++s){if(s==2)settle(state,row.recruitment,p.at("FracFemale_GP_1"));
      const auto fishing=fishing_at_age(rates[s],biology[s]);
      if(s==1){row.ssb=seasonal_ssb(data,state,fishing,biology[s],s+1);
        // Published forecast: deterministic BH recruitment, no future deviations or log-bias adjustment.
        row.recruitment=beverton_holt(row.ssb,r0,b0,p.at("SR_BH_steep"))*std::exp(p.at("SR_regime"));}
      for(int j=0;j<32;++j){double z=data.biology.natural_mortality[j];for(const auto& f:fishing)z+=f[j];
        for(int f=0;f<nfishing;++f)row.yield+=state[j]*rates[s][f]*biology[s].selected_biomass[f][j]*exposure(z,.25);}
      state=advance_season(state,data.biology.natural_mortality,fishing,.25).survivors;
    }
    row.annual_f=annual_f_units4(data,beginning,state);january=age_on_january_first(state);out.push_back(row);
  }return out;
}
}
