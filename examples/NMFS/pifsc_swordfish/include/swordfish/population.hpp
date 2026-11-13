/**
 * @file population.hpp
 * @brief Independent historical population reconstruction.
 *
 * Catches and controls are model inputs; Report.sso trajectories are comparison fixtures only. The
 * hybrid method intentionally performs a finite number of SS3 tuning updates, not an arbitrarily
 * converged catch root. Historical recruitment includes fitted deviations and the supplied bias
 * ramp. Seasonal biology caches are local to an evaluation and cannot retain AD values across
 * tapes.
 */
#pragma once
#include "biology.hpp"

namespace swordfish {
/// Map recruitment years to fitted log recruitment deviations, including early age-structure deviations.
/// @details Keys include early age-structure and historical main recruitment years. Values are additive log deviations. The map itself does not impose the main-period zero-sum constraint; SwordfishModel::decode reconstructs its dependent final year.
template<class T> using Deviations=std::map<int,T>;
/// Per-season selected numbers, selected weight, beginning weights and mature weights for mortality and biomass calculations.
/// @details selectivity is dimensionless number selection, selected_biomass integrates selected kg weight, and weight_begin is beginning-season kg weight. mature_weight is the female maturity/fecundity contribution used at the supported spawning integration time.
template<class T> struct SeasonalBiology {
  std::vector<State<T>> selectivity,selected_biomass;
  State<T> weight_begin{};
  std::array<T,nage> mature_weight{};
};
/// Integrate growth, selectivity and mature weight for one of four seasons; produce all fleet/sex/age factors needed by population dynamics.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active evaluation tape.
/// @param data Validated original-input schema for the supported swordfish assessment.
/// @param p Native canonical-name parameter map; all required fixed and estimated values must be present.
/// @param curves Fleet length-selectivity vectors for the chosen year/block combination.
/// @param season One-based quarter, 1..4.
/// @return The typed result described above, with the module units and indexing conventions.
template<class T> SeasonalBiology<T> seasonal_biology(const InputData& data,const ParameterValues<T>& p,
                        const std::vector<std::vector<T>>& curves,int season) {
  if(season<1 || season>4)throw std::invalid_argument("Season index");
  SeasonalBiology<T> b;b.selectivity.resize(nfishing);b.selected_biomass.resize(nfishing);
  std::array<std::vector<T>,2> weights;
  for(int sex=0;sex<2;++sex)for(size_t bin=0;bin<data.population_length_bins.size();++bin)
    weights[sex].push_back(weight_at_length(p,sex,length_midpoint(data,bin)));
  for(int sex=0;sex<2;++sex)for(int a=0;a<nage;++a) {
    const auto key=age_length_key(data,p,sex,a,(season-1)*0.25+0.125);
    b.weight_begin[cell(sex,a)]=weight_at_age(data,p,sex,a,(season-1)*0.25);
    for(int f=0;f<nfishing;++f) {
      T numbers=T(0),biomass=T(0);
      for(size_t bin=0;bin<key.size();++bin) {
        const T selected=key[bin]*curves.at(f).at(bin);
        numbers+=selected;biomass+=selected*weights[sex][bin];
      }
      b.selectivity[f][cell(sex,a)]=numbers;b.selected_biomass[f][cell(sex,a)]=biomass;
    }
  }
  for(int a=0;a<nage;++a)b.mature_weight[a]=mature_weight_at_age(data,p,a);
  return b;
}
/// Multiply fleet apical annual fishing rates by numbers selectivity to obtain fleet-specific sex/age mortality arrays.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active evaluation tape.
/// @param apical Annual apical fishing rates, indexed by season and fleet where applicable.
/// @param b Seasonal selection/weight integration consistent with the rate vectors.
/// @return The typed result described above, with the module units and indexing conventions.
template<class T> std::vector<State<T>> fishing_at_age(const std::array<T,nfishing>& apical,
                                                      const SeasonalBiology<T>& b) {
  std::vector<State<T>> f(nfishing);
  for(int fleet=0;fleet<nfishing;++fleet)for(int j=0;j<32;++j)f[fleet][j]=apical[fleet]*b.selectivity[fleet][j];
  return f;
}
/// Hybrid-catch solution containing fleet apical rates, predicted catches and mortality-at-age arrays.
/// @details apical_f and fishing are annual mortality rates, not quarterly removals. predicted_catch follows each fleet unit code, so summed mixed-unit catches are only used where SS3 deliberately uses that tuning heuristic.
template<class T> struct HybridResult {
  std::array<T,nfishing> apical_f,predicted_catch;
  std::vector<State<T>> fishing;
  /// Explicitly initialize accumulated catch and apical-rate scalars to numerical zero, including AD types.
  HybridResult() {apical_f.fill(T(0));predicted_catch.fill(T(0));}
};
/// Reproduce SS3 method-3 finite catch tuning, including Pope initialization, smooth caps, mixed catch units and the configured finite adjustment count.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active evaluation tape.
/// @param data Validated original-input schema for the supported swordfish assessment.
/// @param numbers Sex-by-age abundance at the start of the calculation, in thousands of fish.
/// @param biology Seasonal selection, weight and maturity factors from the same native parameters.
/// @param observed Fleet catches in their configured units: tonnes or thousands of fish; zero means no observed catch.
/// @param duration Exposure duration in calendar years; a quarter is 0.25.
/// @return The typed result described above, with the module units and indexing conventions.
template<class T> HybridResult<T> hybrid_catch(const InputData& data,const State<T>& numbers,
                 const SeasonalBiology<T>& biology,const std::array<double,nfishing>& observed,double duration) {
  using std::exp;using std::log;
  if(data.f_method!=3 || data.hybrid_iterations<1 || !(duration>0 && duration<=1))
    throw std::invalid_argument("Hybrid F configuration");
  for(int f=0;f<nfishing;++f)
    if(data.fleets.at(f).catch_multiplier || !std::isfinite(observed[f]) || observed[f]<0)
      throw std::invalid_argument("Unsupported catch multiplier or invalid catch");
  HybridResult<T> result;
  bool any=false;
  for(int f=0;f<nfishing;++f)if(observed[f]>0) {
    any=true;T vulnerable=T(0);
    const auto& selected=data.fleets[f].catch_units==1 ? biology.selected_biomass[f]:biology.selectivity[f];
    for(int j=0;j<32;++j)vulnerable+=numbers[j]*exp(T(-data.biology.natural_mortality[j]*duration/2))*selected[j];
    const T harvest=T(observed[f])/(vulnerable+T(0.1*observed[f]));
    const T join=logistic(T(-30)*(harvest-T(0.95)));
    const T bounded=join*harvest+(T(1)-join)*T(0.95);
    result.apical_f[f]=-log(T(1)-bounded)/duration;
  }
  // SS3 F_Tune=5 means a Pope start plus four adjustments. More iterations
  // would define a different objective, even if catches appeared more exact.
  if(any)for(int tune=1;tune<data.hybrid_iterations;++tune) {
    State<T> z{},integral{};
    for(int j=0;j<32;++j) {
      z[j]=T(data.biology.natural_mortality[j]);
      for(int f=0;f<nfishing;++f)z[j]+=result.apical_f[f]*biology.selectivity[f][j];
      integral[j]=exposure(z[j],duration);
    }
    T predicted=T(0);double target=0;
    for(int f=0;f<nfishing;++f)if(observed[f]>0) {
      const auto& selected=data.fleets[f].catch_units==1 ? biology.selected_biomass[f]:biology.selectivity[f];
      for(int j=0;j<32;++j)predicted+=result.apical_f[f]*numbers[j]*selected[j]*integral[j];
      target+=observed[f];
    }
    // Preserve SS3's aggregate adjustment even when fleet catch units differ.
    const T adjustment=T(target)/(predicted+T(0.0001));
    for(int j=0;j<32;++j) {
      const T adjusted=T(data.biology.natural_mortality[j])+adjustment*(z[j]-T(data.biology.natural_mortality[j]));
      integral[j]=exposure(adjusted,duration);
    }
    for(int f=0;f<nfishing;++f)if(observed[f]>0) {
      T vulnerable=T(0);
      const auto& selected=data.fleets[f].catch_units==1 ? biology.selected_biomass[f]:biology.selectivity[f];
      for(int j=0;j<32;++j)vulnerable+=numbers[j]*selected[j]*integral[j];
      const T candidate=T(observed[f])/(vulnerable+T(0.0001));
      const T join=logistic(T(-30)*(candidate-T(0.95*data.max_f)));
      result.apical_f[f]=join*candidate+(T(1)-join)*T(data.max_f);
    }
  }
  result.fishing=fishing_at_age(result.apical_f,biology);
  State<T> integral{};
  for(int j=0;j<32;++j) {
    T z=T(data.biology.natural_mortality[j]);for(const auto& f:result.fishing)z+=f[j];
    integral[j]=exposure(z,duration);
  }
  for(int f=0;f<nfishing;++f) {
    const auto& selected=data.fleets[f].catch_units==1 ? biology.selected_biomass[f]:biology.selectivity[f];
    for(int j=0;j<32;++j)result.predicted_catch[f]+=numbers[j]*selected[j]*result.apical_f[f]*integral[j];
  }
  return result;
}
/// Evaluate the piecewise recruitment-bias ramp for a historical year from the advanced control settings.
/// @param data Validated original-input schema for the supported swordfish assessment.
/// @param year Calendar year used to select observations, deviations or active control blocks.
/// @return Configured bias-ramp fraction for the year.
inline double recruitment_bias(const InputData& data,int year) {
  const auto& a=data.recruitment_advanced;
  if(year<int(a[0]) || year>data.recdev_end)return 0;
  const double y=year;
  if(y<=a[4] || y>=a[7])return 0;
  if(y<a[5])return a[8]*(y-a[4])/(a[5]-a[4]);
  if(y<=a[6])return a[8];
  return a[8]*(a[7]-y)/(a[7]-a[6]);
}
/// Apply a year's fitted deviation and log-bias correction to expected recruitment without changing the deviation map.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active evaluation tape.
/// @param data Validated original-input schema for the supported swordfish assessment.
/// @param p Native canonical-name parameter map; all required fixed and estimated values must be present.
/// @param dev Year-to-log-recruitment-deviation map; the caller enforces the main-period constraint.
/// @param year Calendar year used to select observations, deviations or active control blocks.
/// @param expected Beverton-Holt expected recruitment in thousands of fish, before deviation/bias adjustment.
/// @return The typed result described above, with the module units and indexing conventions.
template<class T> T adjusted_recruitment(const InputData& data,const ParameterValues<T>& p,
                                        const Deviations<T>& dev,int year,const T& expected) {
  using std::exp;
  const auto it=dev.find(year);
  const T deviation=it==dev.end() ? T(0):it->second;
  const T sigma=p.at("SR_sigmaR");
  return expected*exp(deviation-T(0.5*recruitment_bias(data,year))*sigma*sigma);
}
/// Combine beginning numbers, fishing mortality and mature weights to compute spawning biomass under SS3 v3.30.20 timing.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active evaluation tape.
/// @param data Validated original-input schema for the supported swordfish assessment.
/// @param numbers Sex-by-age abundance at the start of the calculation, in thousands of fish.
/// @param fishing Fleet-specific sex-by-age fishing mortality arrays in annual rate units.
/// @param biology Seasonal selection, weight and maturity factors from the same native parameters.
/// @param season One-based quarter, 1..4.
/// @return Female spawning biomass in tonnes, using legacy SS3 timing.
template<class T> T seasonal_ssb(const InputData& data,const State<T>& numbers,
                      const std::vector<State<T>>& fishing,const SeasonalBiology<T>& biology,int season) {
  using std::exp;
  const double start=(season-1)*0.25;
  const double fraction=((data.dimensions.spawn_month-1)/12-start)/0.25;
  T ssb=T(0);
  for(int a=0;a<nage;++a) {
    T z=T(data.biology.natural_mortality[cell(0,a)]);
    for(const auto& f:fishing)z+=f[cell(0,a)];
    ssb+=numbers[cell(0,a)]*exp(-z*fraction)*biology.mature_weight[a];
  }
  return ssb;
}
/// Per-season beginning population states, equilibrium recruitment and spawning biomass for constant seasonal mortality.
/// @details beginning contains each season state after any settlement occurring at that season start. recruitment is constant annual recruitment in thousands; ssb is the spawning-season mature biomass in tonnes. The plus group includes survivors already at age 15.
template<class T> struct Equilibrium {
  std::array<State<T>,4> beginning{};
  T ssb=T(0),recruitment=T(0);
};
/// Solve the stationary seasonal sex/age population for supplied fishing rates and recruitment, including July settlement and the plus group.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active evaluation tape.
/// @param data Validated original-input schema for the supported swordfish assessment.
/// @param p Native canonical-name parameter map; all required fixed and estimated values must be present.
/// @param biology Seasonal selection, weight and maturity factors from the same native parameters.
/// @param apical Annual apical fishing rates, indexed by season and fleet where applicable.
/// @param recruitment Constant equilibrium recruitment in thousands of fish.
/// @return The typed result described above, with the module units and indexing conventions.
template<class T> Equilibrium<T> equilibrium(const InputData& data,const ParameterValues<T>& p,
                     const std::array<SeasonalBiology<T>,4>& biology,
                     const std::array<std::array<T,nfishing>,4>& apical,const T& recruitment) {
  using std::exp;
  std::array<State<T>,4> z{};
  std::array<std::vector<State<T>>,4> fishing;
  for(int s=0;s<4;++s) {
    fishing[s]=fishing_at_age(apical[s],biology[s]);
    for(int j=0;j<32;++j) {
      z[s][j]=T(data.biology.natural_mortality[j]);for(const auto& f:fishing[s])z[s][j]+=f[j];
    }
  }
  State<T> january=zero_state<T>();
  for(int sex=0;sex<2;++sex) {
    const T share=sex==0 ? p.at("FracFemale_GP_1"):T(1)-p.at("FracFemale_GP_1");
    january[cell(sex,1)]=recruitment*share*exp(T(-0.25)*(z[2][cell(sex,0)]+z[3][cell(sex,0)]));
    for(int a=2;a<nage;++a) {
      T total=T(0);for(int s=0;s<4;++s)total+=z[s][cell(sex,a-1)]*0.25;
      january[cell(sex,a)]=january[cell(sex,a-1)]*exp(-total);
    }
    T total=T(0);for(int s=0;s<4;++s)total+=z[s][cell(sex,15)]*0.25;
    if(!(total>0))throw std::invalid_argument("No finite plus-group equilibrium");
    january[cell(sex,15)]/=T(1)-exp(-total);
  }
  Equilibrium<T> out;out.recruitment=recruitment;
  auto state=january;
  for(int s=0;s<4;++s) {
    if(s==2)settle(state,recruitment,p.at("FracFemale_GP_1"));
    out.beginning[s]=state;
    if(s==1)out.ssb=seasonal_ssb(data,state,fishing[s],biology[s],s+1);
    state=advance_season(state,data.biology.natural_mortality,fishing[s],0.25).survivors;
  }
  return out;
}
/// Virgin and initially fished equilibria plus January starting abundance after early recruitment deviations.
/// @details virgin and initial are distinct stationary populations. january starts from the initially fished equilibrium and then incorporates early recruitment deviations/bias; substituting virgin recruitment would change the starting age structure.
template<class T> struct InitialPopulation {
  Equilibrium<T> virgin,initial;
  State<T> january=zero_state<T>();
};
/// Calculate virgin and initial-fishing Beverton-Holt equilibria, then apply early deviations and bias correction to the initial age structure.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active evaluation tape.
/// @param data Validated original-input schema for the supported swordfish assessment.
/// @param p Native canonical-name parameter map; all required fixed and estimated values must be present.
/// @param deviations Historical fitted log recruitment deviations, including early age-structure effects.
/// @param biology Seasonal selection, weight and maturity factors from the same native parameters.
/// @return The typed result described above, with the module units and indexing conventions.
template<class T> InitialPopulation<T> initialize_population(const InputData& data,const ParameterValues<T>& p,
                  const Deviations<T>& deviations,const std::array<SeasonalBiology<T>,4>& biology) {
  using std::exp;
  std::array<std::array<T,nfishing>,4> rates;
  for(auto& season:rates)season.fill(T(0));
  InitialPopulation<T> out;
  const T r0=exp(p.at("SR_LN(R0)"));
  out.virgin=equilibrium(data,p,biology,rates,r0);
  for(int s=1;s<=4;++s)for(int f=1;f<=nfishing;++f) {
    const auto it=p.find("InitF_seas_"+std::to_string(s)+"_flt_"+std::to_string(f));
    if(it!=p.end())rates[s-1][f-1]=it->second;
  }
  const auto per=equilibrium(data,p,biology,rates,T(1));
  const T h=p.at("SR_BH_steep");
  // The initial fishing equilibrium must also lie on the BH stock-recruit
  // curve; using virgin R0 here would silently alter the starting age structure.
  const T r1=(T(4)*h*r0-(T(1)-h)*out.virgin.ssb/per.ssb)/(T(5)*h-T(1));
  if(!(r1>0))throw std::invalid_argument("Initial fishing has no positive BH equilibrium");
  out.initial=equilibrium(data,p,biology,rates,r1);
  out.january=out.initial.beginning[0];
  for(int sex=0;sex<2;++sex)for(int a=1;a<nage;++a) {
    const int year=data.dimensions.first_year-a;
    if(year>=int(data.recruitment_advanced[0]))
      out.january[cell(sex,a)]=adjusted_recruitment(data,p,deviations,year,out.january[cell(sex,a)]);
  }
  return out;
}
} // namespace swordfish

namespace swordfish {
/// One historical season's beginning and midpoint abundance, hybrid catch, removals, recruitment and biomass summaries.
/// @details beginning/midseason and removals share the abundance scale. ssb uses a negative sentinel outside the spawning season; recruits is the current recruitment quantity, not an independently drawn seasonal recruitment. Hybrid fishing is recomputed from that season observed catches.
template<class T> struct PopulationSeason {
  int year=0,season=0;
  State<T> beginning{},midseason{};
  HybridResult<T> hybrid;
  SeasonResult<T> removals;
  T recruits=T(0),ssb=T(-1),biomass=T(0),summary_biomass=T(0);
};
/// Historical seasonal trajectory and its virgin/initial equilibrium provenance.
/// @details initialization records equilibrium provenance; seasons is chronological year-major, quarter-minor order. population_row returns references into this owning vector, so those references must not outlive the run.
template<class T> struct PopulationRun {
  InitialPopulation<T> initialization;
  std::vector<PopulationSeason<T>> seasons;
};
/// Run all historical seasons from original catches and native parameters, recalculating recruitment and caching biology within this evaluation only.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active evaluation tape.
/// @param data Validated original-input schema for the supported swordfish assessment.
/// @param p Native canonical-name parameter map; all required fixed and estimated values must be present.
/// @param deviations Historical fitted log recruitment deviations, including early age-structure effects.
/// @return The typed result described above, with the module units and indexing conventions.
template<class T> PopulationRun<T> run_population(const InputData& data,const ParameterValues<T>& p,
                                                const Deviations<T>& deviations) {
  using std::exp;
  if(data.catches.empty())throw std::invalid_argument("Original catch observations are required");
  std::map<std::pair<int,int>,std::array<double,nfishing>> catches;
  for(const auto& row:data.catches) {
    if(row.year==-999)continue;
    catches[{row.year,row.season}][row.fleet-1]=row.value;
  }
  PopulationRun<T> out;
  // Biology varies only with the explicitly parsed selectivity blocks.
  // Cache within one evaluation so AD vertices are never reused across tapes.
  std::map<std::vector<int>,std::array<SeasonalBiology<T>,4>> biology_cache;
  auto biology_for=[&](int year)->const std::array<SeasonalBiology<T>,4>& {
    const auto key=selectivity_period(data,year);
    auto it=biology_cache.find(key);
    if(it==biology_cache.end()) {
      const auto curves=length_selectivity(data,p,year);
      std::array<SeasonalBiology<T>,4> value;
      for(int s=0;s<4;++s)value[s]=seasonal_biology(data,p,curves,s+1);
      it=biology_cache.emplace(key,std::move(value)).first;
    }
    return it->second;
  };
  out.initialization=initialize_population(data,p,deviations,biology_for(data.dimensions.first_year));
  auto state=out.initialization.january;
  const T r0=exp(p.at("SR_LN(R0)"));
  T recruits=T(0);
  for(int y=data.dimensions.first_year;y<=data.dimensions.last_year;++y) {
    const auto& biology=biology_for(y);
    for(int s=1;s<=4;++s) {
      if(s==3)settle(state,recruits,p.at("FracFemale_GP_1"));
      PopulationSeason<T> row;row.year=y;row.season=s;row.beginning=state;
      const auto& b=biology[s-1];
      row.hybrid=hybrid_catch(data,state,b,catches[{y,s}],0.25);
      row.removals=advance_season(state,data.biology.natural_mortality,row.hybrid.fishing,0.25);
      row.midseason=advance_season(state,data.biology.natural_mortality,row.hybrid.fishing,0.125).survivors;
      for(int j=0;j<32;++j) {
        row.biomass+=state[j]*b.weight_begin[j];
        if(j%nage>=1)row.summary_biomass+=state[j]*b.weight_begin[j];
      }
      if(s==2) {
        row.ssb=seasonal_ssb(data,state,row.hybrid.fishing,b,s);
        const T expected=beverton_holt(row.ssb,r0,out.initialization.virgin.ssb,p.at("SR_BH_steep"));
        recruits=adjusted_recruitment(data,p,deviations,y,expected*exp(p.at("SR_regime")));
      }
      row.recruits=recruits;
      state=s==4 ? age_on_january_first(row.removals.survivors):row.removals.survivors;
      out.seasons.push_back(std::move(row));
    }
  }
  return out;
}
} // namespace swordfish
