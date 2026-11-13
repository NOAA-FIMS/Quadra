/**
 * @file benchmarks.hpp
 * @brief Equilibrium reference-point roots for a fixed fleet mixture.
 *
 * The object borrows data and native parameter maps, which must outlive it. Rates are scaled
 * proportionally across fleets/seasons; the annual F output is pooled SS3 units 4 over ages 1..10,
 * not the input multiplier. Averaging the 2019..2021 fishing pattern answers a different
 * reference-point question from using terminal-2021 rates.
 */
#pragma once
#include "likelihood.hpp"
#include <functional>
namespace swordfish {
/// Fishing multiplier, raw SPR, SSB in tonnes, recruitment in thousands, dead yield in tonnes and pooled annual F for one equilibrium.
/// @details scale multiplies the reference fleet/season pattern. spr is raw SPR, ssb/yield are tonnes, recruitment is thousands of fish and annual_f is pooled SS3 units 4. Scale and reported annual F are different quantities.
struct BenchmarkResult {double scale,spr,ssb,recruitment,yield,annual_f;};
/// Equilibrium reference-point calculator holding selected historical fleet/season fishing proportions and terminal biology.
/// @details The object borrows data and native parameter maps, which must outlive it. Rates are scaled proportionally across fleets/seasons; the annual F output is pooled SS3 units 4 over ages 1..10, not the input multiplier. Averaging the 2019..2021 fishing pattern answers a different reference-point question from using terminal-2021 rates.
class Benchmarks {
 const InputData& data;const ParameterValues<double>& p;
 std::array<SeasonalBiology<double>,4> biology;
 std::array<std::array<double,nfishing>,4> pattern{};
 double r0,h,phi0;
public:
 /// Construct terminal biology and virgin per-recruit biomass, and average the requested historical fishing-pattern years.
 /// @param d Borrowed validated input schema; must outlive the benchmark calculator.
 /// @param values Borrowed native parameter map, including fixed biology; must outlive this object.
 /// @param deviations Historical fitted log recruitment deviations, including early age-structure effects.
 /// @param first_pattern_year First historical year averaged into the fishing pattern; zero means terminal year.
 /// @param last_pattern_year Last historical year averaged into the fishing pattern; zero means terminal year.
 /// @details data and values are borrowed; do not destroy or modify them while using this calculator. Zero reference-year arguments select the terminal year.
 Benchmarks(const InputData& d,const ParameterValues<double>& values,const Deviations<double>& deviations,int first_pattern_year=0,int last_pattern_year=0):data(d),p(values) {
  const auto curves=length_selectivity(data,p,data.dimensions.last_year);
  const auto run=run_population(data,p,deviations);
  if(first_pattern_year==0)first_pattern_year=data.dimensions.last_year;
  if(last_pattern_year==0)last_pattern_year=data.dimensions.last_year;
  if(first_pattern_year<data.dimensions.first_year||last_pattern_year>data.dimensions.last_year||first_pattern_year>last_pattern_year)throw std::invalid_argument("Benchmark fishing pattern years");
  for(int s=0;s<4;++s){biology[s]=seasonal_biology(data,p,curves,s+1);for(int y=first_pattern_year;y<=last_pattern_year;++y)for(int f=0;f<nfishing;++f)pattern[s][f]+=population_row(run,data,y,s+1).hybrid.apical_f[f]/(last_pattern_year-first_pattern_year+1);}
  auto zero=pattern;for(auto& s:zero)s.fill(0);
  r0=std::exp(p.at("SR_LN(R0)"));h=p.at("SR_BH_steep");phi0=equilibrium(data,p,biology,zero,1.).ssb;
 }
 /// Solve positive Beverton-Holt recruitment, equilibrium abundance, yield and SS3 units-4 annual F for a nonnegative fishing multiplier.
 /// @param scale Nonnegative common multiplier on the reference apical-rate pattern.
 /// @return Equilibrium summaries; recruitment is clamped to zero beyond the positive BH equilibrium.
 BenchmarkResult operator()(double scale)const {
  if(scale<0)throw std::invalid_argument("Negative fishing scale");
  auto rates=pattern;for(auto& s:rates)for(auto& f:s)f*=scale;
  const auto per=equilibrium(data,p,biology,rates,1.);
  // Solve the BH equilibrium analytically after obtaining spawning biomass
  // per recruit. At collapse there is no positive equilibrium to project.
  const double recruitment=std::max(0.,(4*h*r0-(1-h)*r0*phi0/per.ssb)/(5*h-1));
  const auto eq=equilibrium(data,p,biology,rates,recruitment);
  // Yield integrates selected biomass directly: thousands of fish times kg
  // per fish is tonnes. Reported F uses pooled survivors, not mean F-at-age.
  double yield=0,natural=0,actual=0;
  for(int s=0;s<4;++s){const auto fishing=fishing_at_age(rates[s],biology[s]);
   const auto removals=advance_season(eq.beginning[s],data.biology.natural_mortality,fishing,.25);
   for(int j=0;j<32;++j){double z=data.biology.natural_mortality[j];for(const auto& f:fishing)z+=f[j];
    for(int f=0;f<nfishing;++f)yield+=eq.beginning[s][j]*rates[s][f]*biology[s].selected_biomass[f][j]*exposure(z,.25);}
   if(s==3)for(int sex=0;sex<2;++sex)for(int age=1;age<=10;++age)actual+=removals.survivors[cell(sex,age)];
  }
  for(int sex=0;sex<2;++sex)for(int age=1;age<=10;++age){const auto j=cell(sex,age);natural+=eq.beginning[0][j]*std::exp(-data.biology.natural_mortality[j]);}
  return {scale,per.ssb/phi0,eq.ssb,recruitment,yield,recruitment>0?std::log(natural/actual):0};
 }
 /// Use bracket expansion and bisection to find a valid fractional SPR or virgin-biomass target.
 /// @param value Target fraction in (0,1), interpreted by spawning_ratio.
 /// @param spawning_ratio True targets equilibrium SSB divided by unfished SSB; false targets raw SPR.
 /// @return Equilibrium summaries at the bracketed target fishing multiplier.
 BenchmarkResult target(double value,bool spawning_ratio)const {
  if(!(value>0&&value<1&&std::isfinite(value)))throw std::invalid_argument("Benchmark target fraction must be in (0,1)");
  const double b0=(*this)(0).ssb;auto quantity=[&](double f){auto row=(*this)(f);return spawning_ratio?row.ssb/b0:row.spr;};
  double lo=0,hi=1;while(quantity(hi)>value){hi*=2;if(hi>1e5)throw std::runtime_error("Benchmark target not bracketed");}
  for(int i=0;i<75;++i){const double mid=(lo+hi)/2;if(quantity(mid)>value)lo=mid;else hi=mid;}
  return (*this)((lo+hi)/2);
 }
 /// Find maximum equilibrium yield through a bracketed fourth-order numerical derivative root, avoiding flat-objective rounding noise.
 /// @return Equilibrium summaries at the yield-derivative root.
 /// @details The five-point numerical yield derivative avoids comparing almost equal maxima. Its fixed step and bracket tolerance are checked indirectly by the two-step forecast uncertainty validation; this is not a global optimizer for arbitrary selectivity models.
 BenchmarkResult msy()const {
  double lo=0,hi=1;while((*this)(hi*2).yield>(*this)(hi).yield){hi*=2;if(hi>1e5)throw std::runtime_error("MSY not bracketed");}hi*=2;
  // A derivative root avoids the flat-objective rounding noise of direct yield maximization.
  // Five-point differences are fourth order and smooth enough for delta-method derivatives.
  auto slope=[&](double f){const double step=1e-3*(1+f);if(f<2*step)return ((*this)(f+step).yield-(*this)(f).yield)/step;
   return (-(*this)(f+2*step).yield+8*(*this)(f+step).yield-8*(*this)(f-step).yield+(*this)(f-2*step).yield)/(12*step);};
  if(slope(lo)<=0||slope(hi)>=0)throw std::runtime_error("MSY derivative not bracketed");
  for(int i=0;i<60;++i){double mid=(lo+hi)/2;if(slope(mid)>0)lo=mid;else hi=mid;}
  return (*this)((lo+hi)/2);
 }
};
}
