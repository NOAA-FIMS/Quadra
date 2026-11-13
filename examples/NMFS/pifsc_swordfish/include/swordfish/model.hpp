/**
 * @file model.hpp
 * @brief Quadra adapter and bounded independent parameter coordinates.
 *
 * This is joint penalized maximum likelihood: recruitment deviations are optimized parameters, not
 * Laplace-integrated latent variables. The logit is only an optimization reparameterization; no
 * transformation-density Jacobian is added. The dependent final main recruitment deviation is
 * reconstructed after transforming the independent coordinates.
 */
#pragma once
#include "likelihood.hpp"
#include "core/model/quadra_model.hpp"
#include "model/parameter.hpp"

namespace swordfish {
// Joint penalized-likelihood target corresponding to this SS3 setup.
// Recruitment deviations are optimized parameters here, not Laplace-integrated.
/// Quadra model adapter for the 84-coordinate joint penalized-likelihood target; historical recruitment deviations are optimized, not Laplace-integrated.
/// @details This is joint penalized maximum likelihood: recruitment deviations are optimized parameters, not Laplace-integrated latent variables. The logit is only an optimization reparameterization; no transformation-density Jacobian is added. The dependent final main recruitment deviation is reconstructed after transforming the independent coordinates.
class SwordfishModel : public quadra::QuadraModel<SwordfishModel> {
public:
  /// One independent fitting coordinate: canonical name, native bounds/value and optional recruitment year.
  /// @details lo/hi/value are native bounds and the native starting value. year=0 identifies a control parameter; a nonzero year identifies a recruitment deviation. Its position in free defines the covariance/gradient coordinate order.
  struct FreeParameter {std::string name;double lo,hi,value;int year;};
  InputData data;
  ParameterValues<double> baseline;
  Deviations<double> baseline_deviations;
  std::vector<FreeParameter> free;
  /// Build the fitting-coordinate schema from estimated control parameters and independent early/main recruitment deviations.
  /// @param input Input schema transferred into the owning model adapter.
  /// @param p Native canonical-name parameter map; all required fixed and estimated values must be present.
  /// @param deviations Historical fitted log recruitment deviations, including early age-structure effects.
  SwordfishModel(InputData input,ParameterValues<double> p,Deviations<double> deviations)
      :data(std::move(input)),baseline(std::move(p)),baseline_deviations(std::move(deviations)) {
    for(const auto& spec:data.parameters)if(spec.columns[6]>0)
      free.push_back({spec.name,spec.columns[0],spec.columns[1],baseline.at(spec.name),0});
    const double lo=data.recruitment_advanced[10],hi=data.recruitment_advanced[11];
    for(int year=int(data.recruitment_advanced[0]);year<data.recdev_end;++year) {
      const auto it=baseline_deviations.find(year);
      free.push_back({"RecruitmentDev_"+std::to_string(year),lo,hi,it==baseline_deviations.end()?0:it->second,year});
    }
  }
  /// Return canonical names in the same order as fitting coordinates and covariance matrices.
  /// @return Stable names used to label fitting coordinates, gradients and covariance rows.
  std::vector<std::string> parameter_names_impl() const {
    std::vector<std::string> names;for(const auto& x:free)names.push_back(x.name);return names;
  }
  /// Transform interior native starting values to bounded-logit coordinates and reject starts at or beyond bounds.
  /// @return Independent bounded-logit coordinates in parameter_names_impl order.
  std::vector<double> initial_vector() const {
    std::vector<double> x;
    for(const auto& f:free) {
      if(!(f.value>f.lo && f.value<f.hi))throw std::invalid_argument("Fitting initial value at/outside bounds: "+f.name);
      x.push_back(std::log((f.value-f.lo)/(f.hi-f.value)));
    }
    return x;
  }
  /// Create Quadra's parameter vector, treating every independent coordinate as a joint optimized parameter.
  /// @return Quadra ParameterVector with all independent coordinates treated as optimized joint parameters.
  quadra::ParameterVector initial_parameters() const {
    quadra::ParameterVector result;const auto x=initial_vector();
    for(size_t i=0;i<x.size();++i)result.add(quadra::Parameter(free[i].name,x[i]));return result;
  }
  /// Transform fitting coordinates to native parameters/deviations and enforce the zero-sum main-deviation constraint, including the dependent final year.
  /// @tparam T Ordinary or AD scalar; all returned AD values belong to the active evaluation tape.
  /// @param x Independent bounded-logit coordinates in free/parameter_names_impl order; length must match free.size().
  /// @return Native parameter/deviation maps; the final main deviation is dependent.
  /// @details For main years before the final year, deviations are independent bounded coordinates. The final value is minus their sum, so its box bound must be checked separately. No logit Jacobian is added to an MLE objective.
  /// @throws std::invalid_argument On a coordinate count mismatch or an infeasible dependent final deviation.
  template<class T> std::pair<ParameterValues<T>,Deviations<T>> decode(const std::vector<T>& x) const {
    if(x.size()!=free.size())throw std::invalid_argument("Swordfish fitting vector dimension");
    ParameterValues<T> p;Deviations<T> dev;
    for(const auto& item:baseline)p.emplace(item.first,T(item.second));
    for(const auto& item:baseline_deviations)dev.emplace(item.first,T(item.second));
    for(size_t i=0;i<free.size();++i) {
      const auto& spec=free[i];const T value=T(spec.lo)+T(spec.hi-spec.lo)*logistic(x[i]);
      if(spec.year)dev[spec.year]=value;else p.at(spec.name)=value;
    }
    // Main deviations live on a zero-sum hyperplane. Keeping the last year
    // dependent avoids a redundant flat direction in the Hessian/covariance.
    T main_sum=T(0);
    for(int year=data.recdev_start;year<data.recdev_end;++year)main_sum+=dev.at(year);
    dev[data.recdev_end]=-main_sum;
    if(dev.at(data.recdev_end)<data.recruitment_advanced[10] || dev.at(data.recdev_end)>data.recruitment_advanced[11])
      throw std::invalid_argument("Dependent main recruitment deviation exceeds its bounds");
    return {std::move(p),std::move(dev)};
  }
  /// Decode coordinates, run the population and return the complete component-wise likelihood result.
  /// @tparam T Ordinary or AD scalar; all returned AD values belong to the active evaluation tape.
  /// @param x Independent bounded-logit coordinates in free/parameter_names_impl order; length must match free.size().
  /// @return Population-derived components, fits and total joint objective.
  template<class T> LikelihoodResult<T> likelihood(const std::vector<T>& x) const {
    auto decoded=decode(x);
    const auto run=run_population(data,decoded.first,decoded.second);
    return evaluate_likelihood(data,decoded.first,decoded.second,run);
  }
  /// Return only the joint scalar objective for a supplied fitting-coordinate vector.
  /// @tparam T Ordinary or AD scalar; all returned AD values belong to the active evaluation tape.
  /// @param x Independent bounded-logit coordinates in free/parameter_names_impl order; length must match free.size().
  /// @return Joint penalized negative log-likelihood; no Laplace or transform-density term.
  template<class T> T operator()(const std::vector<T>& x) const {return likelihood(x).total;}
  /// Evaluate the scalar objective through Quadra's model interface and register it in the report context.
  /// @tparam T Ordinary or AD scalar; all returned AD values belong to the active evaluation tape.
  /// @param x Independent bounded-logit coordinates in free/parameter_names_impl order; length must match free.size().
  /// @param context Quadra report context receiving objective diagnostics.
  /// @return Joint objective, also registered in context.
  template<class T> T evaluate_impl(const std::vector<T>& x,quadra::ModelReportContext& context) const {
    const auto result=likelihood(x);context.report("objective",result.total);return result.total;
  }
};
} // namespace swordfish
