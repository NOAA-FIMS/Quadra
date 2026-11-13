/**
 * @file control_start.hpp
 * @brief Control-only native initialization and phase-restricted optimization.
 * No report estimates, fitted parameter files or saved covariance are read here.
 */
#pragma once
#include "model.hpp"
namespace swordfish {
/// Build literal native control starts and zero historical recruitment deviations.
/// The optional male-growth override is an explicit fixed-configuration sensitivity,
/// not a source of fitted free parameters; strict mode leaves every control value intact.
inline std::pair<ParameterValues<double>,Deviations<double>> control_start_values(const InputData& data,bool match_fixed_growth=false) {
  auto p=initial_values<double>(data);Deviations<double> dev;
  for(int y=int(data.recruitment_advanced[0]);y<=data.recdev_end;++y)dev.emplace(y,0.);
  if(match_fixed_growth)p.at("L_at_Amin_Mal_GP_1")=83.4;
  return {std::move(p),std::move(dev)};
}
/// Restrict optimization to one control activation phase, retaining inactive values.
/// The full model owns the objective and native bounds. base is a fitting vector,
/// not a native parameter map, so selection introduces no extra transformation.
struct ControlStartPhase : quadra::QuadraModel<ControlStartPhase> {
  swordfish::SwordfishModel& model;
  std::vector<double> base;
  std::vector<size_t> active;
  bool catch_only;
  /// Borrow a model and copy the full coordinate checkpoint for an activation phase.
  ControlStartPhase(swordfish::SwordfishModel& m,std::vector<double> point,int phase,bool feasibility=false):model(m),base(std::move(point)),catch_only(feasibility) {
    for(size_t i=0;i<model.free.size();++i) {
      const auto& f=model.free[i];int activation=0;
      if(f.year)activation=f.year<model.data.recdev_start?int(model.data.recruitment_advanced[1]):model.data.recdev_phase;
      else for(const auto& p:model.data.parameters)if(p.name==f.name){activation=int(p.columns[6]);break;}
      if(activation>0&&activation<=phase)active.push_back(i);
    }
    if(active.empty())throw std::invalid_argument("No active control-start coordinates");
  }
  /// Return active names in optimizer order, preserving the original full ordering.
  std::vector<std::string> parameter_names_impl()const {std::vector<std::string> names;for(auto i:active)names.push_back(model.free[i].name);return names;}
  /// Expand active coordinates into the full fitting vector, keeping inactive values fixed.
  template<class T> std::vector<T> expand(const std::vector<T>& point)const {
    if(point.size()!=active.size())throw std::invalid_argument("Phase coordinate dimension");
    std::vector<T> full;for(double value:base)full.push_back(T(value));
    for(size_t j=0;j<active.size();++j)full[active[j]]=point[j];return full;
  }
  /// Evaluate the full final-phase likelihood through the active-coordinate map.
  template<class T>T operator()(const std::vector<T>& point)const{
    const auto full=expand(point);if(!catch_only)return model(full);
    const auto decoded=model.decode(full);const auto run=swordfish::run_population(model.data,decoded.first,decoded.second);
    T loss=T(0);for(const auto& row:model.data.catches)if(row.year!=-999&&row.value>0)
      loss+=swordfish::stabilized_catch_nll(row.value,swordfish::population_row(run,model.data,row.year,row.season).hybrid.predicted_catch[row.fleet-1],row.log_se)*swordfish::likelihood_weight(model.data,8,row.fleet,5);
    return loss;
  }
  /// Expose the phase objective through the public Quadra model interface.
  template<class T>T evaluate_impl(const std::vector<T>& point,quadra::ModelReportContext&)const{return (*this)(point);}
};

} // namespace swordfish
