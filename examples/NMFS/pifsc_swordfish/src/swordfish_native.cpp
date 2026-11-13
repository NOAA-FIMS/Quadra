/** @file swordfish_native.cpp
 * Reusable native model backend. Receives a validated typed stream on stdin;
 * neither generated data factories nor SS3/reference files are linked or read.
 */
#include "swordfish/native_io.hpp"
#include "swordfish/csv.hpp"
#include "swordfish/model.hpp"
#include "swordfish/report.hpp"
#include "swordfish/reference_uncertainty.hpp"
#include "core/autodiff/model_gradient.hpp"
#include "core/optimizer.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
DECLARE_ADGRAPH();
/// Evaluate or fit a native model state and write its complete derived outputs.
int main(int argc,char **argv){try{
  if(argc!=5&&argc!=6)throw std::invalid_argument("Backend usage: OUTPUT evaluate|fit MAX_ITERATIONS GRADIENT_TOLERANCE; input on stdin");
  const std::filesystem::path output=argv[1];const std::string mode=argv[2];
  const bool uncertainty=argc==6;
  if(uncertainty&&(std::string(argv[5])!="uncertainty"||mode!="fit"))throw std::invalid_argument("Uncertainty requires fit mode");
  if(mode!="fit"&&mode!="evaluate")throw std::invalid_argument("Native mode");
  const int iterations=swordfish::csv_integer(argv[3]);const double tolerance=swordfish::csv_number(argv[4]);
  if(iterations<1||!(tolerance>0))throw std::invalid_argument("Optimizer limits");
  if(std::filesystem::exists(output))throw std::invalid_argument("Output directory already exists");
  auto native=swordfish::read_native_input(std::cin);
  auto parameters=swordfish::initial_values<double>(native.first);
  swordfish::SwordfishModel model(std::move(native.first),std::move(parameters),std::move(native.second));
  auto x=model.initial_vector();const auto initial=quadra::evaluate_gradient(model,x);
  if(!std::isfinite(initial.objective_value_m)||!std::isfinite(initial.gradient_norm_m))throw std::invalid_argument("Nonfinite initial objective/gradient");
  std::filesystem::create_directories(output);
  std::cout<<std::setprecision(15)<<"Native model: "<<x.size()<<" independent coordinates; initial objective "<<initial.objective_value_m<<std::endl;
  bool converged=false;int completed=0;
  if(mode=="fit"){
    auto initial_parameters=model.initial_parameters();
    const auto fit=quadra::optimize_lbfgs(model,initial_parameters,quadra::default_laplace_options(),iterations,tolerance);
    std::ofstream trace(output/"quadra_fit_trace.csv");trace<<"evaluation,objective,gradient_norm\n"<<std::setprecision(17);
    for(const auto &row:fit.trajectory)trace<<row.evaluation<<','<<row.objective<<','<<row.gradient_norm<<'\n';
    trace.flush();if(!trace)throw std::runtime_error("Cannot save native fit trace");
    if(!fit.converged)throw std::runtime_error("Native fit did not converge; fitted outputs were not written");
    x=fit.par;completed=fit.iterations;converged=true;
  }
  const auto checked=quadra::evaluate_gradient(model,x);
  if(!std::isfinite(checked.objective_value_m)||!std::isfinite(checked.gradient_norm_m)||(mode=="fit"&&checked.gradient_norm_m>tolerance))throw std::runtime_error("Native final gradient check failed");
  auto fitted=model.decode(x);const auto likelihood=model.likelihood(x);
  for(const auto &row:likelihood.indices)fitted.first.at("LnQ_fleet_"+std::to_string(row.fleet))=row.log_q;
  swordfish::write_fitted_outputs(model.data,fitted.first,fitted.second,output);
  std::ofstream values(output/"quadra_fitted_parameters.csv");values<<"parameter,value\n"<<std::setprecision(17);
  for(const auto &row:fitted.first)values<<row.first<<','<<row.second<<'\n';
  for(const auto &row:fitted.second)values<<"RecruitmentDev_"<<row.first<<','<<row.second<<'\n';
  std::ofstream gradient(output/"objective_gradient.csv");gradient<<"parameter,unconstrained_value,gradient\n"<<std::setprecision(17);
  for(size_t i=0;i<x.size();++i)gradient<<model.free[i].name<<','<<x[i]<<','<<checked.gradient_m[i]<<'\n';
  if(uncertainty)swordfish::write_reference_point_uncertainty(model,x,output,std::cout);
  std::ofstream summary(output/"quadra_fit_summary.json");summary<<std::setprecision(17)
    <<"{\"mode\":\""<<mode<<"\",\"initialization\":\"native normalized JSON\",\"objective\":"<<checked.objective_value_m
    <<",\"initial_objective\":"<<initial.objective_value_m<<",\"gradient_norm\":"<<checked.gradient_norm_m
    <<",\"converged\":"<<(converged?"true":"false")<<",\"free_parameters\":"<<x.size()<<",\"optimizer_iterations\":"<<completed<<",\"reference_points_calculated\":true,\"reference_uncertainty_calculated\":"<<(uncertainty?"true":"false")<<"}\n";
  values.flush();gradient.flush();summary.flush();if(!values||!gradient||!summary)throw std::runtime_error("Native output write failure");
  std::cout<<"Objective: "<<checked.objective_value_m<<"; gradient norm: "<<checked.gradient_norm_m<<"; reference points written."<<std::endl;
  return 0;
}catch(const std::exception &error){std::cerr<<"Native model error: "<<error.what()<<'\n';return 1;}}
