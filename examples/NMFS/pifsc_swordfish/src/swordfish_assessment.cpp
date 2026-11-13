/**
 * @file swordfish_assessment.cpp
 * @details
 * Public Quadra assessment entry point. Evaluate/gradient modes use versioned reference starts;
 * fit writes a converged native solution, and report reloads that solution. All files in build
 * are scratch outputs; only the isolated release runner publishes a complete accepted bundle.
 */
#include "swordfish/model.hpp"
#include "swordfish/reference_io.hpp"
#include "swordfish/report.hpp"
#include "swordfish_input.hpp"
#include "core/autodiff/model_gradient.hpp"
#include "core/optimizer.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <chrono>

DECLARE_ADGRAPH();
/// Run swordfish assessment; return zero on success and nonzero on a failed check or invalid input.
int main(int argc,char** argv){
  try{
    const std::string mode=argc>1?argv[1]:"--evaluate";
    if(mode!="--evaluate" && mode!="--gradient" && mode!="--fit" && mode!="--report")throw std::invalid_argument("Use --evaluate, --gradient, --fit [max_iterations] or --report");
    auto data=imported_swordfish_data();
    auto baseline=swordfish::read_reference_parameters(data,"build/biology_fitted_parameters.txt","build/population_fitted_deviations.txt");
    if(mode=="--report") {
      baseline=swordfish::read_fitted_parameters(data,"build/quadra_fitted_parameters.csv","build/biology_fitted_parameters.txt","build/population_fitted_deviations.txt");
    }
    swordfish::SwordfishModel model(std::move(data),std::move(baseline.first),std::move(baseline.second));
    auto x=model.initial_vector();
    std::cout<<std::setprecision(12)<<"Free parameters: "<<x.size()<<"; main recruitment deviations constrained to sum zero.\n";
    const auto start=std::chrono::steady_clock::now();
    double objective=0,gradient_norm=std::numeric_limits<double>::quiet_NaN();bool converged=false;int iterations=0;
    if(mode=="--evaluate" || mode=="--report") {
      objective=model(x);
      if(mode=="--report") {
        const auto native=model.decode(x);swordfish::write_fitted_outputs(model.data,native.first,native.second,"build");
      }
    }
    else if(mode=="--gradient") {
      std::cout<<"Evaluating full Quadra gradient..."<<std::endl;
      const auto g=quadra::evaluate_gradient(model,x);objective=g.objective_value_m;gradient_norm=g.gradient_norm_m;
      std::ofstream out("build/objective_gradient.csv");out<<"parameter,unconstrained_value,gradient\n"<<std::setprecision(15);
      for(size_t i=0;i<x.size();++i)out<<model.free[i].name<<','<<x[i]<<','<<g.gradient_m[i]<<'\n';
      if(!out)throw std::runtime_error("Cannot write objective gradient");
    } else {
      auto parameters=model.initial_parameters();const int maximum=argc>2?swordfish::csv_integer(argv[2]):100;
      if(maximum<1)throw std::invalid_argument("Positive iteration limit required");
      std::cout<<"Optimizing joint likelihood from SS3 parameter values..."<<std::endl;
      const auto fit=quadra::optimize_lbfgs(model,parameters,quadra::default_laplace_options(),maximum,1e-4);
      if(!fit.converged||!std::isfinite(fit.value)||!std::isfinite(fit.grad_norm))throw std::runtime_error("Fit did not converge; existing fitted outputs were not replaced");
      x=fit.par;objective=fit.value;gradient_norm=fit.grad_norm;converged=fit.converged;iterations=fit.iterations;
      auto native=model.decode(x);
      const auto final_likelihood=model.likelihood(x);
      for(const auto& row:final_likelihood.indices)native.first.at("LnQ_fleet_"+std::to_string(row.fleet))=row.log_q;
      std::ofstream out("build/quadra_fitted_parameters.csv");out<<"parameter,value\n"<<std::setprecision(15);
      for(const auto& p:native.first)out<<p.first<<','<<p.second<<'\n';
      for(const auto& d:native.second)out<<"RecruitmentDev_"<<d.first<<','<<d.second<<'\n';
      std::ofstream summary("build/quadra_fit_summary.json");summary<<std::setprecision(15)<<"{\n  \"objective\": "<<objective
        <<",\n  \"gradient_norm\": "<<gradient_norm<<",\n  \"iterations\": "<<iterations<<",\n  \"converged\": "<<(converged?"true":"false")
        <<",\n  \"free_parameters\": "<<x.size()<<",\n  \"initialization\": \"supplied SS3 fitted parameters\",\n  \"target\": \"joint penalized likelihood; no Laplace integration\"\n}\n";
      if(!out||!summary)throw std::runtime_error("Cannot write fit output");
      std::ofstream trace("build/quadra_fit_trace.csv");trace<<"evaluation,objective,gradient_norm\n"<<std::setprecision(15);
      for(const auto& r:fit.trajectory)trace<<r.evaluation<<','<<r.objective<<','<<r.gradient_norm<<'\n';
      if(!trace)throw std::runtime_error("Cannot write fit trace");
      swordfish::write_fitted_outputs(model.data,native.first,native.second,"build");
    }
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::cout<<"Objective: "<<objective;
    if(mode=="--evaluate" || mode=="--report")std::cout<<"; gradient not evaluated";
    else std::cout<<"; gradient norm: "<<gradient_norm;
    std::cout<<"; seconds: "<<seconds<<std::endl;
    return mode=="--fit"&&!converged?1:0;
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}
}
