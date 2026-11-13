#include "swordfish/biology.hpp"
#include "swordfish_biology_input.hpp"
#include "core/autodiff/model_gradient.hpp"
#include "core/model/quadra_model.hpp"
#include <iostream>

DECLARE_ADGRAPH();
using namespace swordfish;
/// Test helper for test biology: require. See the typed signature and assertions in the linked test source.
void require(bool condition,const char* message) { if(!condition)throw std::runtime_error(message); }
/// Test helper for test biology: near. See the typed signature and assertions in the linked test source.
void near(double a,double b,double tolerance,const char* message) {
  require(std::isfinite(a) && std::isfinite(b) && std::abs(a-b)<=tolerance,message);
}
/// Test helper for test biology: IntegratedPrediction. See the typed signature and assertions in the linked test source.
class IntegratedPrediction : public quadra::QuadraModel<IntegratedPrediction> {
public:
  InputData data=imported_swordfish_biology_data();
  /// Test helper for test biology: IntegratedPrediction / parameter names impl. See the typed signature and assertions in the linked test source.
  std::vector<std::string> parameter_names_impl() const {return {"fleet2_l50","female_l1"};}
  /// Test helper for test biology: IntegratedPrediction / evaluate impl. See the typed signature and assertions in the linked test source.
  template<class T> T evaluate_impl(const std::vector<T>& x,quadra::ModelReportContext&) const {
    auto p=initial_values<T>(data);
    p.at("SizeSel_fleet_2_p1")=x[0];p.at("L_at_Amin_Fem_GP_1")=x[1];
    const auto curves=length_selectivity(data,p,2000);
    const auto selected=selected_at_age(data,p,curves[1],0,3,0.375);
    // Also instantiate pattern 24 through a path contributing to the objective.
    const auto dome=selected_at_age(data,p,curves[0],0,3,0.375);
    return selected.numbers+T(0.01)*selected.biomass+T(0.001)*mature_weight_at_age(data,p,3)+T(0.1)*dome.numbers;
  }
};
/// Run test biology; return zero on success and nonzero on a failed check or invalid input.
int main() {
  try {
    auto data=imported_swordfish_biology_data();auto p=initial_values<double>(data);
    for(int i=-800;i<=800;++i) {
      const double x=i/100.;
      near(normal_cdf(x),0.5*std::erfc(-x/std::sqrt(2.)),7.5e-8,"Normal CDF approximation bound");
    }
    for(double mean:{-100.,10.,100.,1000.}) {
      const auto key=length_probabilities(data.population_length_bins,mean,10.);
      double sum=0;for(double v:key){require(v>=0 && v<=1,"Probability range");sum+=v;}
      near(sum,1,1e-14,"ALK probability mass");
    }
    near(length_probabilities(data.population_length_bins,-100.,10.).front(),1,1e-14,"Lower open tail");
    near(length_probabilities(data.population_length_bins,1000.,10.).back(),1,1e-14,"Upper open tail");
    near(maturity_at_length(p,p.at("Mat50%_Fem_GP_1")),0.5,1e-14,"Maturity midpoint");
    near(mature_weight_at_age(data,p,1),0,0,"Minimum mature age");
    near(maturity_integration_time(data),0.25,1e-14,"SS3 maturity subseason");
    for(int year:{1993,1994,1997,1998,2003,2004,2009,2010}) {
      const auto curves=length_selectivity(data,p,year);
      require(curves[4]==curves[3],"Fleet5 mirrors fleet4");
      require(curves[11]==curves[0],"Fleet12 mirrors fleet1 block");
      require(curves[24]==std::vector<double>(51,1.),"Fleet25 constant selectivity");
      const bool block=year>=1994 && year<=1997;
      near(effective_selectivity_parameter(data,p,1,1,year),
           p.at(block?"SizeSel_fleet_1_p1_block_2_1":"SizeSel_fleet_1_p1"),0,"Inclusive block endpoints");
    }
    for(int sex=0;sex<2;++sex) {
      const auto curve=std::vector<double>(51,1.);
      const auto selected=selected_at_age(data,p,curve,sex,15,0.375);
      near(selected.numbers,1,1e-14,"Unselected ALK integrates to one");
      near(selected.mean_weight,weight_at_age(data,p,sex,15,0.375),1e-10,"Selected/unselected weight identity");
      require(length_at_age(data,p,sex,15,0.375)>length_at_age(data,p,sex,14,0.375),"Plus group growth");
    }
    IntegratedPrediction model;
    std::vector<double> x={p.at("SizeSel_fleet_2_p1"),p.at("L_at_Amin_Fem_GP_1")};
    const auto ad=quadra::evaluate_gradient(model,x);
    quadra::ModelReportContext ctx;
    for(size_t j=0;j<x.size();++j) {
      const double step=0.0001;auto lo=x,hi=x;lo[j]-=step;hi[j]+=step;
      const double fd=(model.evaluate<double>(hi,ctx)-model.evaluate<double>(lo,ctx))/(2*step);
      near(ad.gradient_m[j],fd,1e-7,"Integrated biological AD derivative");
    }
    std::cout<<"Biology probability limits, weights, maturity, block boundaries, mirrors, plus group and Quadra derivatives: PASS\n";
    return 0;
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
