#include "swordfish/dynamics.hpp"
#include "core/autodiff/model_gradient.hpp"
#include "core/model/quadra_model.hpp"
#include <iostream>
#include <numeric>

DECLARE_ADGRAPH();
using namespace swordfish;
/// Test helper for test dynamics: require. See the typed signature and assertions in the linked test source.
void require(bool ok, const char* message) {
  if (!ok) throw std::runtime_error(message);
}
/// Test helper for test dynamics: near. See the typed signature and assertions in the linked test source.
void near(double value,double expected,double tolerance,const char* message) {
  require(std::isfinite(value) && std::abs(value-expected)<=tolerance,message);
}
/// Test helper for test dynamics: SurvivalModel. See the typed signature and assertions in the linked test source.
class SurvivalModel : public quadra::QuadraModel<SurvivalModel> {
public:
  /// Test helper for test dynamics: SurvivalModel / parameter names impl. See the typed signature and assertions in the linked test source.
  std::vector<std::string> parameter_names_impl() const { return {"log_f"}; }
  template <class T> T evaluate_impl(const std::vector<T>& parameters,
                                    /// Test helper for test dynamics: SurvivalModel / evaluate impl. See the typed signature and assertions in the linked test source.
                                    quadra::ModelReportContext&) const {
    using std::exp;
    State<T> numbers=zero_state<T>(); State<double> m{};
    std::vector<State<T>> f(1);f[0].fill(T(0));
    numbers[cell(0,5)]=T(100); m[cell(0,5)]=0.2;
    f[0][cell(0,5)]=exp(parameters[0]);
    auto r=advance_season(numbers,m,f,0.25);
    // Include both survivor and catch terms, exercising exposure's AD path.
    return r.survivors[cell(0,5)]+T(2)*r.catch_numbers[0][cell(0,5)];
  }
};
/// Test helper for test dynamics: ExposureModel. See the typed signature and assertions in the linked test source.
class ExposureModel : public quadra::QuadraModel<ExposureModel> {
public:
  /// Test helper for test dynamics: ExposureModel / parameter names impl. See the typed signature and assertions in the linked test source.
  std::vector<std::string> parameter_names_impl() const { return {"z"}; }
  /// Test helper for test dynamics: ExposureModel / evaluate impl. See the typed signature and assertions in the linked test source.
  template <class T> T evaluate_impl(const std::vector<T>& p, quadra::ModelReportContext&) const {
    return exposure(p[0],0.25);
  }
};
/// Run test dynamics; return zero on success and nonzero on a failed check or invalid input.
int main() {
  try {
    State<double> numbers{},m{};
    numbers[cell(0,0)]=100; numbers[cell(1,0)]=80;
    auto zero=advance_season(numbers,m,std::vector<State<double>>(2),0.25);
    near(zero.survivors[0],100,0,"Zero mortality survival");
    near(zero.catch_numbers[0][0],0,0,"Zero mortality catch");
    near(exposure(0.,0.25),0.25,0,"Zero mortality exposure");
    numbers[cell(0,14)]=3; numbers[cell(0,15)]=7;
    auto next=age_on_january_first(numbers);
    near(next[cell(0,15)],10,0,"Plus-group accumulation");
    near(next[cell(0,0)],0,0,"Age zero clears in January");
    settle(next,200.,0.6);
    near(next[cell(0,0)],120,0,"Female recruitment");
    near(next[cell(1,0)],80,0,"Male recruitment");
    near(beverton_holt(1000.,500.,1000.,0.9),500,1e-12,"BH at B0");
    near(beverton_holt(200.,500.,1000.,0.9),450,1e-12,"BH steepness definition");
    near(beverton_holt(0.,500.,1000.,0.9),0,0,"BH at zero SSB");
    State<double> spawn_numbers{}, spawn_m{};
    spawn_numbers[cell(0,3)]=100; spawn_m[cell(0,3)]=0.2;
    std::array<double,nage> mw{}; mw[3]=10;
    near(spawning_biomass(spawn_numbers,spawn_m,{},mw,1./12),1000*std::exp(-0.2/12),1e-10,"Calendar spawning survival");
    near(spawning_biomass_ss33020(spawn_numbers,spawn_m,{},mw,1./12,0.25),1000*std::exp(-0.2/3),1e-10,"SS33020 spawning survival");
    ExposureModel exposure_model;
    const auto exposure_gradient=quadra::evaluate_gradient(exposure_model,{0.0});
    near(exposure_gradient.objective_value_m,0.25,1e-15,"Exposure zero limit under AD");
    near(exposure_gradient.gradient_m[0],-0.25*0.25/2,1e-15,"Exposure derivative at zero");
    SurvivalModel model;
    for (double rate : {0.3,1e-8}) {
      auto result=quadra::evaluate_gradient(model,{std::log(rate)});
      const double z=0.2+rate, dt=0.25, surv=100*std::exp(-z*dt);
      const double catch_n=100*rate*(1-std::exp(-z*dt))/z;
      const double catch_derivative=100*rate*(0.2*(1-std::exp(-z*dt))/(z*z)+rate*dt*std::exp(-z*dt)/z);
      near(result.objective_value_m,surv+2*catch_n,1e-10,"AD objective");
      near(result.gradient_m[0],-dt*rate*surv+2*catch_derivative,1e-10,"Analytic AD derivative");
    }
    std::cout << "Dynamics limits, plus group, settlement, BH, and Quadra AD derivatives: PASS\n";
    return 0;
  } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
