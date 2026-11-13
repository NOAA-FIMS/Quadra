/**
 * @file audit_population.cpp
 * @details
 * Independent population audit: initialize and run original catches with fitted native
 * parameters/deviations, then compare computed history to SS3. Do not replace hybrid catch
 * solutions with reported F-at-age in this path.
 */
#include "swordfish/population.hpp"
#include "swordfish_population_input.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>

using namespace swordfish;
/// Accumulate finite benchmark comparison counts and worst absolute/tolerance-scaled errors for this audit.
struct Check {
  size_t n=0;double max_error=0,max_fraction=0;std::string worst;
  /// Compare one finite model/reference pair using absolute plus relative tolerance, updating the audit summary.
  void compare(double a,double b,double atol,double rtol,const std::string& where) {
    if(!std::isfinite(a)||!std::isfinite(b))throw std::runtime_error("Nonfinite comparison");
    const double error=std::abs(a-b),fraction=error/(atol+rtol*std::abs(b));
    ++n;max_error=std::max(max_error,error);if(fraction>max_fraction){max_fraction=fraction;worst=where;}
  }
};
/// Read one typed comparison-fixture token and throw if the fixture is truncated or malformed.
template<class T> void read(std::istream& in,T& value){if(!(in>>value))throw std::runtime_error("Truncated population fixture");}
/// Run audit population; return zero on success and nonzero on a failed check or invalid input.
int main(int argc,char** argv) {
  try {
    if(argc!=7)throw std::runtime_error("Usage: audit_population parameters deviations initial_targets replay audit.csv trajectory.csv");
    const auto data=imported_swordfish_population_data();
    auto p=initial_values<double>(data);
    std::ifstream pin(argv[1]);std::string magic;size_t n;
    read(pin,magic);read(pin,n);
    if(magic!="SWORDFISH_FITTED_PARAMETERS_V1" || n!=p.size())throw std::runtime_error("Fitted parameters header");
    std::set<std::string> names;
    for(size_t i=0;i<n;++i){std::string name;double value;read(pin,name);read(pin,value);
      if(!p.count(name)||!names.insert(name).second||!std::isfinite(value))throw std::runtime_error("Parameter mapping");p.at(name)=value;}
    if(pin>>magic)throw std::runtime_error("Trailing parameter values");
    std::ifstream din(argv[2]);read(din,magic);read(din,n);
    const int early=int(data.recruitment_advanced[0]);
    if(magic!="SWORDFISH_DEVIATIONS_V1" || n!=size_t(data.recdev_end-early+1))throw std::runtime_error("Deviation header");
    Deviations<double> dev;
    for(size_t i=0;i<n;++i){int year;double value;read(din,year);read(din,value);
      if(year!=early+int(i)||!std::isfinite(value))throw std::runtime_error("Deviation mapping");dev.emplace(year,value);}
    if(din>>magic)throw std::runtime_error("Trailing deviations");
    // Compute the complete trajectory before opening any report-derived targets.
    const auto run=run_population(data,p,dev);
    std::map<std::string,Check> checks;
    std::ifstream init(argv[3]);read(init,magic);
    if(magic!="SWORDFISH_INITIAL_TARGETS_V1")throw std::runtime_error("Initialization targets header");
    for(const auto& item:std::vector<std::pair<std::string,double>>{
      {"virgin_recruitment",run.initialization.virgin.recruitment},{"virgin_ssb",run.initialization.virgin.ssb},
      {"initial_recruitment",run.initialization.initial.recruitment},{"initial_ssb",run.initialization.initial.ssb}}) {
      double target;read(init,target);checks[item.first].compare(item.second,target,0.1,0.00002,"equilibrium");
    }
    const auto& i=run.initialization;
    for(const auto& item:std::vector<std::pair<std::string,State<double>>>{
      {"virgin_numbers",i.virgin.beginning[0]},{"initial_equilibrium_numbers",i.initial.beginning[0]},
      {"initial_historical_numbers",i.january}})
      for(int j=0;j<32;++j){double target;read(init,target);checks[item.first].compare(item.second[j],target,0.001,0.00002,std::to_string(j));}
    if(init>>magic)throw std::runtime_error("Trailing initialization targets");
    std::ifstream ref(argv[4]);read(ref,magic);int first,last;read(ref,first);read(ref,last);
    if(magic!="SWORDFISH_CONDITIONAL_REPLAY_V2" || first!=data.dimensions.first_year || last!=data.dimensions.last_year)
      throw std::runtime_error("Replay targets header");
    double ignored;for(int a=0;a<16;++a)read(ref,ignored);
    std::ofstream trajectory(argv[6]);if(!trajectory)throw std::runtime_error("Cannot write trajectory");
    trajectory<<"year,season,total_biomass,age1plus_biomass,spawning_biomass,recruitment_at_settlement,max_relative_catch_error\n"<<std::setprecision(12);
    std::map<std::pair<int,int>,std::array<double,nfishing>> observed;
    for(const auto& c:data.catches)if(c.year!=-999)observed[{c.year,c.season}][c.fleet-1]=c.value;
    for(const auto& row:run.seasons) {
      int year,season;read(ref,year);read(ref,season);
      if(year!=row.year || season!=row.season)throw std::runtime_error("Replay order");
      std::string where=std::to_string(year)+"/"+std::to_string(season);
      for(int sex=0;sex<2;++sex) {
        for(int a=0;a<16;++a){double target;read(ref,target);checks["begin_numbers"].compare(row.beginning[cell(sex,a)],target,0.001,0.00002,where);}
        for(int a=0;a<16;++a){double target;read(ref,target);checks["midseason_numbers"].compare(row.midseason[cell(sex,a)],target,0.001,0.00002,where);}
        for(int a=0;a<16;++a)read(ref,ignored);
        for(int f=0;f<nfishing;++f) {
          for(int a=0;a<16;++a){double target;read(ref,target);checks["f_at_age"].compare(row.hybrid.fishing[f][cell(sex,a)],target,0.000001,0.00005,where);}
          for(int a=0;a<16;++a){double target;read(ref,target);checks["catch_numbers"].compare(row.removals.catch_numbers[f][cell(sex,a)],target,0.0001,0.00003,where);}
        }
      }
      double recruits,bio,summary,ssb;read(ref,recruits);read(ref,bio);read(ref,summary);read(ref,ssb);
      if(season==3)checks["recruitment"].compare(row.recruits,recruits,0.002,0.00002,where);
      checks["total_biomass"].compare(row.biomass,bio,0.7,0.00002,where);
      checks["age1plus_biomass"].compare(row.summary_biomass,summary,0.7,0.00002,where);
      if(ssb>=0)checks["spawning_biomass"].compare(row.ssb,ssb,0.5,0.00002,where);
      double max_relative=0;
      for(int f=0;f<nfishing;++f) {
        double target;read(ref,target);checks["apical_f"].compare(row.hybrid.apical_f[f],target,0.000001,0.00005,where);
        const double obs=observed[{year,season}][f];
        if(obs>0)max_relative=std::max(max_relative,std::abs(row.hybrid.predicted_catch[f]-obs)/obs);
      }
      for(int j=0;j<32;++j) {
        double balance=row.removals.survivors[j]+row.removals.natural_deaths[j];
        for(int f=0;f<nfishing;++f)balance+=row.removals.catch_numbers[f][j];
        checks["mass_balance"].compare(balance,row.beginning[j],1e-10,1e-12,where);
      }
      trajectory<<year<<','<<season<<','<<row.biomass<<','<<row.summary_biomass<<','<<row.ssb<<','
                <<(season==3?row.recruits:0)<<','<<max_relative<<'\n';
    }
    // Consume and validate the remaining growth-only portion of the shared fixture.
    int nb;read(ref,nb);if(nb!=128)throw std::runtime_error("Biology target count");
    for(int k=0;k<nb*5;++k)read(ref,ignored);
    int ng;read(ref,ng);if(ng!=2)throw std::runtime_error("Growth target count");
    for(int k=0;k<ng*10;++k)read(ref,ignored);
    if(ref>>magic)throw std::runtime_error("Trailing population targets");
    std::ofstream out(argv[5]);if(!out)throw std::runtime_error("Cannot write population audit");
    out<<"check,comparisons,max_absolute_error,max_tolerance_fraction,worst_case,passed\n"<<std::setprecision(12);
    bool passed=true;
    for(const auto& [name,c]:checks){const bool ok=c.n>0&&c.max_fraction<=1;passed=passed&&ok;
      out<<name<<','<<c.n<<','<<c.max_error<<','<<c.max_fraction<<','<<c.worst<<','<<(ok?"true":"false")<<'\n';
      std::cout<<name<<": "<<(ok?"PASS":"FAIL")<<" ("<<c.n<<", max error "<<c.max_error<<", tolerance fraction "<<c.max_fraction<<")\n";}
    if(!out||!trajectory)throw std::runtime_error("Could not write population output");
    return passed?0:1;
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}
}
