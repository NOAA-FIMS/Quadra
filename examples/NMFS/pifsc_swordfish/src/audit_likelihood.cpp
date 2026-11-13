/**
 * @file audit_likelihood.cpp
 * @details
 * Likelihood audit at reference native parameters: compare weighted components, profiled index
 * predictions and composition diagnostics. This isolates objective agreement before
 * optimization; assessment fitting is a separate executable mode.
 */
#include "swordfish/likelihood.hpp"
#include "swordfish_input.hpp"
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <set>

using namespace swordfish;
/// Read one typed comparison-fixture token and throw if the fixture is truncated or malformed.
template<class T> void read(std::istream& in,T& x){if(!(in>>x))throw std::runtime_error("Truncated likelihood target/input");}
/// Accumulate finite benchmark comparison counts and worst absolute/tolerance-scaled errors for this audit.
struct Check {
  size_t n=0;double error=0,fraction=0;
  /// Compare one finite model/reference pair using absolute plus relative tolerance, updating the audit summary.
  void compare(double a,double b,double atol,double rtol){
    if(!std::isfinite(a)||!std::isfinite(b))throw std::runtime_error("Nonfinite likelihood comparison");
    double e=std::abs(a-b);error=std::max(error,e);fraction=std::max(fraction,e/(atol+rtol*std::abs(b)));++n;
  }
};
/// Run audit likelihood; return zero on success and nonzero on a failed check or invalid input.
int main(int argc,char** argv){
  try{
    if(argc!=5)throw std::runtime_error("Usage: audit_likelihood parameters deviations targets output_dir");
    const auto data=imported_swordfish_data();auto p=initial_values<double>(data);
    std::ifstream pin(argv[1]);std::string magic;size_t n;read(pin,magic);read(pin,n);
    if(magic!="SWORDFISH_FITTED_PARAMETERS_V1"||n!=p.size())throw std::runtime_error("Parameter header");
    std::set<std::string> seen;
    for(size_t i=0;i<n;++i){std::string name;double v;read(pin,name);read(pin,v);
      if(!p.count(name)||!seen.insert(name).second||!std::isfinite(v))throw std::runtime_error("Parameter mapping");p.at(name)=v;}
    if(pin>>magic)throw std::runtime_error("Trailing parameters");
    std::ifstream din(argv[2]);read(din,magic);read(din,n);const int first=int(data.recruitment_advanced[0]);
    if(magic!="SWORDFISH_DEVIATIONS_V1"||n!=size_t(data.recdev_end-first+1))throw std::runtime_error("Deviation header");
    Deviations<double> dev;
    for(size_t i=0;i<n;++i){int y;double v;read(din,y);read(din,v);if(y!=first+int(i)||!std::isfinite(v))throw std::runtime_error("Deviation mapping");dev.emplace(y,v);}
    if(din>>magic)throw std::runtime_error("Trailing deviations");
    // No report-derived predictions/targets have been opened at this point.
    const auto population=run_population(data,p,dev);
    const auto fit=evaluate_likelihood(data,p,dev,population);
    const std::filesystem::path output(argv[4]);
    std::ofstream components(output/"likelihood_components.csv"),indices(output/"likelihood_index_fits.csv"),
                  lengths(output/"likelihood_length_fits.csv"),bins(output/"likelihood_length_probabilities.csv");
    if(!components||!indices||!lengths||!bins)throw std::runtime_error("Cannot write likelihood outputs");
    components<<"component,value\n"<<std::setprecision(12);
    for(const auto& [name,value]:fit.components)components<<name<<','<<value<<'\n';components<<"TOTAL,"<<fit.total<<'\n';
    indices<<"year,season,fleet,observed,expected,base,q,se,log_residual,nll,lambda\n"<<std::setprecision(12);
    for(const auto& r:fit.indices)indices<<r.year<<','<<r.season<<','<<r.fleet<<','<<r.observation<<','<<r.prediction<<','<<r.base<<','
        <<std::exp(r.log_q)<<','<<r.se<<','<<r.deviation<<','<<r.nll<<','<<likelihood_weight(data,1,r.fleet,5)<<'\n';
    lengths<<"year,season,fleet,sample_size,observed_mean,expected_mean,nll,lambda\n"<<std::setprecision(12);
    bins<<"year,season,fleet,length_midpoint,observed,expected\n"<<std::setprecision(12);
    for(const auto& r:fit.lengths){lengths<<r.year<<','<<r.season<<','<<r.fleet<<','<<r.sample_size<<','<<r.observed_mean<<','
        <<r.expected_mean<<','<<r.nll<<','<<likelihood_weight(data,4,r.fleet,5)<<'\n';
      for(size_t b=0;b<r.observed.size();++b)bins<<r.year<<','<<r.season<<','<<r.fleet<<','<<length_midpoint(data,b)<<','<<r.observed[b]<<','<<r.expected[b]<<'\n';}
    std::ifstream target(argv[3]);read(target,magic);read(target,n);
    if(magic!="SWORDFISH_LIKELIHOOD_TARGETS_V1")throw std::runtime_error("Likelihood targets header");
    std::map<std::string,Check> checks;
    for(size_t i=0;i<n;++i){std::string name;double v;read(target,name);read(target,v);
      double actual=name=="TOTAL"?fit.total:fit.components.at(name);
      checks["component_"+name].compare(actual,v,0.002,0.00001);
      std::cout<<name<<": predicted "<<actual<<", SS3 "<<v<<", difference "<<actual-v<<'\n';}
    read(target,n);if(n!=fit.indices.size())throw std::runtime_error("Index target count");
    std::map<std::tuple<int,int,int>,const IndexFit<double>*> imap;
    for(const auto& r:fit.indices)imap.emplace(std::make_tuple(r.year,r.season,r.fleet),&r);
    for(size_t i=0;i<n;++i){int y,s,f;double base,pred,q,se,dev,nll;read(target,y);read(target,s);read(target,f);
      read(target,base);read(target,pred);read(target,q);read(target,se);read(target,dev);read(target,nll);
      const auto& r=*imap.at({y,s,f});
      checks["index_base"].compare(r.base,base,0.002,0.00003);
      checks["index_prediction"].compare(r.prediction,pred,0.00001,0.00003);
      checks["index_q"].compare(std::exp(r.log_q),q,0.00000001,0.00003);
      checks["index_se"].compare(r.se,se,1e-10,1e-10);
      checks["index_log_residual"].compare(r.deviation,dev,0.00001,0.00003);
      checks["index_nll"].compare(r.nll,nll,0.0001,0.00005);
    }
    read(target,n);if(n!=fit.lengths.size())throw std::runtime_error("Length target count");
    std::map<std::tuple<int,int,int>,const LengthFit<double>*> lmap;
    for(const auto& r:fit.lengths)lmap.emplace(std::make_tuple(r.year,r.season,r.fleet),&r);
    for(size_t i=0;i<n;++i){int y,s,f;double sample,nll,om,em;read(target,y);read(target,s);read(target,f);
      read(target,sample);read(target,nll);read(target,om);read(target,em);const auto& r=*lmap.at({y,s,f});
      checks["length_sample_size"].compare(r.sample_size,sample,1e-10,1e-10);
      checks["length_nll"].compare(r.nll,nll,0.001,0.00005);
      checks["length_observed_mean"].compare(r.observed_mean,om,0.001,0.00001);
      checks["length_expected_mean"].compare(r.expected_mean,em,0.001,0.00001);
    }
    if(target>>magic)throw std::runtime_error("Trailing likelihood targets");
    std::ofstream audit(output/"likelihood_audit.csv");if(!audit)throw std::runtime_error("Cannot write likelihood audit");
    audit<<"check,comparisons,max_absolute_error,max_tolerance_fraction,passed\n"<<std::setprecision(12);bool passed=true;
    for(const auto& [name,c]:checks){bool ok=c.n>0&&c.fraction<=1;passed&=ok;
      audit<<name<<','<<c.n<<','<<c.error<<','<<c.fraction<<','<<(ok?"true":"false")<<'\n';
      std::cout<<name<<": "<<(ok?"PASS":"FAIL")<<" (max error "<<c.error<<", tolerance fraction "<<c.fraction<<")\n";}
    if(!components||!indices||!lengths||!bins||!audit)throw std::runtime_error("Failed writing likelihood output");
    return passed?0:1;
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}
}
