/**
 * @file audit_biology.cpp
 * @details
 * Biological fixture audit: independently recompute growth, age-length probabilities, weights,
 * maturity and selectivity at reported native parameters. Fixtures provide comparison values and
 * parameter seeds, not replacement biological predictions.
 */
#include "swordfish/biology.hpp"
#include "swordfish_biology_input.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>

using namespace swordfish;
/// Biology audit accumulator storing counts, worst errors and the context of the worst comparison.
struct Result {
  size_t n=0; double max_error=0,max_fraction=0; std::string worst;
  /// Validate one finite biological prediction, updating absolute and tolerance-scaled errors and worst-case context.
  void check(double predicted,double expected,double atol,double rtol,const std::string& context) {
    if (!std::isfinite(predicted) || !std::isfinite(expected)) throw std::runtime_error("Nonfinite prediction: "+context);
    const double error=std::abs(predicted-expected), fraction=error/(atol+rtol*std::abs(expected));
    ++n;max_error=std::max(max_error,error);
    if (fraction>max_fraction) {max_fraction=fraction;worst=context;}
  }
};
/// Run audit biology; return zero on success and nonzero on a failed check or invalid input.
int main(int argc,char** argv) {
  try {
    if (argc!=3) throw std::runtime_error("Usage: audit_biology fixtures.txt audit.csv");
    auto data=imported_swordfish_biology_data();
    auto p=initial_values<double>(data);
    std::ifstream in(argv[1]); std::string magic;size_t n;
    if (!(in>>magic>>n) || magic!="SWORDFISH_BIOLOGY_FIXTURE_V1" || n!=p.size()) throw std::runtime_error("Biology fixture header");
    std::map<std::string,bool> seen_parameters;
    for (size_t i=0;i<n;++i) {
      std::string name; double value;
      if (!(in>>name>>value) || !std::isfinite(value) || !p.count(name) || seen_parameters[name]) throw std::runtime_error("Fitted parameter mapping");
      p.at(name)=value;seen_parameters[name]=true;
    }
    if (!(in>>n) || n==0) throw std::runtime_error("Missing biology cases");
    std::map<int,std::vector<std::vector<double>>> curves;
    std::map<std::string,Result> results;
    for (size_t i=0;i<n;++i) {
      std::string kind;int year,season,sub,sex,fleet,age,bin;double expected;
      if (!(in>>kind>>year>>season>>sub>>sex>>fleet>>age>>bin>>expected)) throw std::runtime_error("Truncated biology cases");
      if (sex<1 || sex>2 || age<0 || age>=nage || bin<0 || bin>50) throw std::runtime_error("Fixture index range");
      const double time=(season-1)*0.25;
      double predicted=0,atol=0.001,rtol=0.00003;
      if (kind=="length_begin" || kind=="length_mid" || kind=="sd_begin" || kind=="sd_mid") {
        const double t=time+(kind=="length_mid" || kind=="sd_mid" ? 0.125:0);
        const double mean=length_at_age(data,p,sex-1,age,t);
        predicted=kind.substr(0,2)=="sd" ? length_sd(data,p,sex-1,age,t,mean):mean;
      } else if (kind=="weight_begin" || kind=="weight_mid") {
        predicted=weight_at_age(data,p,sex-1,age,time+(kind=="weight_mid" ? 0.125:0));
      } else if (kind=="mature_weight") predicted=mature_weight_at_age(data,p,age);
      else if (kind=="alk") {
        predicted=age_length_key(data,p,sex-1,age,time+(sub-1)*0.125).at(bin);
        atol=0.0000003;rtol=0.00002;
      } else if (kind=="effective_parameter") {
        predicted=effective_selectivity_parameter(data,p,fleet,bin,year);atol=0.00001;rtol=0.00002;
      } else {
        if (fleet<1 || fleet>27) throw std::runtime_error("Fixture fleet index");
        if (!curves.count(year)) curves.emplace(year,length_selectivity(data,p,year));
        const auto& curve=curves.at(year).at(fleet-1);
        if (kind=="length_selectivity") {predicted=curve.at(bin);atol=0.000001;rtol=0.00004;}
        else {
          const auto selected=selected_at_age(data,p,curve,sex-1,age,time+0.125);
          if (kind=="age_selectivity") {predicted=selected.numbers;atol=0.000001;rtol=0.00005;}
          else if (kind=="selected_weight") predicted=selected.mean_weight;
          else if (kind=="selected_length") predicted=selected.mean_length;
          else throw std::runtime_error("Unknown comparison: "+kind);
        }
      }
      const std::string context=std::to_string(year)+"/"+std::to_string(season)+"/"+std::to_string(sub)+"/"+
          std::to_string(sex)+"/"+std::to_string(fleet)+"/"+std::to_string(age)+"/"+std::to_string(bin);
      results[kind].check(predicted,expected,atol,rtol,context);
    }
    std::string extra;if(in>>extra)throw std::runtime_error("Unread biology fixture data");
    std::ofstream out(argv[2]);if(!out)throw std::runtime_error("Cannot write biology audit");
    out<<"check,comparisons,max_absolute_error,max_tolerance_fraction,worst_case,passed\n"<<std::setprecision(12);
    bool passed=true;
    for (const auto& [name,r]:results) {
      bool ok=r.n>0 && r.max_fraction<=1;passed=passed&&ok;
      out<<name<<','<<r.n<<','<<r.max_error<<','<<r.max_fraction<<','<<r.worst<<','<<(ok?"true":"false")<<'\n';
      std::cout<<name<<": "<<(ok?"PASS":"FAIL")<<" ("<<r.n<<", max error "<<r.max_error<<", tolerance fraction "<<r.max_fraction<<", case "<<r.worst<<")\n";
    }
    if(!out)throw std::runtime_error("Could not write biology audit");
    return passed?0:1;
  } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}
}
