/**
 * @file reference_io.hpp
 * @brief Versioned fixture and native fitted-parameter readers.
 *
 * Canonical names, counts and complete records are required. Fitted parameter files include
 * fixed/native values and all historical deviations; they are not the 84-dimensional optimizer
 * vector. File-format and mapping failures propagate as exceptions. These readers do not authorize
 * changed source data or bypass covariance provenance.
 */
#pragma once
#include "population.hpp"
#include "csv.hpp"
#include <fstream>
#include <set>

namespace swordfish {
/// Read versioned SS3 fitted-value/deviation fixtures with exact mapping, counts and finite values; retain the report/control distinction.
/// @param data Validated original-input schema for the supported swordfish assessment.
/// @param parameter_path Versioned fitted-native-parameter fixture path.
/// @param deviation_path Versioned year/deviation fixture path.
/// @return The typed result described above, with the module units and indexing conventions.
inline std::pair<ParameterValues<double>,Deviations<double>> read_reference_parameters(
    const InputData& data,const std::string& parameter_path,const std::string& deviation_path) {
  auto p=initial_values<double>(data);Deviations<double> dev;
  std::ifstream in(parameter_path);std::string magic;size_t n;
  if(!(in>>magic>>n)||magic!="SWORDFISH_FITTED_PARAMETERS_V1"||n!=p.size())throw std::runtime_error("Parameter file header");
  std::set<std::string> seen;
  for(size_t i=0;i<n;++i){std::string name;double value;
    if(!(in>>name>>value)||!p.count(name)||!seen.insert(name).second||!std::isfinite(value))throw std::runtime_error("Parameter file mapping");
    p.at(name)=value;}
  if(in>>magic)throw std::runtime_error("Trailing parameters");
  std::ifstream din(deviation_path);const int first=int(data.recruitment_advanced[0]);
  if(!(din>>magic>>n)||magic!="SWORDFISH_DEVIATIONS_V1"||n!=size_t(data.recdev_end-first+1))throw std::runtime_error("Deviation file header");
  for(size_t i=0;i<n;++i){int year;double value;
    if(!(din>>year>>value)||year!=first+int(i)||!std::isfinite(value))throw std::runtime_error("Deviation file mapping");
    dev.emplace(year,value);}
  if(din>>magic)throw std::runtime_error("Trailing deviations");
  return {std::move(p),std::move(dev)};
}
/// Read a complete native fitted CSV using strict numeric parsing and canonical recruitment-deviation names; reject unknown, duplicate or missing rows.
/// @param data Validated original-input schema for the supported swordfish assessment.
/// @param path Input file path; records must match the documented version/schema.
/// @param reference_path Reference native-value fixture used to define the expected names.
/// @param deviations_path Reference deviation fixture defining the supported years.
/// @return The typed result described above, with the module units and indexing conventions.
inline std::pair<ParameterValues<double>,Deviations<double>> read_fitted_parameters(
    const InputData& data,const std::string& path,const std::string& reference_path,const std::string& deviations_path) {
  auto native=read_reference_parameters(data,reference_path,deviations_path);
  std::ifstream input(path);std::string line;if(!csv_line(input,line)||line!="parameter,value")throw std::runtime_error("Fitted CSV header");
  std::set<std::string> seen;
  while(csv_line(input,line)){const auto fields=csv_fields(line);if(fields.size()!=2)throw std::runtime_error("Fitted CSV row");
    const auto& name=fields[0];const double value=csv_number(fields[1]);
    if(!seen.insert(name).second)throw std::runtime_error("Duplicate fitted parameter");
    if(name.rfind("RecruitmentDev_",0)==0){const int year=csv_integer(name.substr(15));if(name!="RecruitmentDev_"+std::to_string(year))throw std::runtime_error("Noncanonical deviation name");native.second.at(year)=value;}
    else native.first.at(name)=value;
  }
  if(seen.size()!=native.first.size()+native.second.size())throw std::runtime_error("Incomplete fitted CSV");return native;
}
} // namespace swordfish
