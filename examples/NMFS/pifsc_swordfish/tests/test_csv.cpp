#include "swordfish/scenarios.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace swordfish;
/// Require an invalid numeric token to raise instead of accepting a valid prefix.
template<class Function> void rejects(Function function){bool failed=false;try{function();}catch(const std::exception&){failed=true;}if(!failed)throw std::runtime_error("Malformed value accepted");}
/// Exercise strict records, finite values, trailing fields and scenario contracts.
int main(){try{
 for(const auto& text:{"1x","","1.5","2147483648"," 1"})rejects([&]{csv_integer(text);});
 for(const auto& text:{"1x","","nan","inf","1e999"," 1","1 "})rejects([&]{csv_number(text);});
 if(csv_integer("-2")!=-2||csv_number("1.25e2")!=125||csv_fields("a,b,").size()!=3)throw std::runtime_error("Valid parsing failure");
 rejects([]{csv_fields("\"a\",b");});std::istringstream windows("a,b\r\n");std::string line;if(!csv_line(windows,line)||line!="a,b")throw std::runtime_error("CRLF handling");
 const auto path=std::filesystem::temp_directory_path()/("quadra_scenario_test_"+std::to_string(std::hash<std::string>{}(std::filesystem::current_path().string()))+".csv");
 {std::ofstream out(path);out<<"scenario,name,target,value,reference_first_year,reference_last_year\n1,A,equilibrium_SPR,0.3junk,0,0\n";}rejects([&]{read_projection_scenarios(path.string());});
 std::filesystem::remove(path);std::cout<<"Strict CSV numeric, record, and scenario checks pass\n";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
