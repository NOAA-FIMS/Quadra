/**
 * @file csv.hpp
 * @brief Strict numeric and unquoted CSV records.
 *
 * This dialect accepts one physical record per line and normalizes CRLF endings. Quotes, embedded
 * newlines, numeric suffixes, nonfinite values and overflow are rejected. No locale-dependent
 * separator inference or spreadsheet-style coercion is performed.
 */
#pragma once
#include <charconv>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <istream>
#include <stdexcept>
#include <string>
#include <vector>
namespace swordfish {
/// Read one CSV record, accepting LF/CRLF and reporting stream exhaustion.
/// @param in Input stream; record extraction advances its position.
/// @param line One unquoted CSV record (or output record buffer for csv_line).
/// @return True for an extracted record; false only when no next line can be extracted.
inline bool csv_line(std::istream& in,std::string& line) {
  if(!std::getline(in,line))return false;
  if(!line.empty()&&line.back()=='\r')line.pop_back();return true;
}
/// Split the deliberately unquoted, single-line CSV schema; retain empty fields.
/// @param line One unquoted CSV record (or output record buffer for csv_line).
/// @return All fields in order, including a final empty field after a trailing comma.
inline std::vector<std::string> csv_fields(const std::string& line) {
  if(line.find('"')!=std::string::npos||line.find_first_of("\r\n")!=std::string::npos)throw std::invalid_argument("Quoted/multiline CSV is unsupported");
  std::vector<std::string> out;size_t start=0;
  do {const auto end=line.find(',',start);out.push_back(line.substr(start,end==std::string::npos?end:end-start));if(end==std::string::npos)break;start=end+1;}while(true);
  return out;
}
/// Parse an entire decimal integer token; reject suffixes, overflow and blanks.
/// @param text Entire numeric token; partial parses and coercion are not accepted.
/// @return Exactly parsed decimal integer.
inline int csv_integer(const std::string& text) {
  int value=0;const auto result=std::from_chars(text.data(),text.data()+text.size(),value);
  if(text.empty()||result.ec!=std::errc()||result.ptr!=text.data()+text.size())throw std::invalid_argument("Invalid CSV integer: "+text);return value;
}
/// Parse an entire finite floating-point token, rejecting range errors and whitespace.
/// @param text Entire numeric token; partial parses and coercion are not accepted.
/// @return Exactly parsed finite decimal floating-point value.
inline double csv_number(const std::string& text) {
  if(text.empty()||text.find_first_not_of("0123456789eE+-.")!=std::string::npos)throw std::invalid_argument("Invalid CSV number");
  char* end=nullptr;errno=0;const double value=std::strtod(text.c_str(),&end);
  if(end!=text.c_str()+text.size()||errno==ERANGE||!std::isfinite(value))throw std::invalid_argument("Invalid CSV number: "+text);return value;
}
}
