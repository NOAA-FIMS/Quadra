#pragma once

#include "../quadra/opakapaka_model.hpp"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace opakapaka_example {

std::vector<std::string> split_csv_line_simple(const std::string &line) {
  std::vector<std::string> fields;
  std::stringstream ss(line);
  std::string item;
  while (std::getline(ss, item, ',')) {
    fields.push_back(item);
  }
  return fields;
}

bool finite_double_from_string(const std::string &x, double &out) {
  try {
    std::size_t pos = 0;
    out = std::stod(x, &pos);
    return pos > 0 && std::isfinite(out);
  } catch (...) {
    out = std::numeric_limits<double>::quiet_NaN();
    return false;
  }
}

std::vector<Observation> read_opakapaka_history_csv(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("Could not open Opakapaka CSV: " + path);
  }

  std::string line;
  if (!std::getline(in, line)) {
    throw std::runtime_error("Opakapaka CSV is empty: " + path);
  }

  const auto header = split_csv_line_simple(line);
  int year_col = -1;
  int phase_col = -1;
  int catch_col = -1;
  int index_col = -1;

  for (int i = 0; i < static_cast<int>(header.size()); ++i) {
    if (header[i] == "year")
      year_col = i;
    if (header[i] == "phase")
      phase_col = i;
    if (header[i] == "catch_mt")
      catch_col = i;
    if (header[i] == "index")
      index_col = i;
  }

  if (year_col < 0 || phase_col < 0 || catch_col < 0 || index_col < 0) {
    throw std::runtime_error(
        "Opakapaka CSV must contain year, phase, catch_mt, and index columns");
  }

  std::vector<opakapaka_example::Observation> out;

  while (std::getline(in, line)) {
    if (line.empty())
      continue;
    const auto fields = split_csv_line_simple(line);
    const int max_col =
        std::max(std::max(year_col, phase_col), std::max(catch_col, index_col));
    if (static_cast<int>(fields.size()) <= max_col)
      continue;

    if (fields[phase_col] != "history")
      continue;

    double year_d = 0.0;
    double catch_mt = 0.0;
    double index = 0.0;

    if (!finite_double_from_string(fields[year_col], year_d))
      continue;
    if (!finite_double_from_string(fields[catch_col], catch_mt))
      continue;
    if (!finite_double_from_string(fields[index_col], index))
      continue;

    opakapaka_example::Observation obs;
    obs.year = static_cast<int>(year_d);
    obs.catch_mt = catch_mt;
    obs.index = index;
    out.push_back(obs);
  }

  if (out.empty()) {
    throw std::runtime_error(
        "No usable historical rows found in Opakapaka CSV");
  }

  return out;
}

} // namespace opakapaka_example

// Compatibility aliases for the current Opakapaka driver.
using opakapaka_example::finite_double_from_string;
using opakapaka_example::read_opakapaka_history_csv;
using opakapaka_example::split_csv_line_simple;
