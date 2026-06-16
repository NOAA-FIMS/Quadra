#pragma once

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace pollock_example {

struct PollockDataRow {
  int year = 0;
  double index = 0.0;
  double catch_obs = 0.0;
};

struct PollockData {
  std::vector<PollockDataRow> rows;
};

inline PollockData load_pollock_synthetic_data(const std::string &path) {
  PollockData data;
  std::ifstream in(path);
  if (!in) throw std::runtime_error("Could not open Pollock data file: " + path);

  std::string line;
  std::getline(in, line);

  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::stringstream ss(line);
    std::string item;
    PollockDataRow row;

    std::getline(ss, item, ',');
    row.year = std::stoi(item);
    std::getline(ss, item, ',');
    row.index = std::stod(item);
    std::getline(ss, item, ',');
    row.catch_obs = std::stod(item);

    data.rows.push_back(row);
  }
  return data;
}

}  // namespace pollock_example
