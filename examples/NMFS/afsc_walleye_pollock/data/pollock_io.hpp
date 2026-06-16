#pragma once

#include "pollock_data.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace pollock_example {

std::vector<std::string> split(const std::string &line)
{
  std::vector<std::string> out;
  std::stringstream ss(line);
  std::string x;
  while (std::getline(ss, x, ','))
    out.push_back(x);
  return out;
}

std::vector<Obs> read_obs(const std::string &path)
{
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("could not open " + path);
  std::string line;
  std::getline(in, line);
  std::vector<Obs> rows;
  while (std::getline(in, line))
  {
    if (line.empty())
      continue;
    auto f = split(line);
    Obs o{std::stoi(f[0]), std::stod(f[1]), std::stod(f[2]), {}};
    for (std::size_t i = 3; i < f.size(); ++i)
      o.age.push_back(std::stod(f[i]));
    rows.push_back(o);
  }
  return rows;
}

}  // namespace pollock_example

// Compatibility aliases for current walleye_pollock.cpp call sites.
using pollock_example::read_obs;
using pollock_example::split;
