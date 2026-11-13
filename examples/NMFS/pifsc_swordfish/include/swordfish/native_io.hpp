/** @file native_io.hpp
 * Read the versioned typed transport produced by the native JSON frontend.
 * No SS3 parser, generated input factory or reference fixture is required.
 */
#pragma once
#include "population.hpp"
#include <iomanip>
#include <istream>
#include <set>

namespace swordfish {
/// Read a finite numeric field, rejecting truncation and invalid stream tokens.
template <class T> void native_number(std::istream &in, T &value) {
  if (!(in >> value) || !std::isfinite(double(value)))
    throw std::invalid_argument("Invalid/truncated native numeric input");
}
/// Read a bounded collection length before allocating observation containers.
inline size_t native_count(std::istream &in) {
  long long count;
  native_number(in, count);
  if (count < 0 || count > 1000000)
    throw std::invalid_argument("Native collection size");
  return size_t(count);
}
/// Fill fixed arrays from their statically declared schema width.
template <class T, size_t N>
void native_array(std::istream &in, std::array<T, N> &a) {
  for (auto &value : a)
    native_number(in, value);
}
/// Read a length-prefixed scalar vector.
inline void native_vector(std::istream &in, std::vector<double> &v) {
  v.resize(native_count(in));
  for (auto &value : v)
    native_number(in, value);
}
/// Read a length-prefixed fixed-width settings table.
template <class T, size_t N>
void native_table(std::istream &in, std::vector<std::array<T, N>> &v) {
  v.resize(native_count(in));
  for (auto &row : v)
    native_array(in, row);
}
/// Read a native state and its initial deviations from a validated frontend
/// stream. C++ repeats finite/count/duplicate/shape checks and invokes
/// InputData validation. The public JSON format is documented by
/// native_input.py; transport order is private to this version. The initial
/// parameter map uses each native value.
inline std::pair<InputData, Deviations<double>>
read_native_input(std::istream &in) {
  std::string magic;
  if (!(in >> magic) || magic != "QUADRA_SWORDFISH_NATIVE_V1")
    throw std::invalid_argument("Native input version");
  InputData d{};
  auto &dim = d.dimensions;
  native_number(in, dim.first_year);
  native_number(in, dim.last_year);
  native_number(in, dim.seasons);
  native_number(in, dim.subseasons);
  native_number(in, dim.plus_age);
  native_number(in, dim.areas);
  native_number(in, dim.spawn_month);
  native_array(in, dim.season_months);
  auto &b = d.biology;
  native_array(in, b.natural_mortality);
  native_number(in, b.maturity_l50);
  native_number(in, b.maturity_slope);
  native_number(in, b.female_fraction);
  native_number(in, b.first_mature_age);
  native_number(in, b.settlement_month);
  native_number(in, b.settlement_age);
  for (auto &g : b.growth) {
    native_number(in, g.age1);
    native_number(in, g.age2);
    native_number(in, g.length1);
    native_number(in, g.length2);
    native_number(in, g.k);
    native_number(in, g.cv_young);
    native_number(in, g.cv_old);
    native_number(in, g.weight_coefficient);
    native_number(in, g.weight_power);
  }
  d.fleets.resize(native_count(in));
  for (auto &f : d.fleets) {
    native_number(in, f.id);
    native_number(in, f.type);
    native_number(in, f.timing);
    native_number(in, f.area);
    native_number(in, f.catch_units);
    native_number(in, f.catch_multiplier);
    if (!(in >> std::quoted(f.name)))
      throw std::invalid_argument("Native fleet name");
  }
  d.catches.resize(native_count(in));
  for (auto &r : d.catches) {
    native_number(in, r.year);
    native_number(in, r.season);
    native_number(in, r.fleet);
    native_number(in, r.value);
    native_number(in, r.log_se);
  }
  d.indices.resize(native_count(in));
  for (auto &r : d.indices) {
    native_number(in, r.year);
    native_number(in, r.month);
    native_number(in, r.fleet);
    native_number(in, r.value);
    native_number(in, r.log_se);
    int use;
    native_number(in, use);
    if (use != 0 && use != 1)
      throw std::invalid_argument("Native inclusion flag");
    r.use = use;
  }
  d.compositions.resize(native_count(in));
  for (auto &r : d.compositions) {
    native_number(in, r.year);
    native_number(in, r.month);
    native_number(in, r.fleet);
    native_number(in, r.sex);
    native_number(in, r.partition);
    native_number(in, r.sample_size);
    int use;
    native_number(in, use);
    if (use != 0 && use != 1)
      throw std::invalid_argument("Native inclusion flag");
    r.use = use;
    native_array(in, r.values);
  }
  native_vector(in, d.population_length_bins);
  native_vector(in, d.length_bins);
  native_table(in, d.index_specs);
  native_table(in, d.size_patterns);
  native_table(in, d.age_patterns);
  native_table(in, d.length_controls);
  native_table(in, d.q_setup);
  native_table(in, d.variance_adjustments);
  native_table(in, d.lambda_changes);
  d.blocks.resize(native_count(in));
  for (auto &row : d.blocks)
    native_vector(in, row);
  native_number(in, d.f_method);
  native_number(in, d.hybrid_iterations);
  native_number(in, d.max_lambda_phase);
  native_number(in, d.sd_offset);
  native_number(in, d.max_f);
  native_number(in, d.plus_group_growth);
  native_number(in, d.recdev_start);
  native_number(in, d.recdev_end);
  native_number(in, d.recdev_phase);
  native_array(in, d.recruitment_advanced);
  d.parameters.resize(native_count(in));
  std::set<std::string> names;
  for (auto &p : d.parameters) {
    if (!(in >> std::quoted(p.name)) || !names.insert(p.name).second)
      throw std::invalid_argument("Native parameter name");
    p.source_line = 0;
    p.columns.assign(14, 0.);
    native_number(in, p.columns[0]);
    native_number(in, p.columns[1]);
    native_number(in, p.columns[2]);
    int phase, prior, block, function;
    native_number(in, phase);
    native_number(in, prior);
    native_number(in, block);
    native_number(in, function);
    p.columns[6] = phase;
    p.columns[5] = prior;
    p.columns[12] = block;
    p.columns[13] = function;
    if (prior != 0 || p.columns[0] > p.columns[1] ||
        (phase > 0 &&
         !(p.columns[0] < p.columns[2] && p.columns[2] < p.columns[1])))
      throw std::invalid_argument("Native parameter bounds/prior");
  }
  Deviations<double> dev;
  const auto count = native_count(in);
  for (size_t i = 0; i < count; ++i) {
    int year;
    double value;
    native_number(in, year);
    native_number(in, value);
    if (!dev.emplace(year, value).second)
      throw std::invalid_argument("Native duplicate deviation");
  }
  if (in >> magic)
    throw std::invalid_argument("Trailing native input");
  if (d.parameters.size() != 74 || d.index_specs.size() != 27 ||
      d.size_patterns.size() != 27 || d.age_patterns.size() != 27 ||
      d.length_controls.size() != 27)
    throw std::invalid_argument("Incomplete native settings");
  if (d.recdev_end < int(d.recruitment_advanced[0]) ||
      double(d.recdev_end) - d.recruitment_advanced[0] > 1000)
    throw std::invalid_argument("Native recruitment range");
  double sum = 0.;
  for (int y = int(d.recruitment_advanced[0]); y <= d.recdev_end; ++y) {
    if (!dev.count(y) || dev.at(y) < d.recruitment_advanced[10] ||
        dev.at(y) > d.recruitment_advanced[11])
      throw std::invalid_argument("Native deviation years/bounds");
    if (y >= d.recdev_start)
      sum += dev.at(y);
  }
  if (dev.size() != size_t(d.recdev_end - int(d.recruitment_advanced[0]) + 1) ||
      std::abs(sum) > 1e-8)
    throw std::invalid_argument("Native deviation constraint");
  d.validate();
  return {std::move(d), std::move(dev)};
}
} // namespace swordfish
