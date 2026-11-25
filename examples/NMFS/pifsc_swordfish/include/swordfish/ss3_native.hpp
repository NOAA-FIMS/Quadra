/** @file ss3_native.hpp Strict NP swordfish SS3 3.30.20 intake in C++.
 * Reads original observations/control and optional .par starts, never
 * Report.sso.
 */
#pragma once
#include "native_files.hpp"
#include "native_json.hpp"
#include <iostream>
#include <locale>
#include <regex>

namespace swordfish {
inline double ss_number(const std::string &token) {
  std::istringstream in(token);
  in.imbue(std::locale::classic());
  double value;
  if (!(in >> value) || !std::isfinite(value) ||
      in.peek() != std::char_traits<char>::eof())
    throw std::invalid_argument("Invalid/nonfinite SS3 number: " + token);
  return value;
}
inline int ss_int(double v) {
  if (!std::isfinite(v) || v < INT_MIN || v > INT_MAX || std::trunc(v) != v)
    throw std::invalid_argument("Expected integral SS3 value");
  return int(v);
}
inline Json ss_ints(const Json &values) {
  Json result = Json::array();
  for (const auto &v : values)
    result.push_back(ss_int(json_double(v)));
  return result;
}
class SSReader {
  struct Row {
    std::vector<std::string> tokens;
    int line;
  };
  std::vector<Row> rows;
  size_t position = 0;
  std::string label;

public:
  SSReader(std::string text, std::string source) : label(std::move(source)) {
    if (text.compare(0, 3, "\xef\xbb\xbf") == 0)
      text.erase(0, 3);
    for (size_t pos = 0; (pos = text.find("\\r\\n", pos)) != std::string::npos;)
      text.replace(pos, 4, "\n");
    if (text.rfind("#V3.30.20.00;", 0) != 0)
      throw std::invalid_argument(label + ": unsupported version");
    std::istringstream in(text);
    std::string line;
    int n = 0;
    while (std::getline(in, line)) {
      ++n;
      line = line.substr(0, line.find('#'));
      std::istringstream words(line);
      Row row{{}, n};
      std::string token;
      while (words >> token)
        row.tokens.push_back(token);
      if (!row.tokens.empty())
        rows.push_back(std::move(row));
    }
  }
  std::vector<std::string> take(int width = -1) {
    if (position == rows.size())
      throw std::invalid_argument(label + ": unexpected end of input");
    const auto &row = rows[position++];
    if (width >= 0 && row.tokens.size() != size_t(width))
      throw std::invalid_argument(label + ":" + std::to_string(row.line) +
                                  ": unexpected row width");
    return row.tokens;
  }
  Json numbers(int width = -1) {
    Json result = Json::array();
    for (const auto &token : take(width))
      result.push_back(ss_number(token));
    return result;
  }
  double scalar() { return numbers(1)[0].get<double>(); }
  int integer() { return ss_int(scalar()); }
  void expect(const Json &expected) {
    if (numbers(int(expected.size())) != expected)
      throw std::invalid_argument(label +
                                  ": unsupported setting at active row " +
                                  std::to_string(position));
  }
  Json table(int width) {
    Json result = Json::array();
    for (;;) {
      auto row = numbers(width);
      if (row[0] == -9999)
        return result;
      result.push_back(row);
      if (result.size() > 1000000)
        throw std::invalid_argument("SS3 table limit");
    }
  }
  void finish() {
    expect({999});
    if (position != rows.size())
      throw std::invalid_argument(label + ": unread input after end marker");
  }
};
inline Json ss_record(const std::vector<std::string> &keys,
                      const Json &values) {
  if (keys.size() != values.size())
    throw std::logic_error("SS3 mapping width");
  Json r = Json::object();
  for (size_t i = 0; i < keys.size(); ++i)
    r[keys[i]] = values[i];
  return r;
}
inline Json read_ss_data(const std::string &text, const std::string &label) {
  SSReader r(text, label);
  Json m = Json::object(), dim = Json::object();
  dim["first_year"] = r.integer();
  dim["last_year"] = r.integer();
  dim["seasons"] = r.integer();
  if (dim["seasons"] != 4)
    throw std::invalid_argument("Supported model has four seasons");
  dim["season_months"] = r.numbers(4);
  dim["subseasons"] = r.integer();
  dim["spawn_month"] = r.scalar();
  int sexes = r.integer();
  dim["plus_age"] = r.integer();
  dim["areas"] = r.integer();
  int nfleets = r.integer();
  if (sexes != 2 || dim["areas"] != 1 || nfleets != 27 ||
      dim["plus_age"] != 15 ||
      json_int(dim["first_year"]) > json_int(dim["last_year"]) ||
      dim["season_months"] != Json::array({3, 3, 3, 3}))
    throw std::invalid_argument("Unsupported SS3 dimensions");
  m["dimensions"] = dim;
  m["fleets"] = Json::array();
  for (int f = 1; f <= 27; ++f) {
    auto row = r.take(6);
    Json v = {f};
    for (int i = 0; i < 5; ++i)
      v.push_back(ss_int(ss_number(row[i])));
    v.push_back(row[5]);
    m["fleets"].push_back(ss_record({"id", "type", "timing", "area",
                                     "catch_units", "catch_multiplier", "name"},
                                    v));
  }
  m["catches"] = Json::array();
  for (auto row : r.table(5)) {
    for (int i = 0; i < 3; ++i)
      row[i] = ss_int(json_double(row[i]));
    m["catches"].push_back(
        ss_record({"year", "season", "fleet", "value", "log_se"}, row));
  }
  m["index_specs"] = Json::array();
  for (int f = 1; f <= 27; ++f) {
    auto row = ss_ints(r.numbers(4));
    if (row[0] != f)
      throw std::invalid_argument("Index definitions must follow fleet order");
    m["index_specs"].push_back(row);
  }
  m["indices"] = Json::array();
  for (auto row : r.table(5)) {
    row[0] = ss_int(json_double(row[0]));
    row[2] = ss_int(json_double(row[2]));
    row.push_back(json_int(row[0]) > 0);
    m["indices"].push_back(
        ss_record({"year", "month", "fleet", "value", "log_se", "use"}, row));
  }
  r.expect({0});
  r.expect({0});
  r.expect({2});
  double width = r.scalar(), lo = r.scalar(), hi = r.scalar();
  if (!(width > 0 && hi > lo) || (hi - lo) / width > 1000000 ||
      std::trunc((hi - lo) / width) != (hi - lo) / width)
    throw std::invalid_argument("Invalid population length grid");
  m["population_length_bins"] = Json::array();
  for (int i = 0; i <= ss_int((hi - lo) / width); ++i)
    m["population_length_bins"].push_back(lo + i * width);
  r.expect({1});
  m["length_controls"] = Json::array();
  for (int i = 0; i < 27; ++i)
    m["length_controls"].push_back(r.numbers(7));
  int nbins = r.integer();
  if (nbins != 51)
    throw std::invalid_argument("Expected 51 length bins");
  m["length_bins"] = r.numbers(nbins);
  Json expected = Json::array();
  for (int i = 0; i < 51; ++i)
    expected.push_back(10 + 5 * i);
  if (m["length_bins"] != expected || m["population_length_bins"] != expected)
    throw std::invalid_argument("Unsupported SS3 length grid");
  m["compositions"] = Json::array();
  for (const auto &row : r.table(108)) {
    Json values = Json::array();
    for (int i = 6; i < 108; ++i)
      values.push_back(row[i]);
    m["compositions"].push_back(Json{{"year", ss_int(json_double(row[0]))},
                                     {"month", row[1]},
                                     {"fleet", ss_int(json_double(row[2]))},
                                     {"sex", ss_int(json_double(row[3]))},
                                     {"partition", ss_int(json_double(row[4]))},
                                     {"sample_size", row[5]},
                                     {"use", json_double(row[0]) > 0},
                                     {"values", values}});
  }
  for (int i = 0; i < 7; ++i)
    r.expect({0});
  r.finish();
  return m;
}
inline Json ss_parameter(SSReader &r, const std::string &name,
                         bool extended = true) {
  auto col = r.numbers(extended ? 14 : 7);
  if (!(json_double(col[0]) <= json_double(col[2]) &&
        json_double(col[2]) <= json_double(col[1])))
    throw std::invalid_argument("Parameter initial outside bounds: " + name);
  return Json{{"name", name}, {"columns", col}};
}
inline const std::vector<std::string> &ss_bio_names() {
  static const std::vector<std::string> names = {"L_at_Amin_Fem_GP_1",
                                                 "L_at_Amax_Fem_GP_1",
                                                 "VonBert_K_Fem_GP_1",
                                                 "CV_young_Fem_GP_1",
                                                 "CV_old_Fem_GP_1",
                                                 "Wtlen_1_Fem_GP_1",
                                                 "Wtlen_2_Fem_GP_1",
                                                 "Mat50%_Fem_GP_1",
                                                 "Mat_slope_Fem_GP_1",
                                                 "Eggs/kg_inter_Fem_GP_1",
                                                 "Eggs/kg_slope_wt_Fem_GP_1",
                                                 "L_at_Amin_Mal_GP_1",
                                                 "L_at_Amax_Mal_GP_1",
                                                 "VonBert_K_Mal_GP_1",
                                                 "CV_young_Mal_GP_1",
                                                 "CV_old_Mal_GP_1",
                                                 "Wtlen_1_Mal_GP_1",
                                                 "Wtlen_2_Mal_GP_1",
                                                 "CohortGrowDev",
                                                 "FracFemale_GP_1"};
  return names;
}
struct SSControl {
  Json model, parameters, settlement, ages;
  int mature_age;
};
inline SSControl read_ss_control(const std::string &text,
                                 const std::string &label, Json model) {
  SSReader r(text, label);
  for (int v : {0, 1, 1, 4, 1, 1, 0})
    r.expect({v});
  auto settlement = r.numbers(4);
  if (settlement != Json::array({1, 7, 1, 0}))
    throw std::invalid_argument("Unsupported settlement assignment");
  int nblocks = r.integer();
  if (nblocks < 0 || nblocks > 1000000)
    throw std::invalid_argument("Block count");
  auto counts = ss_ints(r.numbers(nblocks));
  model["blocks"] = Json::array();
  for (const auto &n : counts) {
    int count = json_int(n);
    if (count < 0 || count > 500000)
      throw std::invalid_argument("Block width");
    model["blocks"].push_back(r.numbers(2 * count));
  }
  r.expect({1});
  r.expect({0, 0, 0, 0, 1});
  r.expect({3});
  Json mortality = r.numbers(16);
  for (const auto &x : r.numbers(16))
    mortality.push_back(x);
  r.expect({1});
  auto a1 = r.scalar(), a2 = r.scalar(), plus = r.scalar();
  r.expect({0});
  r.expect({0});
  r.expect({0});
  r.expect({1});
  int mature = r.integer();
  r.expect({1});
  r.expect({0});
  r.expect({1});
  Json parameters = Json::array();
  for (const auto &name : ss_bio_names())
    parameters.push_back(ss_parameter(r, name));
  r.expect({0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
  r.expect({3});
  r.expect({1});
  r.expect({1});
  for (const auto &name :
       {"SR_LN(R0)", "SR_BH_steep", "SR_sigmaR", "SR_regime", "SR_autocorr"})
    parameters.push_back(ss_parameter(r, name));
  r.expect({1});
  Json rec = {{"main_start", r.integer()},
              {"main_end", r.integer()},
              {"phase", r.integer()}};
  r.expect({1});
  for (size_t i = 3; i < recruitment_fields.size(); ++i) {
    const auto &f = recruitment_fields[i];
    double v = r.scalar();
    if (f.type == 'i')
      rec[f.name] = ss_int(v);
    else
      rec[f.name] = v;
  }
  (void)r.scalar();
  (void)r.integer();
  r.expect({3});
  model["f_method"] = 3;
  model["max_f"] = r.scalar();
  model["hybrid_iterations"] = r.integer();
  for (const auto &c : model.at("catches"))
    if (c.at("year") == -999 && json_double(c.at("value")) > 0)
      parameters.push_back(ss_parameter(
          r,
          "InitF_seas_" + std::to_string(json_int(c.at("season"))) + "_flt_" +
              std::to_string(json_int(c.at("fleet"))),
          false));
  model["q_setup"] = r.table(6);
  for (const auto &q : model.at("q_setup"))
    parameters.push_back(ss_parameter(
        r, "LnQ_fleet_" + std::to_string(ss_int(json_double(q[0])))));
  model["size_patterns"] = Json::array();
  model["age_patterns"] = Json::array();
  for (int f = 0; f < 27; ++f)
    model["size_patterns"].push_back(ss_ints(r.numbers(4)));
  for (int f = 0; f < 27; ++f) {
    auto row = ss_ints(r.numbers(4));
    if (row != Json::array({0, 0, 0, 0}))
      throw std::invalid_argument("Unsupported age selectivity");
    model["age_patterns"].push_back(row);
  }
  Json size_params = Json::array();
  for (int f = 1; f <= 27; ++f) {
    const auto &row = model.at("size_patterns")[f - 1];
    int pattern = json_int(row[0]), special = json_int(row[3]);
    if (row[1] != 0 || row[2] != 0 ||
        (pattern != 0 && pattern != 1 && pattern != 15 && pattern != 24) ||
        (pattern == 15 && !(1 <= special && special < f)))
      throw std::invalid_argument("Unsupported size selectivity");
    int n = pattern == 1 ? 2 : pattern == 24 ? 6 : 0;
    for (int p = 1; p <= n; ++p) {
      auto spec = ss_parameter(r, "SizeSel_fleet_" + std::to_string(f) + "_p" +
                                      std::to_string(p));
      parameters.push_back(spec);
      size_params.push_back(spec);
    }
  }
  for (const auto &base : size_params) {
    int block = ss_int(json_double(base.at("columns")[12]));
    if (block) {
      if (block < 1 || block > nblocks)
        throw std::invalid_argument("Invalid block assignment");
      for (int b = 1; b <= json_int(counts[block - 1]); ++b)
        parameters.push_back(
            ss_parameter(r,
                         base.at("name").get<std::string>() + "_block_" +
                             std::to_string(block) + "_" + std::to_string(b),
                         false));
    }
  }
  r.expect({0});
  r.expect({0});
  model["variance_adjustments"] = r.table(3);
  model["max_lambda_phase"] = r.integer();
  model["sd_offset"] = r.integer();
  model["lambda_changes"] = r.table(5);
  r.expect({0});
  r.finish();
  model["recruitment"] = rec;
  model["plus_group_growth"] = plus;
  model["biology"] = {{"natural_mortality", mortality}};
  return {std::move(model), parameters, settlement, Json::array({a1, a2}),
          mature};
}
inline void validate_ss_observations(const Json &m) {
  int first = json_int(m.at("dimensions").at("first_year")),
      last = json_int(m.at("dimensions").at("last_year"));
  std::set<std::tuple<int, int, int>> catches;
  for (const auto &c : m.at("catches")) {
    int year = json_int(c.at("year")), season = json_int(c.at("season")),
        fleet = json_int(c.at("fleet"));
    if (fleet < 1 || fleet > 27 || m.at("fleets")[fleet - 1].at("type") != 1 ||
        season < 1 || season > 4 || json_double(c.at("value")) < 0 ||
        json_double(c.at("log_se")) < 0 ||
        (year != -999 && !(first <= year && year <= last)) ||
        !catches.emplace(year, season, fleet).second)
      throw std::invalid_argument("Invalid/duplicate SS3 catch observation");
  }
  for (const auto &kind : {"indices", "compositions"}) {
    std::set<std::tuple<double, double, int, int, int>> keys;
    for (const auto &x : m.at(kind)) {
      double year = std::abs(double(json_int(x.at("year")))),
             month = json_double(x.at("month"));
      int fleet = json_int(x.at("fleet"));
      if (!(first <= year && year <= last && 1 <= month && month < 13) ||
          fleet < 1 || fleet > 27)
        throw std::invalid_argument("Invalid SS3 observation year/month/fleet");
      int sex = 0, part = 0;
      if (std::string(kind) == "indices") {
        if (json_double(x.at("value")) <= 0 || json_double(x.at("log_se")) <= 0)
          throw std::invalid_argument("Invalid SS3 lognormal index");
      } else {
        sex = json_int(x.at("sex"));
        part = json_int(x.at("partition"));
        double sum = 0;
        for (const auto &v : x.at("values")) {
          double value = json_double(v);
          if (value < 0)
            throw std::invalid_argument("Negative composition count");
          sum += value;
        }
        if (sex < 0 || sex > 3 || part < 0 || part > 2 ||
            json_double(x.at("sample_size")) <= 0 || sum <= 0)
          throw std::invalid_argument("Invalid SS3 composition");
      }
      if (!keys.emplace(year, month, fleet, sex, part).second)
        throw std::invalid_argument("Duplicate SS3 observation");
    }
  }
}
struct SSValues {
  std::map<std::string, double> values;
  Deviations<double> deviations;
  Json metadata = nullptr;
};
inline SSValues read_ss_parameters(const std::string &text,
                                   const SSControl &control) {
  std::istringstream in(text);
  std::string line;
  std::regex header(
      R"(# Number of parameters = (\d+) Objective function value = (\S+)\s+Maximum gradient component = (\S+)\s*)");
  std::smatch match;
  if (!std::getline(in, line) || !std::regex_match(line, match, header))
    throw std::invalid_argument("Unsupported SS3 parameter header");
  if (match[1] != "86")
    throw std::invalid_argument("Unsupported SS3 parameter count");
  SSValues result;
  result.metadata = {{"ss3_parameter_count", 86},
                     {"ss3_objective", ss_number(match[2])},
                     {"ss3_maximum_gradient", ss_number(match[3])}};
  std::map<std::string, std::vector<double>> groups;
  std::string current;
  std::regex group(R"(# ([^:]+):\s*)");
  while (std::getline(in, line)) {
    if (line.find_first_not_of(" \t\r") == std::string::npos)
      continue;
    if (line.rfind("# ", 0) == 0) {
      if (!std::regex_match(line, match, group) || groups.count(match[1]))
        throw std::invalid_argument("Invalid/duplicate SS3 parameter group");
      current = match[1];
      groups.emplace(current, std::vector<double>{});
    } else {
      if (current.empty())
        throw std::invalid_argument("Values outside SS3 parameter group");
      std::istringstream words(line);
      std::string word;
      while (words >> word)
        groups[current].push_back(ss_number(word));
    }
  }
  const auto &specs = control.parameters;
  std::map<std::string, std::vector<Json>> families;
  if (specs.size() != 74)
    throw std::invalid_argument("Unsupported control parameter layout");
  for (size_t i = 0; i < specs.size(); ++i) {
    const auto &p = specs[i];
    std::string name = p.at("name");
    std::string family = i < 20                           ? "MGparm"
                         : i < 25                         ? "SR_parm"
                         : name.rfind("InitF_", 0) == 0   ? "init_F"
                         : name.rfind("LnQ_", 0) == 0     ? "Q_parm"
                         : name.rfind("SizeSel_", 0) == 0 ? "selparm"
                                                          : "unknown";
    families[family].push_back(p);
  }
  std::map<std::string, size_t> expected = {{"recdev_early", 20},
                                            {"recdev1", 37},
                                            {"Fcast_recruitments", 1},
                                            {"dummy_parm", 1},
                                            {"checksum999", 1}};
  const std::map<std::string, size_t> sizes = {{"MGparm", 20},
                                               {"SR_parm", 5},
                                               {"init_F", 1},
                                               {"Q_parm", 8},
                                               {"selparm", 40}};
  if (families.size() != sizes.size())
    throw std::invalid_argument("Unsupported control parameter families");
  for (const auto &family : sizes) {
    if (families[family.first].size() != family.second)
      throw std::invalid_argument("Unsupported control parameter family size");
    for (size_t i = 0; i < family.second; ++i) {
      std::string key = family.first + "[" + std::to_string(i + 1) + "]";
      expected[key] = 1;
      if (!groups.count(key) || groups.at(key).size() != 1)
        throw std::invalid_argument("Missing/incorrectly sized SS3 group: " +
                                    key);
      const auto &p = families[family.first][i];
      double value = groups.at(key)[0];
      const auto &col = p.at("columns");
      if (json_double(col[6]) > 0 &&
          !(json_double(col[0]) < value && value < json_double(col[1])))
        throw std::invalid_argument("Estimated parameter outside open bounds");
      result.values[p.at("name").get<std::string>()] = value;
    }
  }
  if (groups.size() != expected.size())
    throw std::invalid_argument("Unexpected/missing SS3 parameter groups");
  for (const auto &e : expected)
    if (!groups.count(e.first) || groups.at(e.first).size() != e.second)
      throw std::invalid_argument("Incorrect SS3 group size: " + e.first);
  if (groups.at("dummy_parm")[0] != 1 ||
      groups.at("Fcast_recruitments")[0] != 0)
    throw std::invalid_argument("Unsupported dummy/forecast setting");
  const auto &rec = control.model.at("recruitment");
  if (rec.at("early_start") != 1965 || rec.at("main_start") != 1985 ||
      rec.at("main_end") != 2021)
    throw std::invalid_argument("Unsupported .par recruitment mapping");
  double sum = 0;
  for (const auto &span : std::vector<std::pair<int, std::string>>{
           {1965, "recdev_early"}, {1985, "recdev1"}}) {
    int y = span.first;
    for (double value : groups.at(span.second)) {
      if (value < -5 || value > 5)
        throw std::invalid_argument("Recruitment deviation outside bounds");
      result.deviations[y++] = value;
      if (span.first == 1985)
        sum += value;
    }
  }
  if (std::abs(sum) > 1e-8)
    throw std::invalid_argument(
        "Main recruitment deviations do not sum to zero");
  result.metadata.update(Json{{"checksum999", groups.at("checksum999")[0]},
                              {"forecast_recruitment", 0.},
                              {"dummy_parameter", 1.},
                              {"main_deviation_sum", sum},
                              {"mapped_native_parameters", specs.size()},
                              {"historical_deviations", 57}});
  return result;
}
inline Json ss_document(SSControl control, const SSValues &starts) {
  auto &m = control.model;
  std::map<std::string, double> p;
  Json parameters = Json::array();
  for (const auto &spec : control.parameters) {
    const auto &col = spec.at("columns");
    std::string name = spec.at("name");
    double value =
        starts.values.empty() ? json_double(col[2]) : starts.values.at(name);
    p[name] = value;
    parameters.push_back(
        Json{{"name", name},
             {"lower", col[0]},
             {"upper", col[1]},
             {"value", value},
             {"phase", ss_int(json_double(col[6]))},
             {"prior_type", ss_int(json_double(col[5]))},
             {"block", col.size() > 12 ? ss_int(json_double(col[12])) : 0},
             {"block_function",
              col.size() > 13 ? ss_int(json_double(col[13])) : 0}});
  }
  auto &bio = m["biology"];
  bio["maturity_l50"] = p.at("Mat50%_Fem_GP_1");
  bio["maturity_slope"] = p.at("Mat_slope_Fem_GP_1");
  bio["female_fraction"] = p.at("FracFemale_GP_1");
  bio["first_mature_age"] = control.mature_age;
  bio["settlement_month"] = control.settlement[1];
  bio["settlement_age"] = control.settlement[3];
  bio["growth"] = Json::array();
  for (const auto &sex : {"Fem", "Mal"}) {
    Json g = {{"age1", control.ages[0]}, {"age2", control.ages[1]}};
    std::vector<std::string> names = {"L_at_Amin", "L_at_Amax", "VonBert_K",
                                      "CV_young",  "CV_old",    "Wtlen_1",
                                      "Wtlen_2"};
    for (size_t i = 0; i < names.size(); ++i)
      g[growth_fields[i + 2].name] = p.at(names[i] + "_" + sex + "_GP_1");
    bio["growth"].push_back(g);
  }
  Json dev = Json::array();
  const auto &rec = m.at("recruitment");
  int early = json_int(rec.at("early_start")),
      end = json_int(rec.at("main_end"));
  if (double(end) - early > 1000 || early > end)
    throw std::invalid_argument("Unsupported recruitment range");
  for (int y = early; y <= end; ++y)
    dev.push_back(Json{
        {"year", y},
        {"value", starts.deviations.empty() ? 0. : starts.deviations.at(y)}});
  Json doc = {{"schema", "quadra.swordfish.v1"},
              {"units", native_units()},
              {"model", m},
              {"parameters", parameters},
              {"recruitment_deviations", dev}};
  (void)native_from_json(doc);
  return doc;
}
inline fs::path provenance_path(fs::path output) {
  output.replace_extension(".provenance.json");
  return output;
}
inline void convert_ss3(const fs::path &data, const fs::path &control,
                        const fs::path &output,
                        const fs::path &parameters = {}) {
  auto sidecar = provenance_path(output);
  require_new_path(output);
  require_new_path(sidecar);
  if (output == sidecar)
    throw std::invalid_argument("Output path collides with provenance sidecar");
  auto data_text = read_file(data), control_text = read_file(control);
  auto model = read_ss_data(data_text, data.string());
  validate_ss_observations(model);
  auto ctl = read_ss_control(control_text, control.string(), std::move(model));
  SSValues starts;
  std::string parameter_text;
  Json sources = {
      {"data",
       {{"file", data.filename().string()}, {"sha256", sha256(data_text)}}},
      {"control",
       {{"file", control.filename().string()},
        {"sha256", sha256(control_text)}}}};
  if (!parameters.empty()) {
    parameter_text = read_file(parameters);
    starts = read_ss_parameters(parameter_text, ctl);
    sources["parameters"] = {{"file", parameters.filename().string()},
                             {"sha256", sha256(parameter_text)}};
  }
  Json document = ss_document(std::move(ctl), starts);
  std::string native_text = document.dump(2) + '\n';
  if (read_file(data) != data_text || read_file(control) != control_text ||
      (!parameters.empty() && read_file(parameters) != parameter_text))
    throw std::runtime_error("SS3 inputs changed during conversion");
  const auto &m = document.at("model");
  Json provenance = {
      {"schema", "quadra.ss3_to_pifsc_swordfish.provenance.v1"},
      {"native_schema", document.at("schema")},
      {"supported_source", "NP swordfish SS3 3.30.20 assessment layout"},
      {"sources", sources},
      {"initialization",
       parameters.empty() ? "control initial values; historical deviations zero"
                          : "supplied ss.par values"},
      {"saved_ss3_metadata", starts.metadata},
      {"report_sso_used", false},
      {"covariance_used", false},
      {"counts",
       {{"parameters", document.at("parameters").size()},
        {"historical_deviations", document.at("recruitment_deviations").size()},
        {"catches", m.at("catches").size()},
        {"indices", m.at("indices").size()},
        {"length_compositions", m.at("compositions").size()}}},
      {"native_json_sha256", sha256(native_text)}};
  std::string provenance_text = provenance.dump(2) + '\n';
  if (!output.parent_path().empty())
    fs::create_directories(output.parent_path());
  bool owns_output = false;
  try {
    exclusive_file(output, native_text);
    owns_output = true;
    exclusive_file(sidecar, provenance_text);
  } catch (...) {
    if (owns_output) {
      std::error_code ec;
      fs::remove(output, ec);
    }
    throw;
  }
  std::cout << "Native JSON: " << output.string()
            << "\nProvenance: " << sidecar.string() << '\n';
  if (parameters.empty())
    std::cout
        << "Initialization: control values and zero recruitment deviations; "
           "use the separate catch-feasibility workflow before fitting.\n";
  else
    std::cout << "Initialization: supplied .par vector (saved SS3 metadata "
                 "recorded in sidecar).\n";
}
} // namespace swordfish
