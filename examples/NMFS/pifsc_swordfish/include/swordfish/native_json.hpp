/** @file native_json.hpp Strict, versioned native JSON intake (no external
 * process). */
#pragma once
#include "native_io.hpp"
#include "third_party/json.hpp"
#include <climits>
#include <map>
#include <sstream>

namespace swordfish {
using Json = nlohmann::json;
inline Json parse_json(const std::string &text) {
  std::vector<std::set<std::string>> objects;
  return Json::parse(
      text, [&](int depth, Json::parse_event_t event, Json &value) {
        if (depth > 64)
          throw std::invalid_argument("JSON nesting limit");
        if (event == Json::parse_event_t::object_start)
          objects.emplace_back();
        if (event == Json::parse_event_t::key &&
            !objects.back().insert(value.get<std::string>()).second)
          throw std::invalid_argument("Duplicate JSON field: " +
                                      value.get<std::string>());
        if (event == Json::parse_event_t::object_end)
          objects.pop_back();
        return true;
      });
}
inline void json_keys(const Json &j, const std::vector<std::string> &keys) {
  if (!j.is_object() || j.size() != keys.size())
    throw std::invalid_argument("Missing/unknown JSON fields");
  for (const auto &key : keys)
    if (!j.contains(key))
      throw std::invalid_argument("Missing JSON field: " + key);
}
inline const Json &json_rows(const Json &j) {
  if (!j.is_array() || j.size() > 1000000)
    throw std::invalid_argument("Invalid JSON collection size");
  return j;
}
inline double json_double(const Json &j) {
  if (!j.is_number() || !std::isfinite(j.get<double>()))
    throw std::invalid_argument("Expected finite JSON number");
  return j.get<double>();
}
inline int json_int(const Json &j) {
  if (!j.is_number_integer() || json_double(j) < INT_MIN ||
      json_double(j) > INT_MAX)
    throw std::invalid_argument("Expected JSON integer in 32-bit range");
  return j.get<int>();
}
struct JsonField {
  std::string name;
  char type;
  int width = 0;
};
using JsonFields = std::vector<JsonField>;
inline void json_atom(std::ostream &out, const Json &j, char type,
                      int width = 0) {
  if (width) {
    if (!j.is_array() || j.size() != size_t(width))
      throw std::invalid_argument("Invalid JSON fixed array width");
    for (const auto &v : j) {
      json_atom(out, v, type);
      out << ' ';
    }
  } else if (type == 'i')
    out << json_int(j);
  else if (type == 'd')
    out << json_double(j);
  else if (type == 'b') {
    if (!j.is_boolean())
      throw std::invalid_argument("Expected JSON boolean");
    out << int(j.get<bool>());
  } else {
    if (!j.is_string())
      throw std::invalid_argument("Expected JSON string");
    auto s = j.get<std::string>();
    for (unsigned char c : s)
      if (c < 32)
        throw std::invalid_argument("Control character in native string");
    out << std::quoted(s);
  }
}
inline void json_record(std::ostream &out, const Json &j,
                        const JsonFields &fields) {
  std::vector<std::string> keys;
  for (const auto &f : fields)
    keys.push_back(f.name);
  json_keys(j, keys);
  for (const auto &f : fields) {
    json_atom(out, j.at(f.name), f.type, f.width);
    out << ' ';
  }
  out << '\n';
}
inline const JsonFields dimension_fields = {
    {"first_year", 'i', 0},  {"last_year", 'i', 0},    {"seasons", 'i', 0},
    {"subseasons", 'i', 0},  {"plus_age", 'i', 0},     {"areas", 'i', 0},
    {"spawn_month", 'd', 0}, {"season_months", 'd', 4}};
inline const JsonFields biology_fields = {
    {"natural_mortality", 'd', 32}, {"maturity_l50", 'd', 0},
    {"maturity_slope", 'd', 0},     {"female_fraction", 'd', 0},
    {"first_mature_age", 'i', 0},   {"settlement_month", 'd', 0},
    {"settlement_age", 'd', 0}};
inline const JsonFields growth_fields = {
    {"age1", 'd', 0},        {"age2", 'd', 0},
    {"length1", 'd', 0},     {"length2", 'd', 0},
    {"k", 'd', 0},           {"cv_young", 'd', 0},
    {"cv_old", 'd', 0},      {"weight_coefficient", 'd', 0},
    {"weight_power", 'd', 0}};
inline const JsonFields scalar_fields = {{"f_method", 'i', 0},
                                         {"hybrid_iterations", 'i', 0},
                                         {"max_lambda_phase", 'i', 0},
                                         {"sd_offset", 'i', 0},
                                         {"max_f", 'd', 0},
                                         {"plus_group_growth", 'd', 0}};
inline const JsonFields recruitment_fields = {{"main_start", 'i', 0},
                                              {"main_end", 'i', 0},
                                              {"phase", 'i', 0},
                                              {"early_start", 'i', 0},
                                              {"early_phase", 'i', 0},
                                              {"forecast_phase", 'i', 0},
                                              {"forecast_lambda", 'd', 0},
                                              {"last_early_no_bias", 'd', 0},
                                              {"first_full_bias", 'd', 0},
                                              {"last_full_bias", 'd', 0},
                                              {"first_recent_no_bias", 'd', 0},
                                              {"maximum_bias", 'd', 0},
                                              {"cycles", 'i', 0},
                                              {"minimum_deviation", 'd', 0},
                                              {"maximum_deviation", 'd', 0},
                                              {"read_deviations", 'i', 0}};
inline const JsonFields parameter_fields = {
    {"name", 's', 0},  {"lower", 'd', 0},         {"upper", 'd', 0},
    {"value", 'd', 0}, {"phase", 'i', 0},         {"prior_type", 'i', 0},
    {"block", 'i', 0}, {"block_function", 'i', 0}};
inline const JsonFields fleets_fields = {
    {"id", 'i', 0},   {"type", 'i', 0},        {"timing", 'i', 0},
    {"area", 'i', 0}, {"catch_units", 'i', 0}, {"catch_multiplier", 'i', 0},
    {"name", 's', 0}};
inline const JsonFields catches_fields = {{"year", 'i', 0},
                                          {"season", 'i', 0},
                                          {"fleet", 'i', 0},
                                          {"value", 'd', 0},
                                          {"log_se", 'd', 0}};
inline const JsonFields indices_fields = {{"year", 'i', 0},   {"month", 'd', 0},
                                          {"fleet", 'i', 0},  {"value", 'd', 0},
                                          {"log_se", 'd', 0}, {"use", 'b', 0}};
inline const JsonFields compositions_fields = {
    {"year", 'i', 0}, {"month", 'd', 0},     {"fleet", 'i', 0},
    {"sex", 'i', 0},  {"partition", 'i', 0}, {"sample_size", 'd', 0},
    {"use", 'b', 0},  {"values", 'd', 102}};
inline const std::set<std::string> native_parameter_names = {
    "CV_old_Fem_GP_1",
    "CV_old_Mal_GP_1",
    "CV_young_Fem_GP_1",
    "CV_young_Mal_GP_1",
    "CohortGrowDev",
    "Eggs/kg_inter_Fem_GP_1",
    "Eggs/kg_slope_wt_Fem_GP_1",
    "FracFemale_GP_1",
    "InitF_seas_1_flt_6",
    "L_at_Amax_Fem_GP_1",
    "L_at_Amax_Mal_GP_1",
    "L_at_Amin_Fem_GP_1",
    "L_at_Amin_Mal_GP_1",
    "LnQ_fleet_20",
    "LnQ_fleet_21",
    "LnQ_fleet_22",
    "LnQ_fleet_23",
    "LnQ_fleet_24",
    "LnQ_fleet_25",
    "LnQ_fleet_26",
    "LnQ_fleet_27",
    "Mat50%_Fem_GP_1",
    "Mat_slope_Fem_GP_1",
    "SR_BH_steep",
    "SR_LN(R0)",
    "SR_autocorr",
    "SR_regime",
    "SR_sigmaR",
    "SizeSel_fleet_1_p1",
    "SizeSel_fleet_1_p1_block_2_1",
    "SizeSel_fleet_1_p2",
    "SizeSel_fleet_1_p2_block_2_1",
    "SizeSel_fleet_1_p3",
    "SizeSel_fleet_1_p3_block_2_1",
    "SizeSel_fleet_1_p4",
    "SizeSel_fleet_1_p4_block_2_1",
    "SizeSel_fleet_1_p5",
    "SizeSel_fleet_1_p6",
    "SizeSel_fleet_2_p1",
    "SizeSel_fleet_2_p1_block_1_1",
    "SizeSel_fleet_2_p2",
    "SizeSel_fleet_2_p2_block_1_1",
    "SizeSel_fleet_3_p1",
    "SizeSel_fleet_3_p2",
    "SizeSel_fleet_3_p3",
    "SizeSel_fleet_3_p4",
    "SizeSel_fleet_3_p5",
    "SizeSel_fleet_3_p6",
    "SizeSel_fleet_4_p1",
    "SizeSel_fleet_4_p2",
    "SizeSel_fleet_6_p1",
    "SizeSel_fleet_6_p2",
    "SizeSel_fleet_6_p3",
    "SizeSel_fleet_6_p4",
    "SizeSel_fleet_6_p5",
    "SizeSel_fleet_6_p6",
    "SizeSel_fleet_7_p1",
    "SizeSel_fleet_7_p2",
    "SizeSel_fleet_7_p3",
    "SizeSel_fleet_7_p4",
    "SizeSel_fleet_7_p5",
    "SizeSel_fleet_7_p6",
    "SizeSel_fleet_8_p1",
    "SizeSel_fleet_8_p2",
    "SizeSel_fleet_8_p3",
    "SizeSel_fleet_8_p4",
    "SizeSel_fleet_8_p5",
    "SizeSel_fleet_8_p6",
    "VonBert_K_Fem_GP_1",
    "VonBert_K_Mal_GP_1",
    "Wtlen_1_Fem_GP_1",
    "Wtlen_1_Mal_GP_1",
    "Wtlen_2_Fem_GP_1",
    "Wtlen_2_Mal_GP_1"};
inline Json native_units() {
  return Json::parse(
      R"({"numbers": "thousands of fish", "weight": "kg", "biomass": "metric tonnes", "length": "cm EFL", "mortality": "annual instantaneous rate"})");
}
inline void validate_json_layout(const Json &document) {
  const auto &m = document.at("model"), &d = m.at("dimensions"),
             &rec = m.at("recruitment");
  if (d.at("seasons") != 4 || d.at("subseasons") != 2 ||
      d.at("plus_age") != 15 || d.at("areas") != 1 ||
      json_int(d.at("first_year")) < 1 || json_int(d.at("last_year")) > 9999 ||
      json_int(d.at("last_year")) < json_int(d.at("first_year")) ||
      json_int(d.at("last_year")) - json_int(d.at("first_year")) > 500 ||
      d.at("last_year") != rec.at("main_end"))
    throw std::invalid_argument("Unsupported model dimensions/year range");
  const auto &fleets = m.at("fleets");
  if (fleets.size() != 27)
    throw std::invalid_argument("Expected 27 fleets");
  for (int i = 0; i < 27; ++i) {
    const auto &f = fleets.at(i);
    if (f.at("id") != i + 1 || f.at("type") != (i < 19 ? 1 : 3) ||
        f.at("area") != 1 || f.at("catch_multiplier") != 0)
      throw std::invalid_argument("Unsupported fleet structure");
  }
  for (const auto &name :
       {"index_specs", "size_patterns", "age_patterns", "length_controls"})
    if (m.at(name).size() != 27)
      throw std::invalid_argument("Expected 27 fleet settings");
  for (const auto &r : m.at("age_patterns"))
    if (r != Json::array({0, 0, 0, 0}))
      throw std::invalid_argument("Unsupported age selectivity");
  if (m.at("f_method") != 3)
    throw std::invalid_argument("Unsupported catch conditioning");
  for (const auto &name : {"population_length_bins", "length_bins"}) {
    const auto &bins = m.at(name);
    if (bins.size() != 51)
      throw std::invalid_argument("Expected 51 length bins");
    for (size_t i = 1; i < bins.size(); ++i)
      if (!(json_double(bins[i - 1]) < json_double(bins[i])))
        throw std::invalid_argument("Length bins must increase");
  }
  if (m.at("q_setup").size() != 8)
    throw std::invalid_argument("Expected 8 catchability settings");
  for (int i = 0; i < 8; ++i)
    if (m.at("q_setup").at(i).at(0) != i + 20)
      throw std::invalid_argument("Catchability fleet order");
  const int first = json_int(d.at("first_year")),
            last = json_int(d.at("last_year"));
  auto year_ok = [&](const Json &r) {
    int y = json_int(r.at("year"));
    return first <= y && y <= last;
  };
  for (const auto &r : m.at("catches")) {
    int fleet = json_int(r.at("fleet")), season = json_int(r.at("season"));
    if (fleet < 1 || fleet > 19 || season < 1 || season > 4 ||
        json_double(r.at("value")) < 0 || json_double(r.at("log_se")) < 0 ||
        (r.at("year") != -999 && !year_ok(r)))
      throw std::invalid_argument("Invalid catch observation");
  }
  for (const auto &r : m.at("indices")) {
    int fleet = json_int(r.at("fleet"));
    if (fleet < 20 || fleet > 27 ||
        (r.at("use").get<bool>() &&
         (json_double(r.at("value")) <= 0 || json_double(r.at("log_se")) <= 0 ||
          !year_ok(r))))
      throw std::invalid_argument("Invalid abundance index");
  }
  for (const auto &r : m.at("compositions")) {
    int fleet = json_int(r.at("fleet"));
    double sum = 0;
    for (const auto &v : r.at("values")) {
      double x = json_double(v);
      if (x < 0)
        throw std::invalid_argument("Negative composition value");
      sum += x;
    }
    if (fleet < 1 || fleet > 27 || r.at("sex") != 0 || r.at("partition") != 0 ||
        (r.at("use").get<bool>() &&
         (json_double(r.at("sample_size")) <= 0 || sum <= 0 || !year_ok(r))))
      throw std::invalid_argument("Invalid length composition");
  }
}
inline std::pair<InputData, Deviations<double>>
native_from_json(const Json &doc) {
  json_keys(doc, {"schema", "units", "model", "parameters",
                  "recruitment_deviations"});
  if (doc.at("schema") != "quadra.swordfish.v1" ||
      doc.at("units") != native_units())
    throw std::invalid_argument("Unsupported native schema or units");
  const auto &m = doc.at("model");
  json_keys(m, {"dimensions",
                "biology",
                "recruitment",
                "blocks",
                "population_length_bins",
                "length_bins",
                "fleets",
                "catches",
                "indices",
                "compositions",
                "index_specs",
                "size_patterns",
                "age_patterns",
                "length_controls",
                "q_setup",
                "variance_adjustments",
                "lambda_changes",
                "f_method",
                "hybrid_iterations",
                "max_lambda_phase",
                "sd_offset",
                "max_f",
                "plus_group_growth"});
  std::ostringstream out;
  out << std::setprecision(17) << "QUADRA_SWORDFISH_NATIVE_V1\n";
  json_record(out, m.at("dimensions"), dimension_fields);
  auto bio = m.at("biology");
  auto keys = std::vector<std::string>{"growth"};
  for (const auto &f : biology_fields)
    keys.push_back(f.name);
  json_keys(bio, keys);
  const auto growth = bio.at("growth");
  bio.erase("growth");
  json_record(out, bio, biology_fields);
  if (json_rows(growth).size() != 2)
    throw std::invalid_argument("Expected two growth records");
  for (const auto &g : growth)
    json_record(out, g, growth_fields);
  for (const auto &table : std::vector<std::pair<std::string, JsonFields>>{
           {"fleets", fleets_fields},
           {"catches", catches_fields},
           {"indices", indices_fields},
           {"compositions", compositions_fields}}) {
    const auto &rows = json_rows(m.at(table.first));
    out << rows.size() << '\n';
    for (const auto &r : rows)
      json_record(out, r, table.second);
  }
  for (const auto &name : {"population_length_bins", "length_bins"}) {
    const auto &rows = json_rows(m.at(name));
    out << rows.size() << ' ';
    for (const auto &v : rows) {
      json_atom(out, v, 'd');
      out << ' ';
    }
    out << '\n';
  }
  for (const auto &table : std::vector<std::tuple<std::string, int, char>>{
           {"index_specs", 4, 'i'},
           {"size_patterns", 4, 'i'},
           {"age_patterns", 4, 'i'},
           {"length_controls", 7, 'd'},
           {"q_setup", 6, 'd'},
           {"variance_adjustments", 3, 'd'},
           {"lambda_changes", 5, 'd'}}) {
    const auto &rows = json_rows(m.at(std::get<0>(table)));
    out << rows.size() << '\n';
    for (const auto &r : rows) {
      json_atom(out, r, std::get<2>(table), std::get<1>(table));
      out << '\n';
    }
  }
  const auto &blocks = json_rows(m.at("blocks"));
  out << blocks.size() << '\n';
  for (const auto &r : blocks) {
    if (json_rows(r).size() % 2)
      throw std::invalid_argument("Block years must be paired");
    out << r.size() << ' ';
    for (const auto &v : r) {
      json_atom(out, v, 'd');
      out << ' ';
    }
    out << '\n';
  }
  Json scalars = Json::object();
  for (const auto &f : scalar_fields)
    scalars[f.name] = m.at(f.name);
  json_record(out, scalars, scalar_fields);
  json_record(out, m.at("recruitment"), recruitment_fields);
  const auto &specs = json_rows(doc.at("parameters"));
  out << specs.size() << '\n';
  std::set<std::string> names;
  for (const auto &p : specs) {
    json_record(out, p, parameter_fields);
    if (!names.insert(p.at("name").get<std::string>()).second)
      throw std::invalid_argument("Duplicate parameter name");
    int block = json_int(p.at("block")), fun = json_int(p.at("block_function"));
    if (block < 0 || size_t(block) > blocks.size() || fun < 0 || fun > 2)
      throw std::invalid_argument("Unsupported parameter block");
  }
  if (names != native_parameter_names)
    throw std::invalid_argument("Incomplete/unsupported parameter names");
  const auto &rec = m.at("recruitment");
  int early = json_int(rec.at("early_start")),
      start = json_int(rec.at("main_start")),
      end = json_int(rec.at("main_end"));
  if (!(early <= start && start <= end) || double(end) - early > 1000 ||
      rec.at("cycles") != 0 || rec.at("read_deviations") != 0)
    throw std::invalid_argument("Unsupported recruitment range/options");
  const auto &dev = json_rows(doc.at("recruitment_deviations"));
  out << dev.size() << '\n';
  for (const auto &r : dev)
    json_record(out, r, {{"year", 'i'}, {"value", 'd'}});
  validate_json_layout(doc);
  std::istringstream in(out.str());
  return read_native_input(in);
}
} // namespace swordfish
