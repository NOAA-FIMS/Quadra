/** Portable schema and digest tests: no SS3 observations, Python or model fit.
 */
#include "swordfish/native_files.hpp"
#include "swordfish/native_json.hpp"
#include <functional>
#include <iostream>
using namespace swordfish;
namespace {
void check(bool ok) {
  if (!ok)
    throw std::runtime_error("Native JSON test failed");
}
void rejects(const std::function<void()> &f) {
  bool failed = false;
  try {
    f();
  } catch (const std::exception &) {
    failed = true;
  }
  check(failed);
}
Json record(const JsonFields &fields) {
  Json r = Json::object();
  for (const auto &f : fields) {
    Json value = f.type == 'i' ? Json(0) : Json(0.);
    if (f.width) {
      Json a = Json::array();
      for (int i = 0; i < f.width; ++i)
        a.push_back(value);
      value = a;
    }
    r[f.name] = value;
  }
  return r;
}
Json fixture() {
  Json dim = {{"first_year", 1975}, {"last_year", 2021},
              {"seasons", 4},       {"subseasons", 2},
              {"plus_age", 15},     {"areas", 1},
              {"spawn_month", 5.},  {"season_months", {3., 3., 3., 3.}}};
  Json bio = record(biology_fields), growth = record(growth_fields);
  growth.update(Json{{"age1", 2.},
                     {"age2", 15.},
                     {"length1", 80.},
                     {"length2", 180.},
                     {"k", .2},
                     {"cv_young", .1},
                     {"cv_old", .1},
                     {"weight_coefficient", 1e-5},
                     {"weight_power", 3.}});
  bio["growth"] = {growth, growth};
  bio["female_fraction"] = .5;
  Json m = record(scalar_fields);
  m["f_method"] = 3;
  m["dimensions"] = dim;
  m["biology"] = bio;
  for (const auto &name :
       {"fleets", "catches", "indices", "compositions", "blocks",
        "population_length_bins", "length_bins", "index_specs", "size_patterns",
        "age_patterns", "length_controls", "q_setup", "variance_adjustments",
        "lambda_changes"})
    m[name] = Json::array();
  for (int i = 1; i <= 27; ++i) {
    m["fleets"].push_back(
        Json{{"id", i},
             {"type", i <= 19 ? 1 : 3},
             {"timing", 1},
             {"area", 1},
             {"catch_units", 1},
             {"catch_multiplier", 0},
             {"name", u8"Fleet \"\\ 名 " + std::to_string(i)}});
    for (const auto &name : {"index_specs", "size_patterns", "age_patterns"})
      m[name].push_back({0, 0, 0, 0});
    m["length_controls"].push_back({0., 0., 0., 0., 0., 0., 0.});
  }
  for (int i = 0; i < 51; ++i) {
    m["length_bins"].push_back(i * 5. + 10.);
    m["population_length_bins"].push_back(i * 5. + 10.);
  }
  for (int i = 20; i <= 27; ++i)
    m["q_setup"].push_back({i, 1, 0, 0, 0, 1});
  Json rec = record(recruitment_fields);
  rec.update(Json{{"main_start", 1985},
                  {"main_end", 2021},
                  {"early_start", 1965},
                  {"phase", 3},
                  {"early_phase", 4},
                  {"minimum_deviation", -5.},
                  {"maximum_deviation", 5.}});
  m["recruitment"] = rec;
  Json parameters = Json::array(), dev = Json::array();
  for (const auto &name : native_parameter_names)
    parameters.push_back(Json{{"name", name},
                              {"lower", -1000.},
                              {"upper", 1000.},
                              {"value", 0.},
                              {"phase", -1},
                              {"prior_type", 0},
                              {"block", 0},
                              {"block_function", 0}});
  for (int y = 1965; y <= 2021; ++y)
    dev.push_back(Json{{"year", y}, {"value", 0.}});
  return Json{{"schema", "quadra.swordfish.v1"},
              {"units", native_units()},
              {"model", m},
              {"parameters", parameters},
              {"recruitment_deviations", dev}};
}
} // namespace
int main() {
  try {
    check(sha256("") ==
          "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
    check(sha256("abc") ==
          "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    const auto good = fixture();
    auto loaded = native_from_json(parse_json(good.dump()));
    check(loaded.first.fleets[0].name ==
          good.at("model").at("fleets")[0].at("name"));
    check(loaded.first.parameters.size() == 74 && loaded.second.size() == 57);
    for (const auto &text :
         {R"({"a":{},"b":{"x":1,"x":2}})", R"({"a":[{"x":1,"x":2}]})",
          R"({"x":NaN})", R"({"x":1e999})", R"({"x":"\ud800"})"})
      rejects([&] { (void)parse_json(text); });
    check(parse_json(R"({"a":{"x":1},"b":{"x":2}})").size() == 2);
    for (int i = 0; i < 14; ++i) {
      auto bad = good;
      switch (i) {
      case 0:
        bad["model"]["unknown"] = 1;
        break;
      case 1:
        bad["model"]["dimensions"]["first_year"] = true;
        break;
      case 2:
        bad["parameters"].erase(0);
        break;
      case 3:
        bad["parameters"][0]["phase"] = 1;
        bad["parameters"][0]["value"] = 1000.;
        break;
      case 4:
        bad["recruitment_deviations"][56]["value"] = 1.;
        break;
      case 5:
        bad["recruitment_deviations"].erase(0);
        break;
      case 6:
        bad["model"]["fleets"][0]["id"] = 2;
        break;
      case 7:
        bad["model"]["biology"]["natural_mortality"].erase(0);
        break;
      case 8:
        bad["units"]["weight"] = "lbs";
        break;
      case 9:
        bad["parameters"][0]["block"] = 99;
        break;
      case 10:
        bad["model"]["recruitment"]["early_start"] = 1965.5;
        break;
      case 11:
        bad["model"]["dimensions"]["first_year"] = 2147483648ULL;
        break;
      case 12:
        bad["parameters"][1] = bad["parameters"][0];
        break;
      case 13:
        bad["model"]["indices"].push_back(Json{{"year", 2021},
                                               {"month", 5.},
                                               {"fleet", 20},
                                               {"value", 1.},
                                               {"log_se", .1},
                                               {"use", 1}});
        break;
      }
      rejects([&] { (void)native_from_json(bad); });
    }
    std::cout << "Native JSON and SHA-256 checks passed.\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
