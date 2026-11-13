#include "swordfish/benchmarks.hpp"
#include "swordfish/model.hpp"
#include "swordfish/reference_io.hpp"
#include "swordfish/report.hpp"
#include "swordfish_input.hpp"
#include <chrono>
#include <iostream>
/// Run test benchmarks; return zero on success and nonzero on a failed check or
/// invalid input.
int main() {
  auto data = imported_swordfish_data();
  auto p = swordfish::read_reference_parameters(
      data, "build/biology_fitted_parameters.txt",
      "build/population_fitted_deviations.txt");
  swordfish::Benchmarks b(data, p.first, p.second);
  const auto virgin = b(0), spr = b.target(.4, false),
             target = b.target(.3, true), msy = b.msy();
  if (std::abs(virgin.spr - 1) > 1e-12 || virgin.yield != 0 ||
      std::abs(virgin.annual_f) > 1e-12 || std::abs(spr.spr - .4) > 1e-12 ||
      std::abs(target.ssb / virgin.ssb - .3) > 1e-12 ||
      msy.yield < b(msy.scale * .99).yield ||
      msy.yield < b(msy.scale * 1.01).yield)
    throw std::runtime_error("Benchmark root or maximum failure");
  swordfish::SwordfishModel model(data, p.first, p.second);
  const auto points = model.reference_points(model.initial_vector());
  auto larger = p.first;
  larger.at("SR_LN(R0)") += std::log(2.);
  swordfish::SwordfishModel changed(data, larger, p.second);
  const auto changed_points =
      changed.reference_points(changed.initial_vector());
  if (points.size() != 4 || points[0].first != "unfished" ||
      points[3].first != "MSY" ||
      std::abs(changed_points[0].second.ssb / points[0].second.ssb - 2.) >
          1e-12)
    throw std::runtime_error("Model reference points do not follow parameters");
  const auto folder =
      std::filesystem::temp_directory_path() /
      ("quadra-fitted-reference-points-" +
       std::to_string(
           std::chrono::steady_clock::now().time_since_epoch().count()));
  if (!std::filesystem::create_directory(folder))
    throw std::runtime_error(
        "Cannot create reference reporting test directory");
  const auto native = model.decode(model.initial_vector());
  swordfish::write_fitted_outputs(data, native.first, native.second, folder);
  std::ifstream input(folder / "management_reference_points.csv");
  std::string line;
  if (!swordfish::csv_line(input, line) ||
      line != "benchmark,F_multiplier,SPR,SSB_tonnes,recruitment_thousands,"
              "dead_catch_tonnes,annual_F_SS3_units4_ages1_10")
    throw std::runtime_error("Missing automatic model reference-point output");
  for (const auto &point : points) {
    if (!swordfish::csv_line(input, line))
      throw std::runtime_error("Incomplete automatic reference points");
    const auto fields = swordfish::csv_fields(line);
    if (fields.size() != 7 || fields[0] != point.first ||
        std::abs(swordfish::csv_number(fields[3]) / point.second.ssb - 1.) >
            1e-12)
      throw std::runtime_error(
          "Reported reference points differ from model state");
  }
  if (swordfish::csv_line(input, line))
    throw std::runtime_error("Extra reference-point output row");
  input.close();
  std::filesystem::remove_all(folder);
  std::cout << "Equilibrium roots, MSY, parameter dependence and automatic "
               "fitted reference points pass\n";
}
