#include "core/autodiff/model_gradient.hpp"
#include "swordfish/model.hpp"
#include "swordfish/reference_io.hpp"
#include "swordfish_input.hpp"
#include <iostream>
DECLARE_ADGRAPH();
/// Run test model gradient; return zero on success and nonzero on a failed
/// check or invalid input.
int main() {
  try {
    auto data = imported_swordfish_data();
    auto reference = swordfish::read_reference_parameters(
        data, "build/biology_fitted_parameters.txt",
        "build/population_fitted_deviations.txt");
    swordfish::SwordfishModel model(std::move(data), std::move(reference.first),
                                    std::move(reference.second));
    const auto x = model.initial_vector();
    if (x.size() != 84)
      throw std::runtime_error("Independent parameter dimension");
    const auto decoded = model.decode(x);
    double sum = 0;
    for (int year = model.data.recdev_start; year <= model.data.recdev_end;
         ++year)
      sum += decoded.second.at(year);
    if (std::abs(sum) > 1e-12)
      throw std::runtime_error("Main deviation zero-sum constraint");
    const auto ad = quadra::evaluate_gradient(model, x);
    double maximum = 0;
    for (size_t j = 0; j < x.size(); ++j) {
      auto hi = x, lo = x;
      const double step = 2e-6;
      hi[j] += step;
      lo[j] -= step;
      const double fd = (model(hi) - model(lo)) / (2 * step);
      const double error = std::abs(fd - ad.gradient_m[j]);
      maximum = std::max(maximum, error);
      if (!std::isfinite(error) || error > 1e-5) {
        std::cerr << model.free[j].name << ": AD=" << ad.gradient_m[j]
                  << ", FD=" << fd << '\n';
        throw std::runtime_error("Full objective gradient mismatch");
      }
    }
    std::cout << "All 84 full-objective AD derivatives and zero-sum "
                 "constraint: PASS (maximum error "
              << maximum << ")\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
