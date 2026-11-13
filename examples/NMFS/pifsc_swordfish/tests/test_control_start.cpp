#include "core/autodiff/model_gradient.hpp"
#include "swordfish/control_start.hpp"
#include "swordfish_control_input.hpp"
#include <iostream>
DECLARE_ADGRAPH();
/// Verify literal controls/zero deviations, phase masks and the recovery
/// objective gradient.
int main() {
  try {
    const auto data = imported_swordfish_control_data();
    const auto seed = swordfish::control_start_values(data);
    for (const auto &spec : data.parameters)
      if (seed.first.at(spec.name) != spec.columns[2])
        throw std::runtime_error("Non-control native start");
    if (seed.second.size() != 57 || seed.first.at("L_at_Amin_Mal_GP_1") != 83.2)
      throw std::runtime_error("Strict initialization definition");
    for (const auto &d : seed.second)
      if (d.second != 0)
        throw std::runtime_error("Nonzero starting deviation");
    const auto sensitivity = swordfish::control_start_values(data, true);
    for (const auto &row : seed.first)
      if (row.first != "L_at_Amin_Mal_GP_1" &&
          row.second != sensitivity.first.at(row.first))
        throw std::runtime_error("Unexpected fixed-config override");
    swordfish::SwordfishModel model(data, seed.first, seed.second);
    const auto full = model.initial_vector();
    bool empty = false;
    try {
      (void)model(full);
    } catch (const std::invalid_argument &error) {
      if (std::string(error.what()) != "Empty composition")
        throw;
      empty = true;
    }
    if (!empty || sensitivity.first.at("L_at_Amin_Mal_GP_1") != 83.4)
      throw std::runtime_error(
          "Control feasibility or declared fixed override changed");
    for (int phase = 1; phase <= 4; ++phase) {
      swordfish::ControlStartPhase active(model, full, phase);
      const size_t expected = phase == 1   ? 1
                              : phase == 2 ? 8
                              : phase == 3 ? 64
                                           : 84;
      if (active.active.size() != expected)
        throw std::runtime_error("Control phase activation");
      std::vector<double> sub;
      for (auto i : active.active)
        sub.push_back(full[i]);
      if (active.expand(sub) != full)
        throw std::runtime_error("Phase map changed inactive coordinates");
    }
    swordfish::ControlStartPhase feasible(model, full, 1, true);
    const std::vector<double> start{full[feasible.active.at(0)]};
    const auto gradient = quadra::evaluate_gradient(feasible, start);
    auto plus = start, minus = start;
    plus[0] += 1e-5;
    minus[0] -= 1e-5;
    const double finite = (feasible(plus) - feasible(minus)) / 2e-5;
    if (!std::isfinite(gradient.objective_value_m) ||
        std::abs(finite - gradient.gradient_m[0]) >
            1e-4 * std::max(1., std::abs(finite)))
      throw std::runtime_error("Catch-feasibility derivative");
    std::cout << "Control starts, zero deviations, phase masks and feasibility "
                 "gradient: PASS\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
