#include "surplus_production.hpp"

int main() {
  const auto data = quadra_examples::surplus_production::make_demo_data();
  const auto par = quadra_examples::surplus_production::make_demo_parameters();
  quadra_examples::surplus_production::print_report(data, par);
  return 0;
}
