#include <cmath>
#include <iostream>
#include <type_traits>
#include <utility>
#include <vector>

#include "../core/autodiff.hpp"

DECLARE_ADGRAPH();

namespace probe {

template <typename T, typename = void>
struct has_zero_adjoints : std::false_type {};

template <typename T>
struct has_zero_adjoints<
    T, std::void_t<decltype(std::declval<T &>().zero_adjoints())>>
    : std::true_type {};

template <typename T, typename = void>
struct has_reset_adjoints : std::false_type {};

template <typename T>
struct has_reset_adjoints<
    T, std::void_t<decltype(std::declval<T &>().reset_adjoints())>>
    : std::true_type {};

template <typename T, typename = void>
struct has_clear_adjoints : std::false_type {};

template <typename T>
struct has_clear_adjoints<
    T, std::void_t<decltype(std::declval<T &>().clear_adjoints())>>
    : std::true_type {};

template <typename T, typename = void>
struct has_zero_derivatives : std::false_type {};

template <typename T>
struct has_zero_derivatives<
    T, std::void_t<decltype(std::declval<T &>().zero_derivatives())>>
    : std::true_type {};

template <typename T, typename = void>
struct has_reset_derivatives : std::false_type {};

template <typename T>
struct has_reset_derivatives<
    T, std::void_t<decltype(std::declval<T &>().reset_derivatives())>>
    : std::true_type {};

template <typename T, typename = void>
struct has_clear_derivatives : std::false_type {};

template <typename T>
struct has_clear_derivatives<
    T, std::void_t<decltype(std::declval<T &>().clear_derivatives())>>
    : std::true_type {};

template <typename T, typename = void>
struct has_graph_zero_adjoints : std::false_type {};

template <typename T>
struct has_graph_zero_adjoints<
    T, std::void_t<decltype(std::declval<T &>().graph.zero_adjoints())>>
    : std::true_type {};

template <typename T, typename = void>
struct has_graph_reset_adjoints : std::false_type {};

template <typename T>
struct has_graph_reset_adjoints<
    T, std::void_t<decltype(std::declval<T &>().graph.reset_adjoints())>>
    : std::true_type {};

template <typename T, typename = void>
struct has_graph_clear_adjoints : std::false_type {};

template <typename T>
struct has_graph_clear_adjoints<
    T, std::void_t<decltype(std::declval<T &>().graph.clear_adjoints())>>
    : std::true_type {};

} // namespace probe

int main() {
  std::cout << "had_quadra adjoint reset probe\n";

  using TapeContext = quadra::TapeContext;
  using ADScope = quadra::ADScope;

  std::cout << "ADScope has zero_adjoints(): "
            << probe::has_zero_adjoints<ADScope>::value << "\n";
  std::cout << "ADScope has reset_adjoints(): "
            << probe::has_reset_adjoints<ADScope>::value << "\n";
  std::cout << "ADScope has clear_adjoints(): "
            << probe::has_clear_adjoints<ADScope>::value << "\n";
  std::cout << "ADScope has zero_derivatives(): "
            << probe::has_zero_derivatives<ADScope>::value << "\n";
  std::cout << "ADScope has reset_derivatives(): "
            << probe::has_reset_derivatives<ADScope>::value << "\n";
  std::cout << "ADScope has clear_derivatives(): "
            << probe::has_clear_derivatives<ADScope>::value << "\n";

  std::cout << "TapeContext.graph has zero_adjoints(): "
            << probe::has_graph_zero_adjoints<TapeContext>::value << "\n";
  std::cout << "TapeContext.graph has reset_adjoints(): "
            << probe::has_graph_reset_adjoints<TapeContext>::value << "\n";
  std::cout << "TapeContext.graph has clear_adjoints(): "
            << probe::has_graph_clear_adjoints<TapeContext>::value << "\n";

  // Behavioral check: run backward twice on the same graph without mutation.
  // If the second gradient doubles, adjoints are accumulating.
  {
    quadra::TapeContext tape;
    quadra::ADScope scope(tape.graph);

    quadra::AD x = 2.0;
    quadra::AD y = x * x;

    scope.backward(y);
    Eigen::VectorXd g1 = quadra::extract_gradient(std::vector<quadra::AD>{x});

    scope.backward(y);
    Eigen::VectorXd g2 = quadra::extract_gradient(std::vector<quadra::AD>{x});

    std::cout << "repeat backward grad1 = " << g1[0] << "\n";
    std::cout << "repeat backward grad2 = " << g2[0] << "\n";

    const bool accumulates = std::abs(g2[0] - 2.0 * g1[0]) <= 1e-12;

    std::cout << "repeat backward accumulates adjoints = " << accumulates
              << "\n";
  }

  // Behavioral check with mutation from the previous test.
  {
    quadra::TapeContext tape;
    quadra::ADScope scope(tape.graph);

    quadra::AD x = 2.0;
    quadra::AD y = x * x;

    scope.backward(y);
    Eigen::VectorXd g1 = quadra::extract_gradient(std::vector<quadra::AD>{x});

    x.val = 3.0;

    scope.backward(y);
    Eigen::VectorXd g2 = quadra::extract_gradient(std::vector<quadra::AD>{x});

    std::cout << "mutation backward grad1 = " << g1[0] << "\n";
    std::cout << "mutation backward grad2 = " << g2[0] << "\n";
    std::cout << "mutation value y = " << quadra::value_of(y) << "\n";

    if (std::abs(g2[0] - 8.0) <= 1e-12) {
      std::cout << "diagnosis: second backward likely accumulated adjoints on "
                   "stale graph value\n";
    }
  }

  std::cout << "PASS: adjoint reset probe completed\n";
  return 0;
}
