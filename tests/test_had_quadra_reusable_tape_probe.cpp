#include <iostream>
#include <type_traits>
#include <utility>
#include <vector>

#include "../core/autodiff.hpp"

// This is a capability probe, not yet a reusable-tape implementation.
//
// Goal:
//
// Determine whether the current Quadra AD abstraction exposes a way to update
// independent variable values on an existing tape/graph, then run a new
// forward/reverse pass without rebuilding the graph.
//
// If this test reports that mutable independent values are not exposed, the
// next patch should add that capability in had_quadra.hpp/autodiff.hpp before
// we build reusable Laplace evaluators.

namespace probe {

template <typename T, typename = void>
struct has_member_set_value : std::false_type {};

template <typename T>
struct has_member_set_value<
    T,
    std::void_t<decltype(std::declval<T&>().set_value(0.0))>
> : std::true_type {};

template <typename T, typename = void>
struct has_member_setValue : std::false_type {};

template <typename T>
struct has_member_setValue<
    T,
    std::void_t<decltype(std::declval<T&>().setValue(0.0))>
> : std::true_type {};

template <typename T, typename = void>
struct has_assign_from_double : std::false_type {};

template <typename T>
struct has_assign_from_double<
    T,
    std::void_t<decltype(std::declval<T&>() = 0.0)>
> : std::true_type {};

template <typename T, typename = void>
struct has_value_member_assign : std::false_type {};

template <typename T>
struct has_value_member_assign<
    T,
    std::void_t<decltype(std::declval<T&>().value = 0.0)>
> : std::true_type {};

template <typename T, typename = void>
struct has_value_m_member_assign : std::false_type {};

template <typename T>
struct has_value_m_member_assign<
    T,
    std::void_t<decltype(std::declval<T&>().value_m = 0.0)>
> : std::true_type {};

template <typename T, typename = void>
struct has_val_member_assign : std::false_type {};

template <typename T>
struct has_val_member_assign<
    T,
    std::void_t<decltype(std::declval<T&>().val = 0.0)>
> : std::true_type {};

template <typename T, typename = void>
struct has_val_m_member_assign : std::false_type {};

template <typename T>
struct has_val_m_member_assign<
    T,
    std::void_t<decltype(std::declval<T&>().val_m = 0.0)>
> : std::true_type {};

template <typename T, typename = void>
struct has_graph_member : std::false_type {};

template <typename T>
struct has_graph_member<
    T,
    std::void_t<decltype(std::declval<T&>().graph)>
> : std::true_type {};

template <typename T, typename = void>
struct has_graph_m_member : std::false_type {};

template <typename T>
struct has_graph_m_member<
    T,
    std::void_t<decltype(std::declval<T&>().graph_m)>
> : std::true_type {};

template <typename T, typename = void>
struct has_context_set_independent_values : std::false_type {};

template <typename T>
struct has_context_set_independent_values<
    T,
    std::void_t<decltype(
        std::declval<T&>().set_independent_values(
            std::declval<std::vector<double>>()
        )
    )>
> : std::true_type {};

template <typename T, typename = void>
struct has_context_set_values : std::false_type {};

template <typename T>
struct has_context_set_values<
    T,
    std::void_t<decltype(
        std::declval<T&>().set_values(
            std::declval<std::vector<double>>()
        )
    )>
> : std::true_type {};

template <typename T, typename = void>
struct has_context_forward : std::false_type {};

template <typename T>
struct has_context_forward<
    T,
    std::void_t<decltype(std::declval<T&>().forward())>
> : std::true_type {};

} // namespace probe

int main() {
    std::cout << "had_quadra reusable tape capability probe\n";

    using AD = quadra::AD;
    using TapeContext = quadra::TapeContext;

    std::cout << "AD has set_value(double): "
              << probe::has_member_set_value<AD>::value << "\n";

    std::cout << "AD has setValue(double): "
              << probe::has_member_setValue<AD>::value << "\n";

    std::cout << "AD assignable from double: "
              << probe::has_assign_from_double<AD>::value << "\n";

    std::cout << "AD has assignable .value member: "
              << probe::has_value_member_assign<AD>::value << "\n";

    std::cout << "AD has assignable .value_m member: "
              << probe::has_value_m_member_assign<AD>::value << "\n";

    std::cout << "AD has assignable .val member: "
              << probe::has_val_member_assign<AD>::value << "\n";

    std::cout << "AD has assignable .val_m member: "
              << probe::has_val_m_member_assign<AD>::value << "\n";

    std::cout << "TapeContext has .graph member: "
              << probe::has_graph_member<TapeContext>::value << "\n";

    std::cout << "TapeContext has .graph_m member: "
              << probe::has_graph_m_member<TapeContext>::value << "\n";

    std::cout << "TapeContext has set_independent_values(vector<double>): "
              << probe::has_context_set_independent_values<TapeContext>::value
              << "\n";

    std::cout << "TapeContext has set_values(vector<double>): "
              << probe::has_context_set_values<TapeContext>::value << "\n";

    std::cout << "TapeContext has forward(): "
              << probe::has_context_forward<TapeContext>::value << "\n";

    const bool has_direct_update =
        probe::has_member_set_value<AD>::value ||
        probe::has_member_setValue<AD>::value ||
        probe::has_value_member_assign<AD>::value ||
        probe::has_value_m_member_assign<AD>::value ||
        probe::has_val_member_assign<AD>::value ||
        probe::has_val_m_member_assign<AD>::value ||
        probe::has_context_set_independent_values<TapeContext>::value ||
        probe::has_context_set_values<TapeContext>::value;

    std::cout << "Reusable independent-value update appears exposed: "
              << has_direct_update << "\n";

    std::cout << "PASS: capability probe completed\n";

    return 0;
}
