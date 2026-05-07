
#ifndef QUADRA_TEST_COMMON_HPP
#define QUADRA_TEST_COMMON_HPP

#include <cmath>
#include <iostream>
#include <vector>

#include "../core/optimizer/optimizer.hpp"

namespace quadra_tests {

inline void print_banner(const std::string& name) {
    std::cout << "\n==================================================\n";
    std::cout << name << "\n";
    std::cout << "==================================================\n";
}

inline quadra::LaplaceOptions default_test_options() {
    quadra::LaplaceOptions opts;
    opts.use_hutchinson_trace = false; // deterministic for small tests
    opts.hutchinson_probes = 8;
    opts.jitter_initial = 1e-12;
    opts.jitter_max_attempts = 12;
    opts.hessian_drop_tol = 0.0;
    return opts;
}

inline quadra::LaplaceOptions large_test_options() {
    quadra::LaplaceOptions opts;
    opts.use_hutchinson_trace = true;
    opts.hutchinson_probes = 8;
    opts.hutchinson_seed = 12345;
    opts.jitter_initial = 1e-12;
    opts.jitter_max_attempts = 12;
    opts.hessian_drop_tol = 0.0;
    return opts;
}

} // namespace quadra_tests

#endif
