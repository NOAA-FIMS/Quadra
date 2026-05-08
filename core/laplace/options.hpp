#ifndef QUADRA_LAPLACE_OPTIONS_HPP
#define QUADRA_LAPLACE_OPTIONS_HPP
#pragma once

#include "../sparse/factorization.hpp"
#include "../sparse/trace.hpp"

namespace quadra {

/**
 * @brief Options controlling Laplace approximation.
 *
 * Quadra uses
 *
 * \f[
 *   \tilde f(\theta)
 *   =
 *   f(\theta, \hat u(\theta))
 *   +
 *   \frac{1}{2}\log\det H_{uu}
 *   -
 *   \frac{n_u}{2}\log(2\pi),
 * \f]
 *
 * with
 *
 * \f[
 *   H_{uu}
 *   =
 *   \frac{\partial^2 f}{\partial u \partial u}.
 * \f]
 */
struct LaplaceOptions {
    /// Trace strategy for \f$\operatorname{tr}(H^{-1}\dot H)\f$.
    TraceOptions trace;

    /// Sparse factorization stabilization options.
    SparseOptions sparse;

    /// Runtime opt-in for Hdot validation when compiled with QUADRA_VALIDATE_HDOT.
    bool validate_hdot = true;

    /// Sparse Hessian entry drop tolerance.
    double hessian_drop_tol = 0.0;

    // Compatibility fields for current examples/code.
    bool& use_hutchinson_trace = trace.use_hutchinson_trace;
    int& hutchinson_probes = trace.hutchinson_probes;
    unsigned int& hutchinson_seed = trace.hutchinson_seed;
    double& jitter_initial = sparse.jitter_initial;
    int& jitter_max_attempts = sparse.jitter_max_attempts;
};

/**
 * @brief Mutable default Laplace options.
 *
 * Prefer passing explicit options in production paths.
 */
inline LaplaceOptions& default_laplace_options()
{
    static LaplaceOptions options;
    return options;
}

} // namespace quadra

#endif
