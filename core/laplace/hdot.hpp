#ifndef QUADRA_HDOT_HPP
#define QUADRA_HDOT_HPP
#pragma once

/**
 * @file hdot.hpp
 * @brief Exact directional Hessian propagation interface.
 *
 * Laplace log-determinant derivatives require
 *
 * \f[
 *   \dot H_i
 *   =
 *   D H_{uu}
 *   \left[
 *     e_i,\frac{d\hat u}{d\theta_i}
 *   \right].
 * \f]
 *
 * This module owns exact directional Hessian propagation.
 *
 * Current implementation note:
 * the implementation still lives in `laplace.hpp` during the transition.
 */

namespace quadra
{
    // Implementation currently provided by core/laplace/laplace.hpp:
    //   random_hessian_directional_exact(...)
}

#endif
