#ifndef QUADRA_IMPLICIT_DIFF_HPP
#define QUADRA_IMPLICIT_DIFF_HPP
#pragma once

/**
 * @file implicit_diff.hpp
 * @brief Implicit differentiation utilities for optimized random effects.
 *
 * If
 *
 * \f[
 *   f_u(\theta, \hat u(\theta)) = 0,
 * \f]
 *
 * then
 *
 * \f[
 *   \frac{d\hat u}{d\theta_i}
 *   =
 *   -H_{uu}^{-1}H_{u\theta_i}.
 * \f]
 *
 * Current implementation note:
 * the implementation still lives in `laplace.hpp` during the transition.
 */

namespace quadra
{
    // Implementation currently provided by core/laplace/laplace.hpp:
    //   implicit_du_dtheta_i(...)
    //   implicit_du_dtheta_all(...)
}

#endif
