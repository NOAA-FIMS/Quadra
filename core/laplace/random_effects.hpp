#ifndef QUADRA_RANDOM_EFFECTS_HPP
#define QUADRA_RANDOM_EFFECTS_HPP
#pragma once

/**
 * @file random_effects.hpp
 * @brief Random-effect optimization interface.
 *
 * This module owns computation of the conditional mode
 *
 * \f[
 *   \hat u(\theta) = \arg\min_u f(\theta, u).
 * \f]
 *
 * Current implementation note:
 * the implementation still lives in `laplace.hpp` during the transition.
 * This header establishes the public module boundary for the next refactor.
 */

namespace quadra
{
    // Implementation currently provided by core/laplace/laplace.hpp:
    //   solve_random_effects_laplace(...)
}

#endif
