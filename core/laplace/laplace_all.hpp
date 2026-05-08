#ifndef QUADRA_LAPLACE_ALL_HPP
#define QUADRA_LAPLACE_ALL_HPP
#pragma once

/**
 * @file laplace_all.hpp
 * @brief Umbrella header for Laplace approximation machinery.
 *
 * The Laplace layer evaluates
 *
 * \f[
 *   \tilde f(\theta)
 *   =
 *   f(\theta, \hat u(\theta))
 *   +
 *   \frac{1}{2}\log\det H_{uu}
 *   -
 *   \frac{n_u}{2}\log(2\pi).
 * \f]
 */

#include "options.hpp"
#include "evaluation.hpp"
#include "random_effects.hpp"
#include "implicit_diff.hpp"
#include "hdot.hpp"
#include "laplace.hpp"

#endif
