#ifndef QUADRA_HPP
#define QUADRA_HPP
#pragma once

/**
 * @file quadra.hpp
 * @brief Public umbrella header for the Quadra inference framework.
 *
 * Most user-facing code should be able to include this single header.
 */

#include "model/parameter.hpp"

#include "core/autodiff/autodiff.hpp"
#include "core/laplace/laplace_all.hpp"
#include "core/optimizer/optimizer.hpp"
#include "core/inference/inference.hpp"
#include "core/sparse/sparse.hpp"

#endif
