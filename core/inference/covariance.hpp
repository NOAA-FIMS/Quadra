#ifndef QUADRA_COVARIANCE_HPP
#define QUADRA_COVARIANCE_HPP

#pragma once

#include "../eigen/Eigen/Dense"
#include "../../model/parameter.hpp"

namespace quadra
{

    struct CovarianceOptions
    {
        bool include_fixed_effects = true;
        bool include_random_effects = false;
        bool use_laplace_adjusted_hessian = true;
    };

    struct CovarianceResult
    {
        Eigen::MatrixXd covariance;
        std::vector<int> parameter_indices;
    };

    inline void store_standard_errors(ParameterVector &params,
                                      const CovarianceResult &result)
    {
        for (Eigen::Index i = 0; i < result.covariance.rows(); ++i)
        {
            const int pidx = result.parameter_indices[static_cast<std::size_t>(i)];
            params.params[static_cast<std::size_t>(pidx)].covariance_index =
                static_cast<int>(i);

            const double v = result.covariance(i, i);
            params.params[static_cast<std::size_t>(pidx)].std_error =
                (v >= 0.0) ? std::sqrt(v) : std::numeric_limits<double>::quiet_NaN();
        }
    }

} // namespace quadra

#endif
