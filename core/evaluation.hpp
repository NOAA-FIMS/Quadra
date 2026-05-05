#ifndef EVALUATION_HPP
#define EVALUATION_HPP

#pragma once
#include <iostream>
#include <vector>
#include "../model/parameter.hpp"
#include "../math/transforms.hpp"
#include "../core/autodiff.hpp"

namespace quadra
{

    struct ModelContext
    {
        std::vector<double> Y;
    };

    //====================================================
    // 1. Build full parameter vector
    //====================================================
    inline void build_parameters(
        const ParameterVector &params,
        const std::vector<AD> &u,
        std::vector<AD> &x)
    {
        size_t k = 0;

        for (size_t i = 0; i < params.size(); ++i)
        {

            if (params.params[i].is_random)
                x[i] = u[k++];
            else
                x[i] = AD(params.params[i].value);
        }
    }

    //====================================================
    // 2. Transform + Jacobian
    //====================================================
    inline AD apply_transforms(
        const ParameterVector &params,
        const std::vector<AD> &x,
        std::vector<AD> &t)
    {
        AD log_jacobian = 0.0;

        for (size_t i = 0; i < x.size(); ++i)
        {

            t[i] = apply_transform(
                x[i],
                params.params[i].transform);

            log_jacobian += apply_log_jacobian(
                x[i],
                params.params[i].transform);
        }

        return log_jacobian;
    }

    //====================================================
    // 3. SINGLE source of truth for NLL
    //====================================================
    template <typename Model>
    AD evaluate_nll(
        Model &model,
        const ParameterVector &params,
        const std::vector<AD> &u,
        bool include_jacobian = true)
    {

        std::vector<AD> x(params.size());
        std::vector<AD> t(params.size());
        std::cout << "u = ";
        for (auto &x : u)
            std::cout << x.val << " ";
        std::cout << "\n";

        build_parameters(params, u, x);

        AD log_jacobian =
            apply_transforms(params, x, t);
        std::cout << __FUNCTION__ << " called\n";
        AD nll = model(t);

        if (include_jacobian)
            nll -= log_jacobian;

        return nll;
    }

    template <typename F>
    AD evaluate(F &f, const std::vector<AD> &x)
    {
        return f(x);
    }

} // namespace quadra
#endif // EVALUATION_HPP