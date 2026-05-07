#ifndef QUADRA_PARAMETER_HPP
#define QUADRA_PARAMETER_HPP

#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace quadra
{

    enum class Transform
    {
        Identity,
        Log,
        Logit,
        Square
    };

    template <typename T>
    inline T apply_transform(const T &x, Transform transform)
    {
        switch (transform)
        {
        case Transform::Identity:
            return x;
        case Transform::Log:
            return exp(x);
        case Transform::Logit:
            return T(1.0) / (T(1.0) + exp(-x));
        case Transform::Square:
            return x * x;
        default:
            throw std::runtime_error("Unknown parameter transform");
        }
    }

    template <typename T>
    inline T log_jacobian(const T &x, Transform transform)
    {
        switch (transform)
        {
        case Transform::Identity:
            return T(0.0);
        case Transform::Log:
            return x;
        case Transform::Logit:
        {
            T y = T(1.0) / (T(1.0) + exp(-x));
            return log(y) + log(T(1.0) - y);
        }
        case Transform::Square:
            return log(T(2.0) * x);
        default:
            throw std::runtime_error("Unknown parameter transform");
        }
    }

    struct Parameter
    {
        std::string name;
        double value = 0.0;
        Transform transform = Transform::Identity;
        bool is_random = false;

        // Optimization/reporting metadata
        bool active = true;
        bool report = true;

        // Post-fit covariance metadata
        bool estimate_covariance = true;
        double std_error = std::numeric_limits<double>::quiet_NaN();
        int covariance_index = -1;

        Parameter() = default;

        Parameter(const std::string &n,
                  double v,
                  Transform t = Transform::Identity,
                  bool random = false)
            : name(n),
              value(v),
              transform(t),
              is_random(random),
              estimate_covariance(!random)
        {
        }
    };

    struct ParameterVector
    {
        std::vector<Parameter> params;

        void add(const Parameter &p)
        {
            params.push_back(p);
        }

        int size() const
        {
            return static_cast<int>(params.size());
        }

        Parameter &operator[](std::size_t i)
        {
            return params[i];
        }

        const Parameter &operator[](std::size_t i) const
        {
            return params[i];
        }
    };

    // inline std::vector<int> build_fixed_index(const ParameterVector &params)
    // {
    //     std::vector<int> idx;
    //     for (int i = 0; i < params.size(); ++i)
    //     {
    //         if (!params.params[static_cast<std::size_t>(i)].is_random &&
    //             params.params[static_cast<std::size_t>(i)].active)
    //         {
    //             idx.push_back(i);
    //         }
    //     }
    //     return idx;
    // }

    // inline std::vector<int> build_random_index(const ParameterVector &params)
    // {
    //     std::vector<int> idx;
    //     for (int i = 0; i < params.size(); ++i)
    //     {
    //         if (params.params[static_cast<std::size_t>(i)].is_random &&
    //             params.params[static_cast<std::size_t>(i)].active)
    //         {
    //             idx.push_back(i);
    //         }
    //     }
    //     return idx;
    // }

} // namespace quadra

#endif
