#ifndef QUADRA_AUTODIFF_HPP
#define QUADRA_AUTODIFF_HPP

#include <vector>
#include "had/had.h"
#include "eigen/Eigen/Dense"

namespace quadra
{

    using AD = had::AReal;

    struct ADResult
    {
        double value;
        Eigen::VectorXd gradient;
        Eigen::MatrixXd hessian;
    };

    struct ADScope
    {
        had::ADGraph &graph;

        explicit ADScope(had::ADGraph &g) : graph(g)
        {
            had::g_ADGraph = &graph;
            graph.Clear();
        }

        // forward value only
        template <typename T>
        double value(const T &loss) const
        {
            return loss.val;
        }

        // run reverse pass once
        void backward(const had::AReal &loss)
        {
            had::SetAdjoint(loss, 1.0);
            had::PropagateAdjoint();
        }

        // first-order
        double grad(const had::AReal &x) const
        {
            return had::GetAdjoint(x);
        }

        // second-order
        double hess(const had::AReal &x, const had::AReal &y) const
        {
            return had::GetAdjoint(x, y);
        }

        void clear()
        {
            graph.Clear();
        }
    };
    // struct ADEngine
    // {
    //     had::ADGraph &graph;

    //     explicit ADEngine(had::ADGraph &g) : graph(g)
    //     {
    //         had::g_ADGraph = &graph;
    //         graph.Clear();
    //     }

    //     had::AReal eval(const std::function<had::AReal()> &f)
    //     {
    //         return f();
    //     }

    //     void backward(const had::AReal &loss)
    //     {
    //         had::SetAdjoint(loss, 1.0);
    //         had::PropagateAdjoint();
    //     }

    //     Eigen::VectorXd gradient(const std::vector<had::AReal> &vars) const
    //     {
    //         Eigen::VectorXd g(vars.size());
    //         for (size_t i = 0; i < vars.size(); ++i)
    //             g[i] = had::GetAdjoint(vars[i]);
    //         return g;
    //     }

    //     Eigen::MatrixXd hessian(const std::vector<had::AReal> &vars) const
    //     {
    //         Eigen::MatrixXd H(vars.size(), vars.size());
    //         H.setZero();

    //         for (size_t i = 0; i < vars.size(); ++i)
    //             for (size_t j = 0; j < vars.size(); ++j)
    //                 H(i, j) = had::GetAdjoint(vars[i], vars[j]);

    //         return H;
    //     }
    // };

    //--------------------------------------------------
    // Convert Eigen -> std::vector<AD>
    //--------------------------------------------------
    inline std::vector<AD> to_ad(const Eigen::VectorXd &x)
    {
        std::vector<AD> out;
        out.reserve(x.size());

        for (int i = 0; i < x.size(); ++i)
        {
            out.emplace_back(x[i]);
        }

        return out;
    }

    inline std::vector<had::AReal> to_ad(const std::vector<double> &x)
    {
        std::vector<had::AReal> out;
        out.reserve(x.size());

        for (double v : x)
            out.emplace_back(v);

        return out;
    }

    // class ADSession
    // {
    // public:
    //     ADSession()
    //     {
    //         init_graph();
    //     }

    //     ~ADSession()
    //     {
    //         clear_graph();
    //     }

    //     ADSession(const ADSession &) = delete;
    //     ADSession &operator=(const ADSession &) = delete;

    //     template <typename F>
    //     ADResult run(F &&f, const std::vector<double> &x)
    //     {
    //         reset_graph();

    //         // build inputs
    //         std::vector<had::AReal> x_ad;
    //         x_ad.reserve(x.size());

    //         for (double v : x)
    //             x_ad.emplace_back(v);

    //         // forward
    //         had::AReal y = f(x_ad);

    //         // reverse
    //         had::SetAdjoint(y, 1.0);
    //         had::PropagateAdjoint();

    //         // extract
    //         ADResult res;
    //         res.value = y.val;
    //         res.gradient.resize(x.size());

    //         for (size_t i = 0; i < x_ad.size(); i++)
    //         {
    //             res.gradient[i] = had::GetAdjoint(x_ad[i]);
    //         }

    //         return res;
    //     }

    // private:
    //     void init_graph()
    //     {
    //         if (!had::g_ADGraph)
    //         {
    //             static had::ADGraph global_graph;
    //             had::g_ADGraph = &global_graph;
    //         }
    //     }

    //     void reset_graph()
    //     {
    //         had::g_ADGraph->Clear();
    //     }

    //     void clear_graph()
    //     {
    //         if (had::g_ADGraph)
    //             had::g_ADGraph->Clear();
    //     }
    // };

    //--------------------------------------------------
    // Extract gradients AFTER reverse pass
    //--------------------------------------------------
    inline Eigen::VectorXd extract_gradient(const std::vector<AD> &x)
    {
        Eigen::VectorXd g(x.size());

        for (size_t i = 0; i < x.size(); ++i)
        {
            g[i] = had::GetAdjoint(x[i]);
        }

        return g;
    }

    //--------------------------------------------------
    // Extract Hessian element
    //--------------------------------------------------
    inline double get_hessian(const AD &xi, const AD &xj)
    {
        return had::GetAdjoint(xi, xj);
    }
} // namespace quadra
#endif