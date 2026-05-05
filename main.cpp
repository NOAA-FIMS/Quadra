
#include "core/autodiff.hpp"
#include "model/parameter.hpp"
#include "model/objective.hpp"
#include "core/laplace.hpp"
#include "core/optimizer.hpp"
#include "core/optimizer_lbfgs.hpp"

#include <iostream>
#include <chrono>
#include <iomanip>

DECLARE_ADGRAPH();

struct NormalModel
{
    template <typename T>
    T operator()(const std::vector<T> &p) const
    {
        T mu = p[0];
        T sigma = p[1];

        T y = 5.0;

        T nll = 0.0;
        nll += 0.5 * log(2.0 * M_PI);
        nll += log(sigma);
        nll += 0.5 * (y - mu) * (y - mu) / (sigma * sigma);

        return nll;
    }
};

// struct REModel {
//     double y;

//     template <typename T>
//     T operator()(const std::vector<T>& p) const {
//         T mu = p[0];     // fixed
//         T u  = p[1];     // random

//         T sigma = 1.0;

//         T nll = 0.0;

//         // likelihood
//         nll += 0.5 * (y - (mu + u))*(y - (mu + u));

//         // random effect prior
//         nll += 0.5 * u*u;

//         return nll;
//     }
// };

struct GaussianQuadratic
{

    template <typename T>
    T operator()(const std::vector<T> &p) const
    {

        T dx = p[0] - T(1.0);
        T dy = p[1] - T(2.0);

        return T(0.5) * (dx * dx + dy * dy);
    }
};

struct TestFunctor
{
    template <typename T>
    T operator()(const std::vector<T> &p)
    {
        // pelagia::ADContext ctx{};
        // had::g_ADGraph = &ctx.graph;
        //         std::cout << "graph ptr = " << had::g_ADGraph << std::endl;
        //         std::cout << "ctx ptr   = " << &ctx.graph << std::endl;
        T ret = p[0] + p[0];
        std::cout << "graph size = "
                  << had::g_ADGraph->vertices.size()
                  << std::endl;
        for (auto &x : p)
            std::cout << x.val << std::endl;
        // T dx = p[0] - 1.0;
        // T dy = p[1] - 2.0;
        // std::vector<T> grad(2);
        // grad[0] = dx;
        // grad[1] = dy;

        // T ret =  0.5 * (dx*dx + dy*dy);

        return ret;
    }
};

quadra::ParameterVector to_parameters(const Eigen::VectorXd &x)
{
    quadra::ParameterVector p;
    for (size_t i = 0; i < x.size(); ++i)
    {
        p.add({std::string("x") + std::to_string(i), x[i], quadra::Transform::Identity, true});
    }
    return p;
}
//  TestFunctor f;
// // pelagia::ParameterVector params;
// //  params.add({"x1", 0.0});
// //  params.add({"x2", 0.0});
// Eigen::VectorXd x(2);
//  x << 0, 0;
// pelagia::ParameterVector params = to_parameters(x);
// for (auto& p : params.params) {
//     p.is_random = false; // make all parameters fixed for this test
//     std::cout << "is_random=" << p.is_random
//               << " value=" << p.value << "\n";
// }
// // pelagia::LBFGSObjective<double> solver(f, x);

// auto res = pelagia::optimize_lbfgs(f, params);

struct REModel
{
    double y;

    template <typename T>
    T operator()(const std::vector<T> &p) const
    {
        T mu = p[0];
        T u0 = p[1];
        T u1 = p[2];

        T nll = 0.0;

        // likelihood
        T pred = mu + u0 + u1;
        nll += 0.5 * (y - pred) * (y - pred);

        // random effect prior
        nll += 0.5 * u0 * u0;
        nll += 0.5 * u1 * u1;

        return nll;
    }
};

struct REModel2
{
    double y;

    template <typename T>
    T operator()(const std::vector<T> &p) const
    {
        T mu = p[0];
        T u0 = p[1];
        T u1 = p[2];
        for (int i = 0; i < p.size(); ++i)
            std::cout << "p[" << i << "] = " << p[i].val << "\n";

        T nll = 0.0;

        // likelihood
        T pred = mu + u0 + u1;
        nll += 0.5 * (y - pred) * (y - pred);

        // random effect prior
        nll += 0.5 * u0 * u0;
        nll += 0.5 * u1 * u1;

        return nll;
    }
};
struct REModel2Nonlinear
{
    double y;

    template <typename T>
    T operator()(const std::vector<T> &p) const
    {
        T mu = p[0]; // fixed
        T u0 = p[1]; // random
        T u1 = p[2]; // random

        T eta = mu + u0 + 0.1 * sin(u1);
        // T eta = mu + u0 + 0.5 * u1 * u1;

        T nll = 0.0;

        // nonlinear likelihood
        nll += 0.5 * (y - eta) * (y - eta);

        // random effect priors
        nll += 0.5 * u0 * u0;
        nll += 0.5 * u1 * u1;

        return nll;
    }
};

struct ScaledREModel
{
    double y;
    int n;

    template <typename T>
    T operator()(const std::vector<T> &p) const
    {
        T mu = p[0];
        T nll = 0.0;

        T pred = mu;

        for (int i = 0; i < n; ++i)
        {
            T u = p[1 + i];
            pred += u / T(n);
            nll += 0.5 * u * u;
        }

        nll += 0.5 * (y - pred) * (y - pred);

        return nll;
    }
};

struct AR1REModel
{
    int n;    // number of random effects
    double y; // target level

    template <typename T>
    T operator()(const std::vector<T> &p) const
    {
        T mu = p[0]; // fixed effect

        T nll = 0.0;

        // Random effects live at p[1], ..., p[n]
        // Observation contribution: y ~ mu + u_i
        for (int i = 0; i < n; ++i)
        {
            T ui = p[1 + i];
            T err = y - (mu + ui);
            nll += 0.5 * err * err;
        }

        // AR(1)-style process penalty: u_i - u_{i-1}
        for (int i = 1; i < n; ++i)
        {
            T ui = p[1 + i];
            T uim1 = p[1 + i - 1];

            T diff = ui - uim1;
            nll += 0.5 * diff * diff;
        }

        // Weak prior on first state to anchor the chain
        {
            T u0 = p[1];
            nll += 0.5 * u0 * u0;
        }

        return nll;
    }
};

void test_scaled()
{
    std::vector<int> sizes = {1, 2, 5, 10, 25, 50, 100, 200, 500, 1000, 100000};
    int run = 1;
    for (int size : sizes)
    {

        std::cout << "\n"
                  << run++ << ":" << "Testing ScaledREModel with n_random = " << size << "\n\n";
        AR1REModel model{size, 5.0};

        quadra::ParameterVector params;

        // fixed parameter
        params.add({"mu", 0.0});

        // random effects
        for (int i = 0; i < size; ++i)
        {
            params.add({"u" + std::to_string(i),
                        0.0,
                        quadra::Transform::Identity,
                        true});
        }

        auto t0 = std::chrono::high_resolution_clock::now();

        auto res = quadra::optimize_lbfgs(model, params);

        auto t1 = std::chrono::high_resolution_clock::now();

        double ms =
            std::chrono::duration<double, std::milli>(t1 - t0).count();

        std::cout /*<< std::fixed << std::setprecision(12)*/ << "n_random = " << size
                                                             << ", time_ms = " << ms
                                                             << ", fx = " << res.value
                                                             << "\n";
    }
}

int main()
{

    test_scaled();
    // pelagia::ParameterVector params;

    // params.add({"mu", 0.0});                                     // fixed
    // params.add({"u0", 0.0, pelagia::Transform::Identity, true}); // random
    // params.add({"u1", 0.0, pelagia::Transform::Identity, true}); // random

    // REModel2Nonlinear model{5.0};

    // auto res = pelagia::optimize_lbfgs(model, params);

    // std::cout << "final value: " << res.value << "\n";
    // std::cout<< "final par: " << params.params[0].value << "\n";
    // std::cout<< "final par: " << params.params[1].value << "\n";

    return 0;
}