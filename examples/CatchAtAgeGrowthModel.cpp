#include "../core/autodiff.hpp"
#include "../model/parameter.hpp"
#include "../model/objective.hpp"
#include "../core/laplace.hpp"
#include "../core/optimizer.hpp"
#include "../math/distributions.hpp"
#include <iostream>
#include <chrono>
#include <iomanip>

DECLARE_ADGRAPH();

struct VBGFGrowthOnly
{
    std::vector<int> ages;
    std::vector<double> obs_len;

    template <typename T>
    T operator()(const std::vector<T> &p) const
    {
        T Linf = exp(p[0]);
        T K = exp(p[1]);
        T t0 = p[2];
        T sigma_len = exp(p[3]) + T(1e-6);

        T nll = 0.0;

        for (size_t i = 0; i < obs_len.size(); ++i)
        {
            T age = T(static_cast<double>(ages[i]));
            T pred = Linf * (T(1.0) - exp(-K * (age - t0)));
            nll -= quadra::dnorm(T(obs_len[i]), pred, sigma_len, true);
        }

        return nll;
    }
};
struct CatchAtAgeModel
{
    int n_years;
    int n_ages;

    std::vector<double> obs_catch_at_age;
    std::vector<double> obs_len_at_age;

    double M = 0.2;

    template <typename T>
    T operator()(const std::vector<T> &p) const
    {
        T Linf = exp(p[0]);
        T K = exp(p[1]);
        T t0 = p[2];
        T R0 = exp(p[3]);
        T sigma_rec = exp(p[4]);
        T sigma_len = exp(p[5]);
        T log_sigma_caa = p[6];
        T sigma_caa = exp(log_sigma_caa);
        T F = exp(p[7]);

        const int rec_start = 8;

        T nll = 0.0;

        std::vector<T> N(n_years * n_ages, T(0.0));

        // recruitment
        for (int y = 0; y < n_years; ++y)
        {
            T rec_dev = p[rec_start + y];

            nll -= quadra::dnorm(rec_dev, T(0.0), sigma_rec, true);

            N[y * n_ages + 0] =
                R0 * exp(rec_dev - T(0.5) * sigma_rec * sigma_rec);
        }

        // propagate numbers-at-age with total mortality Z = M + F
        T Z = T(M) + F;

        for (int y = 1; y < n_years; ++y)
        {
            for (int a = 1; a < n_ages; ++a)
            {
                N[y * n_ages + a] =
                    N[(y - 1) * n_ages + (a - 1)] * exp(-Z);
            }

            // plus group
            N[y * n_ages + (n_ages - 1)] +=
                N[(y - 1) * n_ages + (n_ages - 1)] * exp(-Z);
        }

        // catch-at-age likelihood using Baranov catch equation
        for (int y = 0; y < n_years; ++y)
        {
            for (int a = 0; a < n_ages; ++a)
            {
                int idx = y * n_ages + a;

                T pred_catch =
                    N[idx] * (F / Z) * (T(1.0) - exp(-Z));

                T obs_log = log(T(obs_catch_at_age[idx]) + T(1e-12));
                T pred_log = log(pred_catch + T(1e-12));

                nll -= quadra::dnorm(obs_log, pred_log, sigma_caa, true);
            }
        }

        // length-at-age likelihood
        for (int y = 0; y < n_years; ++y)
        {
            for (int a = 0; a < n_ages; ++a)
            {
                int idx = y * n_ages + a;
                T age = T(a + 1.0);

                T pred_len =
                    Linf * (T(1.0) - exp(-K * (age - t0)));

                nll -= quadra::dnorm(
                    T(obs_len_at_age[idx]),
                    pred_len,
                    sigma_len,
                    true);
            }
        }

        nll += 0.5 * (log_sigma_caa - log(T(0.3))) *
               (log_sigma_caa - log(T(0.3))) / T(0.25);

        return nll;
    }
};

int main()
{
    CatchAtAgeModel model;
    model.n_years = 10;
    model.n_ages = 5;
    model.obs_catch_at_age = {
        // y0
        1000, 700, 500, 350, 250,

        // y1
        1100, 820, 570, 400, 280,

        // y2
        950, 900, 650, 450, 320,

        // y3
        1200, 780, 710, 520, 360,

        // y4
        1050, 980, 620, 560, 410,

        // y5
        1300, 860, 800, 510, 460,

        // y6
        1150, 1080, 710, 660, 430,

        // y7
        1400, 950, 900, 590, 540,

        // y8
        1250, 1170, 800, 740, 500,

        // y9
        1500, 1030, 980, 680, 620};

    model.obs_len_at_age = {
        20.5, 34.1, 45.8, 54.2, 60.3,
        21.3, 35.0, 46.4, 55.1, 61.2,
        19.8, 33.5, 45.1, 53.7, 60.1,
        20.9, 34.7, 46.0, 54.6, 61.0,
        21.1, 35.2, 46.8, 55.0, 61.5,
        20.2, 33.9, 45.4, 54.0, 60.6,
        20.7, 34.5, 46.1, 54.8, 61.1,
        21.0, 35.1, 46.7, 55.3, 61.7,
        20.4, 34.0, 45.5, 54.1, 60.8,
        20.8, 34.6, 46.2, 54.9, 61.3};

    quadra::ParameterVector params;

    params.add({"log_Linf", std::log(75.0)});
    params.add({"log_K", std::log(0.25)});
    params.add({"t0", -0.2});
    params.add({"log_R0", std::log(1000.0)});
    params.add({"log_sigma_rec", std::log(0.4)});
    params.add({"log_sigma_len", std::log(3.0)});
    params.add({"log_sigma_caa", std::log(0.3)});
    params.add({"log_F", std::log(0.2)});

    for (int y = 0; y < model.n_years; ++y)
    {
        params.add({"rec_dev_" + std::to_string(y),
                    0.0,
                    quadra::Transform::Identity,
                    true});
    }

    auto t0 = std::chrono::high_resolution_clock::now();

    auto res = quadra::optimize_lbfgs(model, params);

    auto t1 = std::chrono::high_resolution_clock::now();

    double ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count();

    std::cout << "final parameters:\n";
    for (int i = 0; i < params.size(); ++i)
    {
        std::cout << "Parameter index " << i << ": " << params.params[i].name
                  << " = " << params.params[i].value
                  << "\n";
    }

    std::cout << "Linf = " << std::exp(params.params[0].value) << "\n";
    std::cout << "K = " << std::exp(params.params[1].value) << "\n";
    std::cout << "t0 = " << params.params[2].value << "\n";
    std::cout << "sigma_len = "
              << std::exp(params.params[5].value)
              << "\n";
}