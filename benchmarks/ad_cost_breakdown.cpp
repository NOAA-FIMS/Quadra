#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include "../core/model/quadra_model.hpp"
#include "../core/laplace/random_effect_hessian.hpp"

DECLARE_ADGRAPH();

class RandomInterceptModel : public quadra::QuadraModel<RandomInterceptModel> {
public:
    explicit RandomInterceptModel(std::vector<double> y)
        : y_m(std::move(y))
    {
        parameters_m.add("mu", 0.0, quadra::ParameterTransform::Identity, false);
        parameters_m.add("u", 0.0, quadra::ParameterTransform::Identity, true);
    }

    std::vector<std::string> parameter_names_impl() const {
        return parameters_m.names();
    }

    const quadra::ParameterSet& parameters() const {
        return parameters_m;
    }

    template <typename Type>
    Type evaluate_impl(
        const std::vector<Type>& p,
        quadra::ModelReportContext&
    ) const {
        Type mu = p[0];
        Type u  = p[1];

        Type nll = Type(0.0);

        for (double yi : y_m) {
            Type r = Type(yi) - (mu + u);
            nll += Type(0.5) * r * r;
        }

        nll += Type(0.5) * u * u;

        return nll;
    }

private:
    std::vector<double> y_m;
    quadra::ParameterSet parameters_m;
};

std::vector<double> simulate_data(size_t n) {
    std::mt19937 rng(1234);
    std::normal_distribution<double> dist(5.25, 1.0);

    std::vector<double> y(n);

    for (size_t i = 0; i < n; ++i) {
        y[i] = dist(rng);
    }

    return y;
}

template <typename F>
double time_ms(F&& f) {
    auto start = std::chrono::high_resolution_clock::now();
    f();
    auto end = std::chrono::high_resolution_clock::now();

    return std::chrono::duration<double, std::milli>(end - start).count();
}

int main() {

    std::cout << "\nQuadra AD cost breakdown benchmark\n\n";

    std::cout
        << std::setw(8)  << "n"
        << std::setw(16) << "double eval"
        << std::setw(16) << "tape build"
        << std::setw(16) << "reverse grad"
        << std::setw(16) << "Hessian"
        << std::setw(16) << "total AD"
        << "\n";

    std::cout << std::string(88, '-') << "\n";

    for (size_t n : std::vector<size_t>{10, 100, 1000, 5000, 10000}) {

        RandomInterceptModel model(simulate_data(n));

        quadra::ModelReportContext ctx;

        std::vector<double> p = {5.25, 0.0};

        double obj = 0.0;

        double eval_ms = time_ms([&]() {
            obj = model.evaluate<double>(p, ctx);
        });

        quadra::TapeContext tape;
        quadra::ADScope scope(tape.graph);

        std::vector<quadra::AD> ad_params;

        double tape_ms = time_ms([&]() {
            ad_params = quadra::to_ad(p);
            quadra::AD joint = model.evaluate<quadra::AD>(ad_params, ctx);
            (void)joint;
        });

        double reverse_ms = time_ms([&]() {
            quadra::TapeContext tape2;
            quadra::ADScope scope2(tape2.graph);

            auto ad2 = quadra::to_ad(p);

            quadra::AD joint =
                model.evaluate<quadra::AD>(ad2, ctx);

            scope2.backward(joint);

            auto g = quadra::extract_gradient(ad2);
            (void)g;
        });

        auto partition =
            quadra::partition_parameters(model.parameters());

        double hessian_ms = time_ms([&]() {

            auto h =
                quadra::evaluate_random_effect_hessian(
                    model,
                    {5.25},
                    {0.0},
                    partition
                );

            (void)h;
        });

        double total_ad = tape_ms + reverse_ms + hessian_ms;

        std::cout
            << std::setw(8)  << n
            << std::setw(16) << std::fixed << std::setprecision(3) << eval_ms
            << std::setw(16) << tape_ms
            << std::setw(16) << reverse_ms
            << std::setw(16) << hessian_ms
            << std::setw(16) << total_ad
            << "\n";
    }

    return 0;
}
