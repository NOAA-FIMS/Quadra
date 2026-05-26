#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "../core/laplace/laplace_objective_cached.hpp"
#include "../core/model/quadra_model.hpp"

DECLARE_ADGRAPH();

class CorrelatedRandomInterceptModel
    : public quadra::QuadraModel<CorrelatedRandomInterceptModel> {
public:
  CorrelatedRandomInterceptModel(std::vector<double> y, std::vector<int> group,
                                 int n_groups)
      : y_m(std::move(y)), group_m(std::move(group)), n_groups_m(n_groups) {
    parameters_m.add("mu", 0.0, quadra::ParameterTransform::Identity, false);

    for (int g = 0; g < n_groups_m; ++g) {
      parameters_m.add("u_" + std::to_string(g), 0.0,
                       quadra::ParameterTransform::Identity, true);
    }
  }

  std::vector<std::string> parameter_names_impl() const {
    return parameters_m.names();
  }

  const quadra::ParameterSet &parameters() const { return parameters_m; }

  template <typename Type>
  Type evaluate_impl(const std::vector<Type> &p,
                     quadra::ModelReportContext &) const {
    Type mu = p[0];

    Type nll = Type(0.0);

    for (size_t i = 0; i < y_m.size(); ++i) {
      const int g = group_m[i];
      Type u = p[1 + g];

      Type r = Type(y_m[i]) - (mu + u);
      nll += Type(0.5) * r * r;
    }

    Type u0 = p[1];
    nll += Type(0.5) * u0 * u0;

    const double rho = 0.8;

    for (int g = 1; g < n_groups_m; ++g) {
      Type ug = p[1 + g];
      Type up = p[1 + g - 1];

      Type diff = ug - Type(rho) * up;
      nll += Type(0.5) * diff * diff;
    }

    return nll;
  }

private:
  std::vector<double> y_m;
  std::vector<int> group_m;
  int n_groups_m;
  quadra::ParameterSet parameters_m;
};

struct SimData {
  std::vector<double> y_m;
  std::vector<int> group_m;
};

SimData simulate_grouped_data(int G, int m) {
  std::mt19937 rng(1234);
  std::normal_distribution<double> e_dist(0.0, 1.0);

  SimData out;
  out.y_m.reserve(static_cast<size_t>(G * m));
  out.group_m.reserve(static_cast<size_t>(G * m));

  for (int g = 0; g < G; ++g) {
    const double ug = 0.25 * std::sin(0.1 * static_cast<double>(g));

    for (int i = 0; i < m; ++i) {
      out.y_m.push_back(5.0 + ug + e_dist(rng));
      out.group_m.push_back(g);
    }
  }

  return out;
}

template <typename F> double time_ms(F &&f) {
  const auto start = std::chrono::high_resolution_clock::now();
  f();
  const auto end = std::chrono::high_resolution_clock::now();

  return std::chrono::duration<double, std::milli>(end - start).count();
}

int main() {
  std::cout << "\nQuadra cached Laplace objective benchmark\n\n";

  std::cout << std::setw(8) << "G" << std::setw(10) << "m" << std::setw(12)
            << "evals" << std::setw(18) << "uncached ms" << std::setw(18)
            << "cached ms" << std::setw(14) << "speedup" << std::setw(14)
            << "nnz Huu" << "\n";

  std::cout << std::string(100, '-') << "\n";

  const int m = 20;
  const int evals = 50;

  for (int G : std::vector<int>{25, 50, 100, 250, 500}) {
    SimData data = simulate_grouped_data(G, m);

    CorrelatedRandomInterceptModel model(data.y_m, data.group_m, G);

    std::vector<double> theta = {4.5};
    std::vector<double> u0(static_cast<size_t>(G), 0.0);

    quadra::LaplaceObjectiveResult last_uncached;

    double uncached_ms = time_ms([&]() {
      std::vector<double> u_start = u0;

      for (int i = 0; i < evals; ++i) {
        theta[0] = 4.5 + 0.01 * static_cast<double>(i);

        last_uncached = quadra::evaluate_laplace_objective(
            model, theta, u_start, model.parameters());

        u_start = last_uncached.u_hat_m;
      }
    });

    quadra::CachedLaplaceObjectiveState cache_state;
    quadra::LaplaceObjectiveResult last_cached;

    double cached_ms = time_ms([&]() {
      std::vector<double> u_start = u0;

      for (int i = 0; i < evals; ++i) {
        theta[0] = 4.5 + 0.01 * static_cast<double>(i);

        last_cached = quadra::evaluate_laplace_objective_cached(
            model, theta, u_start, model.parameters(), cache_state);

        u_start = last_cached.u_hat_m;
      }
    });

    const double speedup = uncached_ms / cached_ms;

    std::cout << std::setw(8) << G << std::setw(10) << m << std::setw(12)
              << evals << std::setw(18) << std::fixed << std::setprecision(3)
              << uncached_ms << std::setw(18) << cached_ms << std::setw(14)
              << speedup << std::setw(14)
              << last_cached.hessian_random_m.nonZeros() << "\n";
  }

  return 0;
}
