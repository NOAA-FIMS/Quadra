#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "../core/laplace/reusable_total_hdot_tape.hpp"
#include "../core/laplace/third_order_dense_hdot.hpp"
#include "../examples/big/catch_at_age_shared.hpp"
#include "../math/distributions.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

struct DenseDynamicObjective {
  int years;

  template <class T> T operator()(const std::vector<T> &x) const {
    constexpr int ages = 8;
    const T log_initial = x[0];
    const T mortality = exp(x[1]);
    const T slope = exp(x[2]);
    const T sigma = T(0.03) + exp(x[3]);
    const T fishing = exp(log_initial - T(2.0));
    const T a50 = T(1.0) + T(ages - 1) /
                                   (T(1.0) + exp(-log_initial));
    const T oldest_selectivity =
        T(1.0) / (T(1.0) + exp(-slope * (T(ages) - a50)));

    std::vector<T> state(ages);
    std::vector<T> next(ages);
    for (int age = 0; age < ages; ++age) {
      state[static_cast<std::size_t>(age)] =
          exp(log_initial) * exp(-mortality * T(age));
    }
    T nll = T(0.0);
    for (int year = 0; year < years; ++year) {
      const T random = x[10 + year];
      T prediction = T(0.0);
      std::vector<T> composition_weight(ages);
      for (int age = 0; age < ages; ++age) {
        const T selectivity =
            (T(1.0) /
             (T(1.0) + exp(-slope * (T(age + 1) - a50)))) /
            oldest_selectivity;
        const T total_mortality = mortality + fishing * selectivity;
        composition_weight[static_cast<std::size_t>(age)] =
            state[static_cast<std::size_t>(age)] *
            (fishing * selectivity / total_mortality) *
            (T(1.0) - exp(-total_mortality));
        prediction += composition_weight[static_cast<std::size_t>(age)];
      }
      prediction += T(1.0e-8);
      const T residual =
          (log(prediction) - T(1.2 + 0.08 * year)) / sigma;
      nll += T(0.5) * residual * residual + log(sigma) +
             T(0.4) * random * random;
      const std::vector<int> counts{6 + year, 14, 28, 42, 45, 32, 20, 13};
      nll -= quadra::ddirichlet_multinomial_robust(
          counts, composition_weight, T(35.0), true);

      std::fill(next.begin(), next.end(), T(0.0));
      next[0] = exp(log_initial + random - T(0.5) * sigma * sigma);
      for (int age = 1; age < ages; ++age) {
        const T selectivity =
            (T(1.0) / (T(1.0) + exp(-slope * (T(age) - a50)))) /
            oldest_selectivity;
        next[static_cast<std::size_t>(age)] =
            state[static_cast<std::size_t>(age - 1)] *
            exp(-(mortality + fishing * selectivity));
      }
      next[ages - 1] += state[ages - 1] *
                        exp(-(mortality + fishing));
      state.swap(next);
    }
    return nll;
  }
};

template <class Objective>
Eigen::MatrixXd random_hessian(const Objective &objective,
                               const Eigen::VectorXd &theta,
                               const Eigen::VectorXd &random) {
  had::ADGraph graph;
  had::g_ADGraph = &graph;
  std::vector<had::AReal> x;
  x.reserve(static_cast<std::size_t>(theta.size() + random.size()));
  for (int i = 0; i < theta.size(); ++i) x.emplace_back(theta[i]);
  for (int i = 0; i < random.size(); ++i) x.emplace_back(random[i]);
  had::AReal y = objective(x);
  had::SetAdjoint(y, 1.0);
  had::PropagateAdjoint();

  Eigen::MatrixXd hessian(random.size(), random.size());
  for (int row = 0; row < random.size(); ++row) {
    for (int col = 0; col < random.size(); ++col) {
      hessian(row, col) = had::GetAdjoint(
          x[static_cast<std::size_t>(theta.size() + row)],
          x[static_cast<std::size_t>(theta.size() + col)]);
    }
  }
  had::g_ADGraph = nullptr;
  return hessian;
}

void run_catch_at_age_log_m_test() {
  example::CatchAtAgeLaplaceModel model;
  model.data.n_years = 6;
  auto objective = [&model](const auto &x) { return model(x); };
  constexpr int theta_dim = 10;
  constexpr int random_dim = 6;
  Eigen::VectorXd theta(theta_dim);
  theta << std::log(900.0), std::log(0.25), std::log(0.15), std::log(0.18),
      0.0, std::log(1.25), std::log(0.20), std::log(0.15), std::log(0.35),
      std::log(40.0);
  const Eigen::VectorXd random = Eigen::VectorXd::Zero(random_dim);
  quadra::laplace::RandomHessianPattern pattern;
  for (int row = 0; row < random_dim; ++row) {
    for (int col = 0; col <= row; ++col) pattern.emplace_back(row, col);
  }
  quadra::laplace::ReusableTotalHdotTape<decltype(objective)> tape(
      objective, theta_dim, random_dim, pattern, {1}, theta, random);
  const auto analytic = tape.compute(theta, random, [random_dim](int) {
    return Eigen::VectorXd::Zero(random_dim);
  });
  const double step = 1.0e-5;
  Eigen::VectorXd plus = theta;
  Eigen::VectorXd minus = theta;
  plus[1] += step;
  minus[1] -= step;
  const Eigen::MatrixXd finite_difference =
      (random_hessian(objective, plus, random) -
       random_hessian(objective, minus, random)) /
      (2.0 * step);
  std::vector<double> combined(theta.size() + random.size());
  std::vector<double> total_direction(combined.size(), 0.0);
  std::vector<int> random_indices(random.size());
  for (int i = 0; i < theta.size(); ++i) combined[i] = theta[i];
  for (int i = 0; i < random.size(); ++i) {
    combined[theta.size() + i] = random[i];
    random_indices[i] = theta.size() + i;
  }
  total_direction[1] = 1.0;
  const Eigen::MatrixXd third_order =
      quadra::laplace::dense_hdot_third_order_polarized(
          objective, combined, total_direction, random_indices);
  const double relative_error =
      (Eigen::MatrixXd(analytic[1]) - finite_difference).norm() /
      std::max(1.0e-12, finite_difference.norm());
  Eigen::Index max_row = 0;
  Eigen::Index max_col = 0;
  const Eigen::MatrixXd difference =
      Eigen::MatrixXd(analytic[1]) - finite_difference;
  difference.cwiseAbs().maxCoeff(&max_row, &max_col);
  std::cout << "catch-at-age log_M relative Hdot error: " << relative_error
            << ", max entry (" << max_row << ", " << max_col
            << ") analytic=" << Eigen::MatrixXd(analytic[1])(max_row, max_col)
            << " fd=" << finite_difference(max_row, max_col) << "\n";
  std::cout << "  third-order polarization relative error: "
            << (third_order - finite_difference).norm() /
                   std::max(1.0e-12, finite_difference.norm()) << "\n";
  const double third_order_error =
      (third_order - finite_difference).norm() /
      std::max(1.0e-12, finite_difference.norm());
  if (third_order_error > 1.0e-6)
    throw std::runtime_error("exact dense Hdot differs from validation FD");
}

void run_test() {
  constexpr int theta_dim = 10;
  constexpr int random_dim = 6;
  DenseDynamicObjective objective{random_dim};
  Eigen::VectorXd theta(theta_dim);
  theta << std::log(2.0), std::log(0.35), std::log(0.8), std::log(0.2), 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0;
  Eigen::VectorXd random = Eigen::VectorXd::LinSpaced(random_dim, -0.2, 0.25);

  quadra::laplace::RandomHessianPattern pattern;
  for (int row = 0; row < random_dim; ++row) {
    for (int col = 0; col <= row; ++col) pattern.emplace_back(row, col);
  }
  std::vector<int> directions{1};
  quadra::laplace::ReusableTotalHdotTape<DenseDynamicObjective> tape(
      objective, theta_dim, random_dim, pattern, directions, theta, random);
  const auto analytic = tape.compute(theta, random, [](int) {
    return Eigen::VectorXd::Zero(random_dim);
  });

  for (int j : directions) {
    const double step = 1.0e-5 * std::max(1.0, std::abs(theta[j]));
    Eigen::VectorXd plus = theta;
    Eigen::VectorXd minus = theta;
    plus[j] += step;
    minus[j] -= step;
    const Eigen::MatrixXd finite_difference =
        (random_hessian(objective, plus, random) -
         random_hessian(objective, minus, random)) /
        (2.0 * step);
    const Eigen::MatrixXd analytic_dense = Eigen::MatrixXd(analytic[j]);
    const double relative_error =
        (analytic_dense - finite_difference).norm() /
        std::max(1.0e-12, finite_difference.norm());
    std::cout << "direction " << j
              << " relative Hdot error: " << relative_error << "\n";
    if (relative_error > 1.0e-6) {
      throw std::runtime_error(
          "dense dynamic analytic Hdot differs from finite difference");
    }
  }
}

} // namespace

int main() {
  run_test();
  run_catch_at_age_log_m_test();
  std::cout << "PASS: exact dense analytic Hdot validation\n";
  return 0;
}
