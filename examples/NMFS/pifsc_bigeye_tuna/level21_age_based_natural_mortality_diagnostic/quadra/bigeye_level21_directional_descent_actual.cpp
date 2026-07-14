#include "../../../../../core/optimizer.hpp"
#include "../objective/bigeye_quadra_objective.hpp"
#include "../quadra/bigeye_age_structured.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using pifsc_bigeye_tuna::BigeyeQuadraObjective;
using pifsc_bigeye_tuna::read_fleet_observations;

namespace {

std::vector<std::string> split_csv(const std::string &line) {
  std::vector<std::string> out;
  std::string cur;
  for (char c : line) {
    if (c == ',') {
      out.push_back(cur);
      cur.clear();
    } else {
      cur.push_back(c);
    }
  }
  out.push_back(cur);
  return out;
}

struct Row {
  std::string name;
  double value = 0.0;
  double gradient = 0.0;
  std::string group;
};

std::vector<Row> read_gradient_rows(const std::string &path) {
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("Could not open gradient CSV: " + path);

  std::string line;
  std::getline(in, line); // header
  std::vector<Row> rows;

  while (std::getline(in, line)) {
    if (line.empty())
      continue;
    const auto x = split_csv(line);
    if (x.size() < 6)
      continue;

    Row r;
    r.name = x[1];
    r.value = std::stod(x[2]);
    r.gradient = std::stod(x[3]);
    r.group = x[5];
    rows.push_back(r);
  }
  return rows;
}

std::vector<double> fixed_vector_from_rows(const std::vector<Row> &rows) {
  std::vector<double> x;
  x.reserve(rows.size());
  for (const auto &r : rows)
    x.push_back(r.value);
  return x;
}

std::vector<double> full_vector_from_fit(const std::vector<double> &fixed,
                                         const Eigen::VectorXd &u_hat) {
  std::vector<double> full = fixed;
  full.reserve(fixed.size() + static_cast<std::size_t>(u_hat.size()));
  for (Eigen::Index i = 0; i < u_hat.size(); ++i)
    full.push_back(u_hat[i]);
  return full;
}

double eval_joint(BigeyeQuadraObjective &obj, const std::vector<double> &fixed,
                  const Eigen::VectorXd &u_hat) {
  auto full = full_vector_from_fit(fixed, u_hat);
  return obj(full);
}

bool in_direction(const Row &r, const std::string &direction) {
  if (direction == "all")
    return true;
  return r.group == direction;
}

} // namespace

int main() {
  const std::string data_path =
      "examples/NMFS/pifsc_bigeye_tuna/"
      "level21_age_based_natural_mortality_diagnostic/data/"
      "synthetic_bigeye_level21_age_based_natural_mortality_observations.csv";

  const std::string grad_path =
      "examples/NMFS/pifsc_bigeye_tuna/"
      "level21_age_based_natural_mortality_diagnostic/outputs/"
      "bigeye_level21_gradient_by_parameter.csv";

  const auto observations = read_fleet_observations(data_path);
  BigeyeQuadraObjective objective(observations);

  auto rows = read_gradient_rows(grad_path);
  if (rows.empty()) {
    throw std::runtime_error("No gradient rows read from " + grad_path);
  }

  auto fixed0 = fixed_vector_from_rows(rows);

  quadra::ParameterVector params;
  for (std::size_t i = 0; i < fixed0.size(); ++i) {
    params.fixed.push_back(static_cast<int>(i));
  }

  const std::size_t n_random = objective.n_years();
  for (std::size_t i = 0; i < n_random; ++i) {
    params.random.push_back(static_cast<int>(fixed0.size() + i));
  }

  Eigen::VectorXd theta(static_cast<Eigen::Index>(fixed0.size()));
  for (std::size_t i = 0; i < fixed0.size(); ++i)
    theta[static_cast<Eigen::Index>(i)] = fixed0[i];

  Eigen::VectorXd u0 =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(n_random));
  auto u_hat =
      quadra::solve_random_effects_laplace(objective, params, theta, u0);

  const double f0 = eval_joint(objective, fixed0, u_hat);

  const std::vector<std::string> directions = {"all",
                                               "age_m",
                                               "base_scale",
                                               "purse_seine_selectivity",
                                               "longline_selectivity",
                                               "initial_numbers"};

  const std::vector<double> alphas = {1e-8, 3e-8, 1e-7, 3e-7, 1e-6, 3e-6,
                                      1e-5, 3e-5, 1e-4, 3e-4, 1e-3};

  const std::string csv_path =
      "examples/NMFS/pifsc_bigeye_tuna/"
      "level21_age_based_natural_mortality_diagnostic/outputs/"
      "bigeye_level21_actual_directional_descent.csv";
  const std::string txt_path =
      "examples/NMFS/pifsc_bigeye_tuna/"
      "level21_age_based_natural_mortality_diagnostic/outputs/"
      "bigeye_level21_actual_directional_descent.txt";

  std::ofstream csv(csv_path);
  std::ofstream txt(txt_path);

  csv << "direction,alpha,n_params,grad_sq,predicted_delta,actual_objective,"
         "actual_delta,"
         "actual_minus_predicted,max_abs_step,max_grad_param,max_abs_grad\n";

  txt << "Bigeye Level 21 Actual Directional Descent Diagnostic\n";
  txt << "=====================================================\n\n";
  txt << std::setprecision(15);
  txt << "Baseline joint objective at fixed fit and re-solved u_hat: " << f0
      << "\n\n";
  txt << "direction,alpha,n_params,grad_sq,predicted_delta,actual_delta,actual_"
         "minus_predicted,max_abs_step,max_grad_param,max_abs_grad\n";

  for (const auto &direction : directions) {
    for (double alpha : alphas) {
      auto trial = fixed0;
      double grad_sq = 0.0;
      double max_abs_step = 0.0;
      double max_abs_grad = 0.0;
      std::string max_grad_param = "";
      int n_params = 0;

      for (std::size_t i = 0; i < rows.size(); ++i) {
        if (!in_direction(rows[i], direction))
          continue;
        const double g = rows[i].gradient;
        const double step = alpha * g;
        trial[i] -= step;

        grad_sq += g * g;
        max_abs_step = std::max(max_abs_step, std::abs(step));
        if (std::abs(g) > max_abs_grad) {
          max_abs_grad = std::abs(g);
          max_grad_param = rows[i].name;
        }
        ++n_params;
      }

      if (n_params == 0)
        continue;

      Eigen::VectorXd trial_theta(static_cast<Eigen::Index>(trial.size()));
      for (std::size_t i = 0; i < trial.size(); ++i) {
        trial_theta[static_cast<Eigen::Index>(i)] = trial[i];
      }

      auto trial_u = quadra::solve_random_effects_laplace(objective, params,
                                                          trial_theta, u_hat);
      const double f_trial = eval_joint(objective, trial, trial_u);
      const double predicted_delta = -alpha * grad_sq;
      const double actual_delta = f_trial - f0;

      csv << direction << "," << alpha << "," << n_params << "," << grad_sq
          << "," << predicted_delta << "," << f_trial << "," << actual_delta
          << "," << (actual_delta - predicted_delta) << "," << max_abs_step
          << "," << max_grad_param << "," << max_abs_grad << "\n";

      txt << direction << "," << alpha << "," << n_params << "," << grad_sq
          << "," << predicted_delta << "," << actual_delta << ","
          << (actual_delta - predicted_delta) << "," << max_abs_step << ","
          << max_grad_param << "," << max_abs_grad << "\n";
    }
  }

  std::cout << "wrote: " << txt_path << "\n";
  std::cout << "wrote: " << csv_path << "\n";
  return 0;
}
