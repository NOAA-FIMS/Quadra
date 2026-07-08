#include "../objective/bigeye_quadra_objective.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct GradRow {
  int rank = 0;
  std::string name;
  double value = 0.0;
  double grad = 0.0;
  double abs_grad = 0.0;
};

std::vector<std::string> split_csv_line(const std::string &line) {
  std::vector<std::string> out;
  std::string cur;
  bool in_quotes = false;
  for (char c : line) {
    if (c == '"') {
      in_quotes = !in_quotes;
    } else if (c == ',' && !in_quotes) {
      out.push_back(cur);
      cur.clear();
    } else {
      cur.push_back(c);
    }
  }
  out.push_back(cur);
  return out;
}

std::vector<GradRow> read_gradient_csv(const std::string &path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("could not open gradient CSV: " + path);

  std::vector<GradRow> rows;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    auto v = split_csv_line(line);
    if (v.size() < 5 || v[0] == "rank") continue;
    GradRow r;
    r.rank = std::stoi(v[0]);
    r.name = v[1];
    r.value = std::stod(v[2]);
    r.grad = std::stod(v[3]);
    r.abs_grad = std::stod(v[4]);
    rows.push_back(r);
  }
  std::sort(rows.begin(), rows.end(),
            [](const GradRow &a, const GradRow &b) { return a.rank < b.rank; });
  return rows;
}

bool is_age_m(const std::string &name) {
  return name == "log_m_young_offset" || name == "log_m_old_offset";
}

bool is_base_scale(const std::string &name) {
  return name == "log_r0" || name == "log_fbar" || name == "log_q_purse_seine";
}

std::vector<double> values_from_rows(const std::vector<GradRow> &rows) {
  std::vector<double> x(rows.size());
  for (std::size_t i = 0; i < rows.size(); ++i) x[i] = rows[i].value;
  return x;
}

std::vector<double> direction_from_rows(const std::vector<GradRow> &rows,
                                        const std::string &direction) {
  std::vector<double> g(rows.size(), 0.0);
  for (std::size_t i = 0; i < rows.size(); ++i) {
    const auto &name = rows[i].name;
    if (direction == "full") {
      g[i] = rows[i].grad;
    } else if (direction == "age_m" && is_age_m(name)) {
      g[i] = rows[i].grad;
    } else if (direction == "base_scale" && is_base_scale(name)) {
      g[i] = rows[i].grad;
    }
  }
  return g;
}

double norm2(const std::vector<double> &g) {
  double s = 0.0;
  for (double x : g) s += x * x;
  return std::sqrt(s);
}

std::vector<double> step(const std::vector<double> &x,
                         const std::vector<double> &g,
                         double alpha) {
  std::vector<double> y = x;
  for (std::size_t i = 0; i < y.size(); ++i) y[i] -= alpha * g[i];
  return y;
}

} // namespace

int main() {
  using namespace pifsc_bigeye_tuna;

  const std::string out_dir =
      "examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs";
  const std::string wf_dir =
      "examples/NMFS/pifsc_bigeye_tuna/workflow";
  const std::string grad_csv = out_dir + "/bigeye_level21_gradient_by_parameter.csv";
  const std::string report_csv = wf_dir + "/bigeye_level21_actual_directional_descent.csv";
  const std::string report_txt = wf_dir + "/bigeye_level21_actual_directional_descent.txt";

  const auto rows = read_gradient_csv(grad_csv);
  if (rows.empty()) throw std::runtime_error("no gradient rows parsed");

  // API-sensitive line:
  // If this constructor is wrong, mirror the object construction from the normal driver.
  BigeyeQuadraObjective<double> objective;

  const std::vector<double> theta = values_from_rows(rows);
  const double f0 = objective(theta);

  const std::vector<double> alphas = {0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0};
  const std::vector<std::string> directions = {"full", "age_m", "base_scale"};

  std::ofstream csv(report_csv);
  csv << std::setprecision(15);
  csv << "direction,alpha,f0,f_alpha,delta_f,gradient_norm,expected_first_order_delta\n";

  std::ofstream txt(report_txt);
  txt << std::setprecision(15);
  txt << "Bigeye Level 21 Actual Directional Descent Check\n";
  txt << "===============================================\n\n";
  txt << "f0," << f0 << "\n\n";
  txt << "direction,alpha,f_alpha,delta_f,gradient_norm,expected_first_order_delta\n";

  for (const auto &direction : directions) {
    const auto g = direction_from_rows(rows, direction);
    const double gn = norm2(g);
    const double g2 = gn * gn;

    for (const double alpha : alphas) {
      const auto candidate = step(theta, g, alpha);
      const double fa = objective(candidate);
      const double delta = fa - f0;
      const double expected = -alpha * g2;

      csv << direction << "," << alpha << "," << f0 << "," << fa << ","
          << delta << "," << gn << "," << expected << "\n";

      txt << direction << "," << alpha << "," << fa << "," << delta << ","
          << gn << "," << expected << "\n";
    }
  }

  std::cout << "wrote: " << report_txt << "\n";
  std::cout << "wrote: " << report_csv << "\n";
  return 0;
}
