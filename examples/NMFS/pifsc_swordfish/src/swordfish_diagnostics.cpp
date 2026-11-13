/**
 * @file swordfish_diagnostics.cpp
 * @details
 * Local fit, curvature and benchmark diagnostics. Perturbed starts use a
 * constant covariance- based linear coordinate change and check the
 * objective/gradient chain rule. Covariance uses central differences of exact
 * first-order AD gradients, not an exact AD Hessian.
 */
#include "core/autodiff/model_gradient.hpp"
#include "core/optimizer.hpp"
#include "swordfish/benchmarks.hpp"
#include "swordfish/model.hpp"
#include "swordfish/reference_io.hpp"
#include "swordfish/report.hpp"
#include "swordfish_input.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
DECLARE_ADGRAPH();
// Constant linear change of coordinates: the likelihood and parameter bounds
// are unchanged.
/// Constant linear coordinate adapter using the fitted covariance Cholesky
/// factor; the likelihood and native bounds are unchanged.
struct RescaledModel : quadra::QuadraModel<RescaledModel> {
  const swordfish::SwordfishModel &model;
  std::vector<double> center;
  Eigen::MatrixXd factor;
  /// Store the underlying model, fitted coordinate center and constant linear
  /// rescaling matrix.
  RescaledModel(const swordfish::SwordfishModel &m, std::vector<double> c,
                Eigen::MatrixXd f)
      : model(m), center(std::move(c)), factor(std::move(f)) {}
  /// Return the underlying model coordinate names in unchanged order.
  std::vector<std::string> parameter_names_impl() const {
    return model.parameter_names_impl();
  }
  /// Map rescaled optimizer coordinates back to original bounded-logit
  /// coordinates as center plus factor times displacement.
  template <class T> std::vector<T> original(const std::vector<T> &z) const {
    std::vector<T> x;
    for (size_t i = 0; i < center.size(); ++i) {
      T value = T(center[i]);
      for (size_t j = 0; j < z.size(); ++j)
        value += T(factor(i, j)) * z[j];
      x.push_back(value);
    }
    return x;
  }
  template <class T> T operator()(const std::vector<T> &z) const {
    return model(original(z));
  }
  /// Evaluate the original joint objective after the constant coordinate
  /// transformation.
  template <class T>
  T evaluate_impl(const std::vector<T> &z, quadra::ModelReportContext &) const {
    return (*this)(z);
  }
};
/// Run swordfish diagnostics; return zero on success and nonzero on a failed
/// check or invalid input.
int main(int argc, char **argv) {
  try {
    auto data = imported_swordfish_data();
    auto native = swordfish::read_fitted_parameters(
        data, "build/quadra_fitted_parameters.csv",
        "build/biology_fitted_parameters.txt",
        "build/population_fitted_deviations.txt");
    std::string line;
    swordfish::SwordfishModel model(data, native.first, native.second);
    auto x = model.initial_vector();
    const std::string mode = argc > 1 ? argv[1] : "curvature";
    if (mode == "derived") {
      const std::vector<std::string> labels = {
          "SSB_unfished", "SSB_Btgt",       "annF_Btgt",      "Dead_Catch_Btgt",
          "SSB_SPR",      "annF_SPR",       "Dead_Catch_SPR", "SSB_MSY",
          "annF_MSY",     "Dead_Catch_MSY", "SSB_2021"};
      auto calculate = [&](const std::vector<double> &point) {
        const auto decoded = model.decode(point);
        swordfish::Benchmarks b(data, decoded.first, decoded.second);
        const auto virgin = b(0), bt = b.target(.2, true),
                   spr = b.target(.2, false), msy = b.msy();
        const auto run =
            swordfish::run_population(data, decoded.first, decoded.second);
        return std::vector<double>{
            virgin.ssb,
            bt.ssb,
            bt.annual_f,
            bt.yield,
            spr.ssb,
            spr.annual_f,
            spr.yield,
            msy.ssb,
            msy.annual_f,
            msy.yield,
            swordfish::population_row(run, data, 2021, 2).ssb};
      };
      const double step = argc > 2 ? swordfish::csv_number(argv[2]) : 1e-4;
      if (!(step > 0 && step < .01))
        throw std::invalid_argument("Derived step outside (0,.01)");
      const auto center = calculate(x);
      std::ofstream estimates("build/derived_estimates.csv");
      estimates << "quantity,estimate\n" << std::setprecision(17);
      for (size_t i = 0; i < labels.size(); ++i)
        estimates << labels[i] << ',' << center[i] << '\n';
      std::ofstream out("build/derived_jacobian.csv");
      out << "parameter";
      for (const auto &label : labels)
        out << ',' << label;
      out << '\n' << std::setprecision(17);
      for (size_t j = 0; j < x.size(); ++j) {
        auto plus = x, minus = x;
        plus[j] += step;
        minus[j] -= step;
        auto a = calculate(plus), b = calculate(minus);
        out << model.free[j].name;
        for (size_t i = 0; i < labels.size(); ++i)
          out << ',' << (a[i] - b[i]) / (2 * step);
        out << '\n';
      }
    } else if (mode == "benchmarks") {
      swordfish::write_reference_points(data, native.first, native.second,
                                        "build");
    } else if (mode == "starts") {
      std::ifstream cov("build/coordinate_covariance.csv");
      if (!std::getline(cov, line))
        throw std::runtime_error("Run uncertainty before robustness");
      Eigen::MatrixXd covariance(x.size(), x.size());
      for (size_t j = 0; j < x.size(); ++j) {
        if (!std::getline(cov, line))
          throw std::runtime_error("Incomplete covariance");
        std::stringstream row(line);
        std::string cell;
        std::getline(row, cell, ',');
        if (cell != model.free[j].name)
          throw std::runtime_error("Covariance order mismatch");
        for (size_t k = 0; k < x.size(); ++k) {
          std::getline(row, cell, ',');
          covariance(j, k) = std::stod(cell);
        }
      }
      Eigen::LLT<Eigen::MatrixXd> cholesky(covariance);
      if (cholesky.info() != Eigen::Success)
        throw std::runtime_error("Covariance factorization failed");
      Eigen::MatrixXd factor = cholesky.matrixL();
      RescaledModel scaled(model, x, factor);
      std::ofstream out("build/robustness.csv");
      out << "seed,logit_sd,initial_objective,final_objective,gradient_norm,"
             "converged,iterations,max_coordinate_distance\n"
          << std::setprecision(16);
      for (unsigned seed : {1701u, 1702u, 1703u}) {
        const double sd = seed == 1701 ? .02 : (seed == 1702 ? .05 : .1);
        std::mt19937 rng(seed);
        std::normal_distribution<double> noise(0, sd);
        Eigen::VectorXd displacement(x.size());
        for (size_t i = 0; i < x.size(); ++i)
          displacement[i] = noise(rng);
        const Eigen::VectorXd z =
            factor.triangularView<Eigen::Lower>().solve(displacement);
        quadra::ParameterVector parameters;
        std::vector<double> start;
        for (size_t i = 0; i < x.size(); ++i) {
          start.push_back(z[i]);
          parameters.add(quadra::Parameter(model.free[i].name, z[i]));
        }
        const double initial = scaled(start);
        const auto direct_gradient =
                       quadra::evaluate_gradient(model, scaled.original(start)),
                   scaled_gradient = quadra::evaluate_gradient(scaled, start);
        if (std::abs(initial - direct_gradient.objective_value_m) > 1e-9)
          throw std::runtime_error("Rescaled objective mismatch");
        for (size_t j = 0; j < x.size(); ++j) {
          double expected = 0;
          for (size_t i = 0; i < x.size(); ++i)
            expected += factor(i, j) * direct_gradient.gradient_m[i];
          if (std::abs(expected - scaled_gradient.gradient_m[j]) >
              1e-7 * (1 + std::abs(expected)))
            throw std::runtime_error("Rescaled gradient chain rule mismatch");
        }
        std::cout << "Fitting seed " << seed << ", logit SD " << sd
                  << ", initial " << initial << std::endl;
        const auto fit = quadra::optimize_lbfgs(
            scaled, parameters, quadra::default_laplace_options(), 500, 1e-7);
        const auto fitted = scaled.original(fit.par);
        const auto gradient = quadra::evaluate_gradient(model, fitted);
        double distance = 0;
        for (size_t i = 0; i < x.size(); ++i)
          distance = std::max(distance, std::abs(fitted[i] - x[i]));
        const bool converged = fit.converged && gradient.gradient_norm_m < 1e-4;
        out << seed << ',' << sd << ',' << initial << ',' << fit.value << ','
            << gradient.gradient_norm_m << ',' << converged << ','
            << fit.iterations << ',' << distance << '\n';
        out.flush();
        std::cout << "Finished " << seed << ": " << fit.value
                  << ", original gradient " << gradient.gradient_norm_m
                  << ", converged " << converged << std::endl;
        if (!converged)
          throw std::runtime_error("Perturbed fit failed convergence");
      }
    } else if (mode == "curvature") {
      // Two central-difference steps of exact AD gradients expose numerical
      // sensitivity.
      for (double step : {1e-4, 5e-5}) {
        const std::string suffix = step == 1e-4 ? "coarse" : "fine";
        std::ofstream out("build/hessian_" + suffix + ".csv");
        out << "parameter";
        for (const auto &p : model.free)
          out << ',' << p.name;
        out << '\n' << std::setprecision(17);
        for (size_t j = 0; j < x.size(); ++j) {
          auto plus = x, minus = x;
          plus[j] += step;
          minus[j] -= step;
          const auto gp = quadra::evaluate_gradient(model, plus),
                     gm = quadra::evaluate_gradient(model, minus);
          out << model.free[j].name;
          for (size_t i = 0; i < x.size(); ++i)
            out << ',' << (gp.gradient_m[i] - gm.gradient_m[i]) / (2 * step);
          out << '\n';
          out.flush();
          if (j % 14 == 0)
            std::cout << suffix << " Hessian column " << j << '/' << x.size()
                      << std::endl;
        }
      }
      std::ofstream out("build/uncertainty_coordinates.csv");
      out << "parameter,coordinate,native,lower,upper\n"
          << std::setprecision(17);
      for (size_t i = 0; i < x.size(); ++i)
        out << model.free[i].name << ',' << x[i] << ',' << model.free[i].value
            << ',' << model.free[i].lo << ',' << model.free[i].hi << '\n';
    } else
      throw std::invalid_argument(
          "Use starts, curvature, derived [step], or benchmarks");
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
