#ifndef QUADRA_REPORT_HPP
#define QUADRA_REPORT_HPP
#pragma once

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "../eigen/Eigen/Dense"
#include "../autodiff/autodiff.hpp"
#include "../../model/parameter.hpp"
#include "covariance.hpp"

namespace quadra {

struct ReportEntry {
    std::string name;
    double value = std::numeric_limits<double>::quiet_NaN();
    bool estimate_uncertainty = false;
    Eigen::VectorXd gradient;
    double std_error = std::numeric_limits<double>::quiet_NaN();
};

template <typename Scalar>
class Report {
public:
    std::vector<ReportEntry> entries;

    void add(const std::string& name, const Scalar& value) {
        add_impl(name, value, false);
    }

    void estimate(const std::string& name, const Scalar& value) {
        add_impl(name, value, true);
    }

    std::size_t size() const { return entries.size(); }

    void print(std::ostream& os = std::cout) const {
        os << "Quadra report\n-------------\n";
        for (const auto& e : entries) {
            os << e.name << " = " << e.value;
            if (std::isfinite(e.std_error)) os << " (SE = " << e.std_error << ")";
            os << "\n";
        }
    }

private:
    void add_impl(const std::string& name, const Scalar& value, bool estimate) {
        ReportEntry e;
        e.name = name;
        e.value = value_of_arithmetic_or_ad(value);
        e.estimate_uncertainty = estimate;
        entries.push_back(std::move(e));
    }
};

template <typename Scalar>
class ReportTape {
public:
    std::vector<std::string> names;
    std::vector<Scalar> values;
    std::vector<bool> estimate_uncertainty;

    void add(const std::string& name, const Scalar& value) {
        names.push_back(name);
        values.push_back(value);
        estimate_uncertainty.push_back(false);
    }

    void estimate(const std::string& name, const Scalar& value) {
        names.push_back(name);
        values.push_back(value);
        estimate_uncertainty.push_back(true);
    }
};

template <typename Model, typename Scalar, typename ReportLike>
auto call_report_impl(Model& model, const std::vector<Scalar>& p, ReportLike& report, int)
    -> decltype(model.report(p, report), void()) {
    model.report(p, report);
}

template <typename Model, typename Scalar, typename ReportLike>
void call_report_impl(Model&, const std::vector<Scalar>&, ReportLike&, long) {
    throw std::runtime_error("Model does not define report(p, report)");
}

template <typename Model, typename Scalar, typename ReportLike>
void call_report(Model& model, const std::vector<Scalar>& p, ReportLike& report) {
    call_report_impl(model, p, report, 0);
}

template <typename Scalar>
std::vector<Scalar> full_parameter_vector(const ParameterVector& params) {
    std::vector<Scalar> p;
    p.reserve(static_cast<std::size_t>(params.size()));
    for (int i = 0; i < params.size(); ++i) {
        p.emplace_back(params.params[static_cast<std::size_t>(i)].value);
    }
    return p;
}

template <typename Model>
Report<double> evaluate_report_values(Model& model, const ParameterVector& params) {
    auto p = full_parameter_vector<double>(params);
    Report<double> report;
    call_report(model, p, report);
    return report;
}

template <typename Model>
Eigen::VectorXd evaluate_single_report_gradient(
    Model& model,
    const ParameterVector& params,
    std::size_t entry_index,
    const std::vector<int>& covariance_parameter_indices) {

    had::ADGraph graph;
    ADScope scope(graph);

    std::vector<AD> p;
    p.reserve(static_cast<std::size_t>(params.size()));
    for (int i = 0; i < params.size(); ++i) {
        p.emplace_back(AD(params.params[static_cast<std::size_t>(i)].value));
    }

    ReportTape<AD> tape;
    call_report(model, p, tape);

    if (entry_index >= tape.values.size()) {
        throw std::out_of_range("evaluate_single_report_gradient: entry index out of range");
    }

    AD y = tape.values[entry_index];
    scope.backward(y);

    Eigen::VectorXd g(static_cast<Eigen::Index>(covariance_parameter_indices.size()));
    for (std::size_t i = 0; i < covariance_parameter_indices.size(); ++i) {
        const int pidx = covariance_parameter_indices[i];
        g[static_cast<Eigen::Index>(i)] = scope.grad(p[static_cast<std::size_t>(pidx)]);
    }

    return g;
}

template <typename Model>
Report<double> evaluate_report_with_uncertainty(
    Model& model,
    const ParameterVector& params,
    const CovarianceResult& covariance) {

    auto report = evaluate_report_values(model, params);

    for (std::size_t e = 0; e < report.entries.size(); ++e) {
        if (!report.entries[e].estimate_uncertainty) continue;

        Eigen::VectorXd g = evaluate_single_report_gradient(
            model, params, e, covariance.parameter_indices);

        report.entries[e].gradient = g;

        const double v = (g.transpose() * covariance.covariance * g)(0, 0);
        report.entries[e].std_error =
            (v >= 0.0 && std::isfinite(v))
                ? std::sqrt(v)
                : std::numeric_limits<double>::quiet_NaN();
    }

    return report;
}

} // namespace quadra
#endif
