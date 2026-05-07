#ifndef QUADRA_REPORT_HPP
#define QUADRA_REPORT_HPP
#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../eigen/Eigen/Dense"
#include "../autodiff/autodiff.hpp"
#include "../../model/parameter.hpp"
#include "covariance.hpp"

namespace quadra {

struct ReportMetadata {
    std::string units;
    std::string description;
    std::string group;
};

struct ReportEntry {
    std::string name;
    std::string path;
    ReportMetadata metadata;

    double value = std::numeric_limits<double>::quiet_NaN();
    bool estimate_uncertainty = false;

    Eigen::VectorXd gradient;
    double std_error = std::numeric_limits<double>::quiet_NaN();
};

inline std::string report_group_from_path(const std::string& path) {
    const auto pos = path.find_last_of('/');
    if (pos == std::string::npos) return "";
    return path.substr(0, pos);
}

inline std::string report_name_from_path(const std::string& path) {
    const auto pos = path.find_last_of('/');
    if (pos == std::string::npos) return path;
    return path.substr(pos + 1);
}

template <typename Scalar>
class Report {
public:
    std::vector<ReportEntry> entries;

    void clear() { entries.clear(); }

    void add(const std::string& path, const Scalar& value) {
        add_impl(path, value, false, {});
    }

    void estimate(const std::string& path, const Scalar& value) {
        add_impl(path, value, true, {});
    }

    void add(const std::string& path,
             const Scalar& value,
             const ReportMetadata& metadata) {
        add_impl(path, value, false, metadata);
    }

    void estimate(const std::string& path,
                  const Scalar& value,
                  const ReportMetadata& metadata) {
        add_impl(path, value, true, metadata);
    }

    std::size_t size() const { return entries.size(); }

    const ReportEntry& operator[](std::size_t i) const { return entries[i]; }
    ReportEntry& operator[](std::size_t i) { return entries[i]; }

    void print(std::ostream& os = std::cout) const {
        os << "Quadra report\n-------------\n";
        for (const auto& e : entries) {
            os << e.path << " = " << e.value;
            if (std::isfinite(e.std_error)) os << " (SE = " << e.std_error << ")";
            if (!e.metadata.units.empty()) os << " [" << e.metadata.units << "]";
            os << "\n";
        }
    }

    void to_csv(const std::string& filename) const {
        std::ofstream out(filename);
        if (!out) {
            throw std::runtime_error("Report::to_csv: failed to open output file");
        }

        out << "path,group,name,value,std_error,estimate_uncertainty,units,description\n";

        for (const auto& e : entries) {
            out << csv_escape(e.path) << ","
                << csv_escape(e.metadata.group) << ","
                << csv_escape(e.name) << ","
                << e.value << ",";

            if (std::isfinite(e.std_error)) out << e.std_error;
            out << "," << (e.estimate_uncertainty ? "true" : "false") << ","
                << csv_escape(e.metadata.units) << ","
                << csv_escape(e.metadata.description) << "\n";
        }
    }

private:
    static std::string csv_escape(const std::string& s) {
        bool needs_quotes = false;
        for (char c : s) {
            if (c == ',' || c == '"' || c == '\n' || c == '\r') {
                needs_quotes = true;
                break;
            }
        }

        if (!needs_quotes) return s;

        std::string out = "\"";
        for (char c : s) {
            if (c == '"') out += "\"\"";
            else out += c;
        }
        out += "\"";
        return out;
    }

    void add_impl(const std::string& path,
                  const Scalar& value,
                  bool estimate,
                  ReportMetadata metadata) {
        ReportEntry e;
        e.path = path;
        e.name = report_name_from_path(path);

        if (metadata.group.empty()) {
            metadata.group = report_group_from_path(path);
        }

        e.metadata = std::move(metadata);
        e.value = value_of_arithmetic_or_ad(value);
        e.estimate_uncertainty = estimate;
        entries.push_back(std::move(e));
    }
};

template <typename Scalar>
class ReportTape {
public:
    std::vector<std::string> paths;
    std::vector<ReportMetadata> metadata;
    std::vector<Scalar> values;
    std::vector<bool> estimate_uncertainty;

    void add(const std::string& path, const Scalar& value) {
        add_impl(path, value, false, {});
    }

    void estimate(const std::string& path, const Scalar& value) {
        add_impl(path, value, true, {});
    }

    void add(const std::string& path,
             const Scalar& value,
             const ReportMetadata& md) {
        add_impl(path, value, false, md);
    }

    void estimate(const std::string& path,
                  const Scalar& value,
                  const ReportMetadata& md) {
        add_impl(path, value, true, md);
    }

private:
    void add_impl(const std::string& path,
                  const Scalar& value,
                  bool estimate,
                  ReportMetadata md) {
        if (md.group.empty()) md.group = report_group_from_path(path);
        paths.push_back(path);
        metadata.push_back(std::move(md));
        values.push_back(value);
        estimate_uncertainty.push_back(estimate);
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

// One AD pass creates all AD report values. Then each estimated value gets a
// reverse sweep. This avoids re-running model.report(...) for every quantity.
template <typename Model>
Report<double> evaluate_report_with_uncertainty(
    Model& model,
    const ParameterVector& params,
    const CovarianceResult& covariance) {

    had::ADGraph graph;
    ADScope scope(graph);

    std::vector<AD> p;
    p.reserve(static_cast<std::size_t>(params.size()));
    for (int i = 0; i < params.size(); ++i) {
        p.emplace_back(AD(params.params[static_cast<std::size_t>(i)].value));
    }

    ReportTape<AD> tape;
    call_report(model, p, tape);

    Report<double> report;
    report.entries.reserve(tape.values.size());

    for (std::size_t e = 0; e < tape.values.size(); ++e) {
        ReportEntry entry;
        entry.path = tape.paths[e];
        entry.name = report_name_from_path(tape.paths[e]);
        entry.metadata = tape.metadata[e];
        entry.value = value_of_arithmetic_or_ad(tape.values[e]);
        entry.estimate_uncertainty = tape.estimate_uncertainty[e];

        report.entries.push_back(std::move(entry));
    }

    for (std::size_t e = 0; e < tape.values.size(); ++e) {
        if (!tape.estimate_uncertainty[e]) continue;

        // had's reverse sweep mutates adjoints/Hessian storage, so recompute
        // the graph for each report entry. This still avoids re-running the
        // model twice per quantity and keeps the API explicit. A later ADGraph
        // reset/adjoint-only mode can make this cheaper.
        had::ADGraph g_entry;
        ADScope scope_entry(g_entry);

        std::vector<AD> p_entry;
        p_entry.reserve(static_cast<std::size_t>(params.size()));
        for (int i = 0; i < params.size(); ++i) {
            p_entry.emplace_back(AD(params.params[static_cast<std::size_t>(i)].value));
        }

        ReportTape<AD> tape_entry;
        call_report(model, p_entry, tape_entry);

        AD y = tape_entry.values[e];
        scope_entry.backward(y);

        Eigen::VectorXd grad(static_cast<Eigen::Index>(covariance.parameter_indices.size()));

        for (std::size_t i = 0; i < covariance.parameter_indices.size(); ++i) {
            const int pidx = covariance.parameter_indices[i];
            grad[static_cast<Eigen::Index>(i)] =
                scope_entry.grad(p_entry[static_cast<std::size_t>(pidx)]);
        }

        report.entries[e].gradient = grad;

        const double v = (grad.transpose() * covariance.covariance * grad)(0, 0);

        report.entries[e].std_error =
            (v >= 0.0 && std::isfinite(v))
                ? std::sqrt(v)
                : std::numeric_limits<double>::quiet_NaN();
    }

    return report;
}

} // namespace quadra

#endif
