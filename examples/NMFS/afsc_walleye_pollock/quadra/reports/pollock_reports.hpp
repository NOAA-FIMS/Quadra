#pragma once

#include "../../../../../core/diagnostics/functional_analysis.hpp"

#include <string>

namespace pollock_example {

inline void
write_pollock_markdown_report(const std::string &md_path,
                              const std::string &functional_csv_path,
                              const std::string &structure_txt_path) {
  quadra::diagnostics::MarkdownReportConfig config;
  config.title = "Synthetic Walleye Pollock Functional Analysis";
  config.subtitle =
      "Synthetic and public-data-safe. Not an official assessment.";
  config.output_path = md_path;
  config.functional_csv_path = functional_csv_path;
  config.structure_txt_path = structure_txt_path;
  config.fixed_effects = "2";
  config.total_estimated = "22";
  config.effective_entries_95 = "58";
  config.effective_bandwidth_95 = "1";

  quadra::diagnostics::write_markdown_report(config);
}

} // namespace pollock_example
