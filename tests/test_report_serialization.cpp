#include "../core/inference/report_serialization.hpp"
#include "../core/model/parameter_transform.hpp"

#include <Eigen/Dense>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

int main() {
  quadra::LaplaceProfiledDerivedReport report;
  report.success_m = true;
  report.message_m = "ok";
  report.names_m = {"depletion", "ssb"};

  report.delta_m.success_m = true;
  report.delta_m.estimate_m = Eigen::VectorXd(2);
  report.delta_m.std_error_m = Eigen::VectorXd(2);
  report.delta_m.cv_m = Eigen::VectorXd(2);

  report.delta_m.estimate_m << 0.55, 1234.0;
  report.delta_m.std_error_m << 0.04, 100.0;
  report.delta_m.cv_m << 0.0727272727, 0.0810372771;

  report.delta_m.covariance_m = Eigen::MatrixXd(2, 2);
  report.delta_m.correlation_m = Eigen::MatrixXd(2, 2);

  report.delta_m.covariance_m << 0.0016, 2.0, 2.0, 10000.0;

  report.delta_m.correlation_m << 1.0, 0.5, 0.5, 1.0;

  const std::string report_path =
      "tests/profiled_report_serialization_test.csv";

  const std::string matrix_path = "tests/profiled_report_covariance_test.csv";

  quadra::write_profiled_derived_report_csv(report_path, report);

  quadra::write_matrix_csv(matrix_path, report.delta_m.covariance_m);

  std::ifstream report_in(report_path);

  if (!report_in) {
    std::cerr << "FAIL: could not read report CSV\\n";
    return 1;
  }

  std::vector<std::string> lines;
  std::string line;

  while (std::getline(report_in, line)) {
    lines.push_back(line);
  }

  if (lines.size() != 3) {
    std::cerr << "FAIL: wrong number of report lines: " << lines.size()
              << "\\n";
    return 1;
  }

  if (lines[0] != "quantity,estimate,std_error,cv") {
    std::cerr << "FAIL: bad report header: " << lines[0] << "\\n";
    return 1;
  }

  if (lines[1].find("depletion,") != 0 || lines[2].find("ssb,") != 0) {
    std::cerr << "FAIL: bad report rows\\n";
    return 1;
  }

  std::ifstream matrix_in(matrix_path);

  if (!matrix_in) {
    std::cerr << "FAIL: could not read matrix CSV\\n";
    return 1;
  }

  std::vector<std::string> matrix_lines;

  while (std::getline(matrix_in, line)) {
    matrix_lines.push_back(line);
  }

  if (matrix_lines.size() != 2) {
    std::cerr << "FAIL: wrong number of matrix rows\\n";
    return 1;
  }

  if (matrix_lines[0].find(",") == std::string::npos) {
    std::cerr << "FAIL: matrix row does not contain comma\\n";
    return 1;
  }

  std::cout << "PASS: report serialization\\n";
  std::cout << "  report: " << report_path << "\\n";
  std::cout << "  matrix: " << matrix_path << "\\n";

  return 0;
}
