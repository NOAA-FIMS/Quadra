#include <Eigen/Dense>

#include <cmath>
#include <iostream>

struct SimpleModel {
  double operator()(double theta, double u) const {
    return 0.5 * (u - theta) * (u - theta) + 0.25 * u * u;
  }
};

double solve_mode(double theta) { return theta / 1.5; }

double finite_difference_du_dtheta(double theta, double h = 1.0e-6) {
  return (solve_mode(theta + h) - solve_mode(theta - h)) / (2.0 * h);
}

int main() {
  const double theta = 2.5;
  const double uhat = solve_mode(theta);

  SimpleModel model;

  const double h = 1.0e-5;

  const double f_upp = model(theta, uhat + h);
  const double f_u00 = model(theta, uhat);
  const double f_umm = model(theta, uhat - h);

  const double H_uu = (f_upp - 2.0 * f_u00 + f_umm) / (h * h);

  const double f_pp = model(theta + h, uhat + h);
  const double f_pm = model(theta + h, uhat - h);
  const double f_mp = model(theta - h, uhat + h);
  const double f_mm = model(theta - h, uhat - h);

  const double H_u_theta = (f_pp - f_pm - f_mp + f_mm) / (4.0 * h * h);

  const double ift = -H_u_theta / H_uu;

  const double expected = finite_difference_du_dtheta(theta);

  const double error = std::abs(ift - expected);

  if (!std::isfinite(ift)) {
    std::cerr << "FAIL: non-finite IFT derivative\n";
    return 1;
  }

  if (error > 1.0e-6) {
    std::cerr << "FAIL: IFT derivative mismatch\n";
    std::cerr << "IFT: " << ift << "\n";
    std::cerr << "FD : " << expected << "\n";
    std::cerr << "err: " << error << "\n";
    return 1;
  }

  std::cout << "PASS: dense IFT validation\n";
  std::cout << "  uhat: " << uhat << "\n";
  std::cout << "  H_uu: " << H_uu << "\n";
  std::cout << "  H_u_theta: " << H_u_theta << "\n";
  std::cout << "  IFT du/dtheta: " << ift << "\n";
  std::cout << "  FD  du/dtheta: " << expected << "\n";

  return 0;
}
