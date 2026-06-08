#include "../core/laplace/hessian_structure.hpp"
#include "../core/laplace/laplace_backend.hpp"

#include <Eigen/Sparse>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using quadra::laplace::BandedBackend;
using quadra::laplace::LogDetSparseLDLT;

void expect_close(const double a, const double b, const char *msg) {
  const double diff = std::abs(a - b);
  const double scale = 1.0 + std::max(std::abs(a), std::abs(b));

  if (diff > 1e-9 * scale) {
    std::cerr << msg << ": a=" << a << " b=" << b << " diff=" << diff << "\n";
    throw std::runtime_error(msg);
  }
}

Eigen::SparseMatrix<double> make_banded(const int n, const int bandwidth) {
  std::vector<Eigen::Triplet<double>> t;
  t.reserve(static_cast<std::size_t>(n * (2 * bandwidth + 1)));

  for (int i = 0; i < n; ++i) {
    double diag = 10.0 + 0.001 * i;

    for (int d = 1; d <= bandwidth; ++d) {
      if (i - d < 0)
        continue;

      const double e = ((d % 2 == 0) ? 0.015 : -0.025) / static_cast<double>(d);

      t.emplace_back(i, i - d, e);
      t.emplace_back(i - d, i, e);

      diag += 2.0 * std::abs(e);
    }

    t.emplace_back(i, i, diag);
  }

  Eigen::SparseMatrix<double> H(n, n);
  H.setFromTriplets(t.begin(), t.end());
  H.makeCompressed();
  return H;
}

void test_banded_backend(const int n, const int bandwidth) {
  const Eigen::SparseMatrix<double> H = make_banded(n, bandwidth);

  BandedBackend backend(bandwidth);
  backend.factorize(H);

  if (!backend.is_spd()) {
    throw std::runtime_error("banded backend reported non-SPD");
  }

  expect_close(backend.logdet(), LogDetSparseLDLT(H),
               "banded backend structured-value logdet");
}

int main() {
  test_banded_backend(8, 2);
  test_banded_backend(100, 2);
  test_banded_backend(100, 5);
  test_banded_backend(250, 10);

  std::cout << "banded backend structured-value tests passed\n";
  return 0;
}
