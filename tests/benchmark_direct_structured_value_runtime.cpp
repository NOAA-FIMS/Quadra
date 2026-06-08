#include "../core/laplace/persistent_structured_runtime.hpp"
#include "../core/laplace/structured_value_backend.hpp"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using quadra::laplace::BandedValues;
using quadra::laplace::DiagonalValues;
using quadra::laplace::PersistentStructuredRuntimeState;
using quadra::laplace::TridiagonalValues;

using Clock = std::chrono::high_resolution_clock;

double ms_between(const Clock::time_point &a, const Clock::time_point &b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

std::vector<int> parse_lengths(const std::string &s) {
  std::vector<int> out;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, ',')) {
    if (!item.empty())
      out.push_back(std::stoi(item));
  }
  return out;
}

DiagonalValues make_diag(const int n, const double scale) {
  DiagonalValues H;
  H.diag = Eigen::VectorXd::Zero(n);
  for (int i = 0; i < n; ++i) {
    H.diag[i] = scale * (2.0 + 0.001 * i);
  }
  return H;
}

TridiagonalValues make_tri(const int n, const double scale) {
  TridiagonalValues H;
  H.diag = Eigen::VectorXd::Zero(n);
  H.offdiag = Eigen::VectorXd::Zero(std::max(0, n - 1));
  for (int i = 0; i < n; ++i) {
    H.diag[i] = scale * (4.0 + 0.001 * i);
    if (i > 0) {
      H.offdiag[i - 1] = scale * (-0.20 + 0.00001 * (i % 17));
    }
  }
  return H;
}

BandedValues make_banded(const int n, const int bandwidth, const double scale) {
  BandedValues H;
  H.bandwidth = bandwidth;
  H.diag = Eigen::VectorXd::Zero(n);
  H.lower_bands.resize(static_cast<std::size_t>(bandwidth));

  for (int d = 1; d <= bandwidth; ++d) {
    H.lower_bands[static_cast<std::size_t>(d - 1)] =
        Eigen::VectorXd::Zero(std::max(0, n - d));
  }

  for (int i = 0; i < n; ++i) {
    double diag = scale * (10.0 + 0.001 * i);
    for (int d = 1; d <= bandwidth; ++d) {
      const int j = i - d;
      if (j < 0)
        continue;
      const double e =
          scale * (((d % 2 == 0) ? 0.015 : -0.025) / static_cast<double>(d));
      H.lower_bands[static_cast<std::size_t>(d - 1)][j] = e;
      diag += 2.0 * std::abs(e);
    }
    H.diag[i] = diag;
  }
  return H;
}

struct Result {
  std::string name;
  int n = 0;
  int band = 0;
  double construct_ms = 0.0;
  double direct_update_ms = 0.0;
  double logdet = 0.0;
};

template <class Maker>
Result bench_case(const std::string &name, const int n, const int band,
                  const int reps, Maker maker) {
  volatile double acc = 0.0;

  const auto c0 = Clock::now();
  for (int r = 0; r < reps; ++r) {
    const double scale = 1.0 + 1e-6 * static_cast<double>(r + 1);
    const auto H = maker(n, band, scale);
    acc += static_cast<double>(H.diag.size());
  }
  const auto c1 = Clock::now();

  PersistentStructuredRuntimeState state;
  const auto H0 = maker(n, band, 1.0);
  state.update_direct(H0);

  const auto u0 = Clock::now();
  for (int r = 0; r < reps; ++r) {
    const double scale = 1.0 + 1e-6 * static_cast<double>(r + 1);
    const auto H = maker(n, band, scale);
    state.update_direct(H);
    acc += state.logdet();
  }
  const auto u1 = Clock::now();

  (void)acc;

  Result out;
  out.name = name;
  out.n = n;
  out.band = band;
  out.construct_ms = ms_between(c0, c1) / static_cast<double>(reps);
  out.direct_update_ms = ms_between(u0, u1) / static_cast<double>(reps);
  out.logdet = state.logdet();
  return out;
}

int main(int argc, char **argv) {
  int reps = 20;
  std::vector<int> lengths = {100, 500, 1000, 5000};

  if (argc > 1)
    reps = std::stoi(argv[1]);
  if (argc > 2)
    lengths = parse_lengths(argv[2]);

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "Direct structured value runtime benchmark\n";
  std::cout << "reps per case = " << reps << "\n\n";

  std::cout << std::setw(14) << "case" << std::setw(8) << "n" << std::setw(8)
            << "band" << std::setw(16) << "construct_ms" << std::setw(18)
            << "direct_update_ms" << std::setw(14) << "logdet" << "\n";

  for (const int n : lengths) {
    const Result d =
        bench_case("diagonal", n, 0, reps, [](int n_, int, double scale) {
          return make_diag(n_, scale);
        });
    const Result t =
        bench_case("tridiagonal", n, 1, reps, [](int n_, int, double scale) {
          return make_tri(n_, scale);
        });
    const Result b2 =
        bench_case("banded2", n, 2, reps, [](int n_, int band_, double scale) {
          return make_banded(n_, band_, scale);
        });
    const Result b5 =
        bench_case("banded5", n, 5, reps, [](int n_, int band_, double scale) {
          return make_banded(n_, band_, scale);
        });
    const Result b10 = bench_case("banded10", n, 10, reps,
                                  [](int n_, int band_, double scale) {
                                    return make_banded(n_, band_, scale);
                                  });

    for (const auto &r : {d, t, b2, b5, b10}) {
      std::cout << std::setw(14) << r.name << std::setw(8) << r.n
                << std::setw(8) << r.band << std::setw(16) << r.construct_ms
                << std::setw(18) << r.direct_update_ms << std::setw(14)
                << r.logdet << "\n";
    }
  }

  return 0;
}
