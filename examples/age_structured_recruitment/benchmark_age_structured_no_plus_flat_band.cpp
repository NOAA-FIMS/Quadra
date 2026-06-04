#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using Clock = std::chrono::high_resolution_clock;

namespace {

double ms_between(const Clock::time_point &a, const Clock::time_point &b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

std::vector<int> parse_lengths(const std::string &s) {
  std::vector<int> out;
  std::stringstream ss_in(s);
  std::string item;
  while (std::getline(ss_in, item, ',')) {
    if (!item.empty())
      out.push_back(std::stoi(item));
  }
  return out;
}

double logistic(const double x) { return 1.0 / (1.0 + std::exp(-x)); }

double normal_nll(const double residual, const double sigma) {
  const double z = residual / sigma;
  return std::log(sigma) + 0.5 * std::log(2.0 * M_PI) + 0.5 * z * z;
}

struct Data {
  int n_years = 0;
  int n_ages = 0;
  std::vector<double> index_obs;
};

struct Parameters {
  double log_R0 = std::log(1000.0);
  double log_M = std::log(0.20);
  double log_q = std::log(0.001);
  double log_sigma_R = std::log(0.35);
  double log_sigma_index = std::log(0.10);
};

struct Derived {
  double R0 = 0.0;
  double M = 0.0;
  double survival = 0.0;
  double q = 0.0;
  double sigma_R = 0.0;
  double sigma_index = 0.0;
  double inv_var_R = 0.0;
  double inv_var_index = 0.0;
};

Derived transform(const Parameters &par) {
  Derived d;
  d.R0 = std::exp(par.log_R0);
  d.M = std::exp(par.log_M);
  d.survival = std::exp(-d.M);
  d.q = std::exp(par.log_q);
  d.sigma_R = std::exp(par.log_sigma_R);
  d.sigma_index = std::exp(par.log_sigma_index);
  d.inv_var_R = 1.0 / (d.sigma_R * d.sigma_R);
  d.inv_var_index = 1.0 / (d.sigma_index * d.sigma_index);
  return d;
}

std::vector<double> selectivity(const int n_ages) {
  std::vector<double> sel(static_cast<std::size_t>(n_ages));
  for (int a = 0; a < n_ages; ++a) {
    sel[static_cast<std::size_t>(a)] =
        logistic((static_cast<double>(a + 1) - 4.0) / 0.8);
  }
  return sel;
}

std::vector<double> simulate_index(const int n_years, const int n_ages,
                                   const Parameters &par) {
  const Derived p = transform(par);
  const std::vector<double> sel = selectivity(n_ages);

  std::vector<double> N(static_cast<std::size_t>(n_ages), 0.0);
  for (int a = 0; a < n_ages; ++a) {
    N[static_cast<std::size_t>(a)] =
        p.R0 * std::exp(-p.M * static_cast<double>(a));
  }

  std::vector<double> index(static_cast<std::size_t>(n_years), 0.0);

  for (int y = 0; y < n_years; ++y) {
    double vulnerable = 0.0;
    for (int a = 0; a < n_ages; ++a) {
      vulnerable +=
          sel[static_cast<std::size_t>(a)] * N[static_cast<std::size_t>(a)];
    }

    const double obs_error =
        0.04 * std::sin(2.0 * M_PI * static_cast<double>(y) / 13.0) +
        0.02 * std::cos(2.0 * M_PI * static_cast<double>(y) / 7.0);

    index[static_cast<std::size_t>(y)] = p.q * vulnerable * std::exp(obs_error);

    std::vector<double> N_next(static_cast<std::size_t>(n_ages), 0.0);
    const double recruitment_dev =
        0.15 * std::sin(2.0 * M_PI * static_cast<double>(y) / 17.0);
    N_next[0] = p.R0 * std::exp(recruitment_dev);

    for (int a = 1; a < n_ages; ++a) {
      N_next[static_cast<std::size_t>(a)] =
          N[static_cast<std::size_t>(a - 1)] * p.survival;
    }

    N.swap(N_next);
  }

  return index;
}

Data make_data(const int n_years, const int n_ages, const Parameters &par) {
  Data data;
  data.n_years = n_years;
  data.n_ages = n_ages;
  data.index_obs = simulate_index(n_years, n_ages, par);
  return data;
}

struct EvalAll {
  double objective = 0.0;
  Eigen::VectorXd gradient;
  Eigen::SparseMatrix<double> hessian;
};

struct ActiveCohort {
  int k = -1;         // recruitment deviation index
  int age_index = -1; // zero-based age slot currently occupied
  double N = 0.0;     // abundance contribution
  double dN = 0.0;    // d abundance / dx_k
  double ddN = 0.0;   // d2 abundance / dx_k^2
};

struct PairKey {
  int i;
  int j;
};

EvalAll eval_all_flat_band(const Data &data, const Parameters &par,
                           const Eigen::VectorXd &x) {
  const Derived p = transform(par);
  const std::vector<double> sel = selectivity(data.n_ages);

  const int n = data.n_years;
  const int A = data.n_ages;

  EvalAll out;
  out.gradient = Eigen::VectorXd::Zero(n);

  double nll = 0.0;

  // Initial equilibrium-like numbers are fixed, not dependent on x.
  std::vector<double> fixed_N(static_cast<std::size_t>(A), 0.0);
  for (int a = 0; a < A; ++a) {
    fixed_N[static_cast<std::size_t>(a)] =
        p.R0 * std::exp(-p.M * static_cast<double>(a));
  }

  std::vector<ActiveCohort> active;
  active.reserve(static_cast<std::size_t>(A));

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(static_cast<std::size_t>(n * (2 * A + 1)));

  // Direct flat band storage. For this no-plus model, max distance is A.
  const int bw = A;
  std::vector<double> Hband(static_cast<std::size_t>((2 * bw + 1) * n), 0.0);

  auto add_H = [&](int i, int j, double value) {
    if (i < 0 || j < 0 || i >= n || j >= n)
      return;
    const int d = j - i;
    if (std::abs(d) > bw) {
      throw std::runtime_error("flat-band Hessian write outside bandwidth");
    }
    const int band = d + bw;
    Hband[static_cast<std::size_t>(band * n + i)] += value;
  };

  for (int y = 0; y < n; ++y) {
    // Observation before current-year recruitment enters.
    double vulnerable = 0.0;

    for (int a = 0; a < A; ++a) {
      vulnerable += sel[static_cast<std::size_t>(a)] *
                    fixed_N[static_cast<std::size_t>(a)];
    }

    for (const auto &c : active) {
      vulnerable += sel[static_cast<std::size_t>(c.age_index)] * c.N;
    }

    std::vector<int> keys;
    std::vector<double> dV;
    std::vector<double> ddV;

    keys.reserve(active.size());
    dV.reserve(active.size());
    ddV.reserve(active.size());

    for (const auto &c : active) {
      keys.push_back(c.k);
      const double s = sel[static_cast<std::size_t>(c.age_index)];
      dV.push_back(s * c.dN);
      ddV.push_back(s * c.ddN);
    }

    const double obs_resid =
        std::log(data.index_obs[static_cast<std::size_t>(y)]) -
        (std::log(p.q) + std::log(vulnerable));

    nll += normal_nll(obs_resid, p.sigma_index);

    for (std::size_t ii = 0; ii < keys.size(); ++ii) {
      const int i = keys[ii];
      const double de_i = -dV[ii] / vulnerable;
      out.gradient[i] += obs_resid * de_i * p.inv_var_index;

      for (std::size_t jj = ii; jj < keys.size(); ++jj) {
        const int j = keys[jj];
        const double de_j = -dV[jj] / vulnerable;

        const double Vij = (i == j) ? ddV[ii] : 0.0;
        const double d2e =
            -(Vij / vulnerable - (dV[ii] * dV[jj]) / (vulnerable * vulnerable));

        const double hij = (de_i * de_j + obs_resid * d2e) * p.inv_var_index;
        add_H(i, j, hij);
      }
    }

    // Recruitment prior for current x[y].
    nll += normal_nll(x[y], p.sigma_R);
    out.gradient[y] += x[y] * p.inv_var_R;
    add_H(y, y, p.inv_var_R);

    // Age fixed initial population one year. No plus group.
    std::vector<double> fixed_next(static_cast<std::size_t>(A), 0.0);
    fixed_next[0] = 0.0; // current recruitment handled as active cohort
    for (int a = 1; a < A; ++a) {
      fixed_next[static_cast<std::size_t>(a)] =
          fixed_N[static_cast<std::size_t>(a - 1)] * p.survival;
    }
    fixed_N.swap(fixed_next);

    // Age active cohorts, dropping those beyond terminal age.
    std::vector<ActiveCohort> next_active;
    next_active.reserve(static_cast<std::size_t>(A));

    for (const auto &c : active) {
      if (c.age_index + 1 < A) {
        ActiveCohort nc = c;
        nc.age_index += 1;
        nc.N *= p.survival;
        nc.dN *= p.survival;
        nc.ddN *= p.survival;
        next_active.push_back(nc);
      }
    }

    // Add current recruitment as age 0 for future observations.
    const double R_y = p.R0 * std::exp(x[y]);
    next_active.push_back(ActiveCohort{y, 0, R_y, R_y, R_y});

    active.swap(next_active);
  }

  out.objective = nll;

  for (int i = 0; i < n; ++i) {
    for (int d = -bw; d <= bw; ++d) {
      const int j = i + d;
      if (j < 0 || j >= n)
        continue;

      const double v = Hband[static_cast<std::size_t>((d + bw) * n + i)];
      if (std::abs(v) > 1e-12) {
        triplets.emplace_back(i, j, v);
        if (i != j) {
          triplets.emplace_back(j, i, v);
        }
      }
    }
  }

  out.hessian.resize(n, n);
  out.hessian.setFromTriplets(triplets.begin(), triplets.end());

  return out;
}

class ObjX {
public:
  ObjX(const Data &data, const Parameters &par) : data_(data), par_(par) {}

  double operator()(const Eigen::VectorXd &x, Eigen::VectorXd &grad) {
    const EvalAll e = eval_all_flat_band(data_, par_, x);
    grad = e.gradient;
    return e.objective;
  }

private:
  const Data &data_;
  const Parameters &par_;
};

Eigen::VectorXd optimize_x(const Data &data, const Parameters &par) {
  Eigen::VectorXd x = Eigen::VectorXd::Zero(data.n_years);

  double f = eval_all_flat_band(data, par, x).objective;

  for (int iter = 0; iter < 50; ++iter) {
    const EvalAll e = eval_all_flat_band(data, par, x);
    f = e.objective;

    const double gnorm = e.gradient.norm();
    if (gnorm < 1e-6)
      break;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
    ldlt.compute(e.hessian);
    if (ldlt.info() != Eigen::Success) {
      throw std::runtime_error("Newton Hessian factorization failed");
    }

    const Eigen::VectorXd dx = ldlt.solve(e.gradient);
    const double step_norm = dx.norm();
    if (step_norm < 1e-10 * (1.0 + x.norm())) {
      break;
    }
    if (ldlt.info() != Eigen::Success || !dx.allFinite()) {
      throw std::runtime_error("Newton Hessian solve failed");
    }

    double step = 1.0;
    bool accepted = false;

    for (int ls = 0; ls < 30; ++ls) {
      const Eigen::VectorXd candidate = x - step * dx;
      const double f_candidate =
          eval_all_flat_band(data, par, candidate).objective;

      const double rel_tol = 1e-12 * (1.0 + std::abs(f));
      if (std::isfinite(f_candidate) &&
          (f_candidate < f ||
           (std::abs(f_candidate - f) <= rel_tol && gnorm < 1e-6))) {
        x = candidate;
        f = f_candidate;
        accepted = true;
        break;
      }

      step *= 0.5;
    }

    if (!accepted) {
      if (gnorm < 1e-5)
        break;
      throw std::runtime_error("Newton line search failed");
    }
  }

  return x;
}

double sparse_logdet_ldlt(const Eigen::SparseMatrix<double> &H) {
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
  ldlt.compute(H);

  if (ldlt.info() != Eigen::Success) {
    throw std::runtime_error("sparse LDLT failed");
  }

  const auto &D = ldlt.vectorD();
  double logdet = 0.0;

  for (int i = 0; i < D.size(); ++i) {
    if (!(D[i] > 0.0)) {
      throw std::runtime_error("sparse Hxx is not positive definite");
    }
    logdet += std::log(D[i]);
  }

  return logdet;
}

struct EvalResult {
  double marginal = 0.0;
  double joint = 0.0;
  double logdet = 0.0;
  double grad_norm = 0.0;
  int nnz = 0;
};

EvalResult eval_laplace(const Data &data, const Parameters &par) {
  const Eigen::VectorXd xhat = optimize_x(data, par);
  const EvalAll e = eval_all_flat_band(data, par, xhat);

  EvalResult out;
  out.joint = e.objective;
  out.grad_norm = e.gradient.norm();
  out.nnz = static_cast<int>(e.hessian.nonZeros());
  out.logdet = sparse_logdet_ldlt(e.hessian);

  const double n_x = static_cast<double>(xhat.size());
  out.marginal =
      out.joint + 0.5 * out.logdet - 0.5 * n_x * std::log(2.0 * M_PI);

  return out;
}

} // namespace

int main(int argc, char **argv) {
  int reps = 10;
  std::vector<int> lengths = {25, 50, 100, 250, 500, 1000};
  int n_ages = 10;

  if (argc > 1)
    reps = std::stoi(argv[1]);
  if (argc > 2)
    lengths = parse_lengths(argv[2]);
  if (argc > 3)
    n_ages = std::stoi(argv[3]);

  const Parameters par;

  std::cout << std::fixed << std::setprecision(6);
  std::cout
      << "Quadra no-plus age-structured flat-band analytic Laplace benchmark\n";
  std::cout << "reps per n = " << reps << ", ages = " << n_ages << "\n\n";

  std::cout << std::setw(8) << "n" << std::setw(14) << "objective"
            << std::setw(14) << "joint" << std::setw(14) << "logdet"
            << std::setw(14) << "nnz" << std::setw(14) << "grad_norm"
            << std::setw(14) << "avg_ms" << "\n";

  for (const int n : lengths) {
    const Data data = make_data(n, n_ages, par);
    EvalResult last = eval_laplace(data, par);

    const auto t0 = Clock::now();
    for (int r = 0; r < reps; ++r) {
      last = eval_laplace(data, par);
    }
    const auto t1 = Clock::now();

    const double avg_ms = ms_between(t0, t1) / static_cast<double>(reps);

    std::cout << std::setw(8) << n << std::setw(14) << last.marginal
              << std::setw(14) << last.joint << std::setw(14) << last.logdet
              << std::setw(14) << last.nnz << std::setw(14) << last.grad_norm
              << std::setw(14) << avg_ms << "\n";
  }

  return 0;
}
