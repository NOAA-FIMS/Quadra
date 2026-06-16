#include "red_snapper_age_structured.hpp"

#include "../../../../core/optimizer.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace sefsc_red_snapper
{

  template <class T>
  T exp_t(const T &x)
  {
    using std::exp;
    return exp(x);
  }

  template <class T>
  T log_t(const T &x)
  {
    using std::log;
    return log(x);
  }

  template <class T>
  T invlogit_t(const T &x)
  {
    return T(1.0) / (T(1.0) + exp_t(-x));
  }

  template <class T>
  T max_t(const T &x, double floor)
  {
    return x > T(floor) ? x : T(floor);
  }

  template <class T>
  T square_t(const T &x)
  {
    return x * x;
  }

  template <class T>
  T logistic_selectivity_t(const T &age, const T &a50, const T &slope)
  {
    return T(1.0) / (T(1.0) + exp_t(-slope * (age - a50)));
  }

  template <class T>
  T age_comp_nll(const std::array<double, kAges> &observed,
                 const std::array<T, kAges> &predicted,
                 double effective_n,
                 double floor = 1.0e-12)
  {
    T nll = T(0.0);
    for (int a = 0; a < kAges; ++a)
    {
      const auto i = static_cast<std::size_t>(a);
      const double obs = std::max(observed[i], 0.0);
      if (obs > 0.0)
      {
        nll = nll - T(effective_n * obs) * log_t(max_t(predicted[i], floor));
      }
    }
    return nll;
  }

  class RedSnapperQuadraObjective
  {
  public:
    explicit RedSnapperQuadraObjective(std::vector<Observation> observations)
        : observations_(std::move(observations)) {}

    template <class T>
    T operator()(const std::vector<T> &par) const
    {
      if (par.size() < 5 + observations_.size())
      {
        throw std::runtime_error(
            "RedSnapperQuadraObjective expected parameters: log_r0, log_fbar, log_q");
      }

      const T log_r0 = par[0];
      const T log_fbar = par[1];
      const T log_q = par[2];
      const T logit_sel_a50 = par[3];
      const T log_sel_slope = par[4];

      const T r0 = exp_t(log_r0);
      const T m = T(0.18);
      const T fbar = exp_t(log_fbar);
      const T q = exp_t(log_q);
      const T sel_a50 = T(1.0) + T(9.0) * invlogit_t(logit_sel_a50);
      const T sel_slope = exp_t(log_sel_slope);

      const T sigma_log_index = T(0.20);
      const T sigma_log_catch = T(0.15);

      const T sigma_rec_dev = T(0.35);
      const double age_comp_effective_n = 2.0;
      const double min_positive = 1.0e-12;

      const auto weight = default_weight_at_age();
      const auto maturity = default_maturity_at_age();

      std::array<T, kAges> selectivity{};
      for (int a = 0; a < kAges; ++a)
      {
        selectivity[static_cast<std::size_t>(a)] =
            logistic_selectivity_t(T(a + 1), sel_a50, sel_slope);
      }

      std::array<T, kAges> n{};
      n[0] = r0;
      for (int a = 1; a < kAges; ++a)
      {
        n[static_cast<std::size_t>(a)] =
            n[static_cast<std::size_t>(a - 1)] * exp_t(-m);
      }
      n[static_cast<std::size_t>(kAges - 1)] =
          n[static_cast<std::size_t>(kAges - 1)] /
          (T(1.0) - exp_t(-m));

      T nll = T(0.0);
      T fixed_prior_nll = T(0.0);
      T rec_prior_nll = T(0.0);
      T index_nll = T(0.0);
      T catch_nll = T(0.0);
      T age_comp_nll_total = T(0.0);

      auto normal_prior = [](const T &x, double mean, double sd)
      {
        const T z = (x - T(mean)) / T(sd);
        return T(0.5) * z * z;
      };

      fixed_prior_nll = fixed_prior_nll + normal_prior(log_r0, std::log(1200.0), 1.0);
      fixed_prior_nll = fixed_prior_nll + normal_prior(log_fbar, std::log(0.025), 0.75);
      fixed_prior_nll = fixed_prior_nll + normal_prior(log_q, std::log(0.00005), 1.0);
      fixed_prior_nll = fixed_prior_nll + normal_prior(sel_a50, 4.0, 0.75);
      fixed_prior_nll = fixed_prior_nll + normal_prior(log_sel_slope, std::log(1.2), 0.35);

      nll = nll + fixed_prior_nll;

      for (std::size_t t = 0; t < observations_.size(); ++t)
      {

        const auto &obs = observations_[t];

        const T rec_dev = par[5 + t];

        {
          T term = T(0.5) * square_t(rec_dev / sigma_rec_dev);
          rec_prior_nll = rec_prior_nll + term;
          nll = nll + term;
        }
        T biomass = T(0.0);
        for (int a = 0; a < kAges; ++a)
        {
          biomass = biomass +
                    n[static_cast<std::size_t>(a)] *
                        T(weight[static_cast<std::size_t>(a)]);
        }

        T catch_hat = T(0.0);
        for (int a = 0; a < kAges; ++a)
        {
          const auto i = static_cast<std::size_t>(a);
          const T f_a = fbar * selectivity[i];
          const T z_a = m + f_a;
          const T harvest_rate =
              (f_a / z_a) * (T(1.0) - exp_t(-z_a));
          catch_hat = catch_hat + n[i] * T(weight[i]) * harvest_rate;
        }

        const T index_hat = q * biomass;

        if (obs.index > 0.0)
        {
          const T z = (log_t(T(obs.index)) -
                       log_t(max_t(index_hat, min_positive))) /
                      sigma_log_index;
          {
            T term = T(0.5) * square_t(z);
            index_nll = index_nll + term;
            nll = nll + term;
          }
        }

        if (obs.catch_mt > 0.0)
        {
          const T z = (log_t(T(obs.catch_mt)) -
                       log_t(max_t(catch_hat, min_positive))) /
                      sigma_log_catch;
          {
            T term = T(0.5) * square_t(z);
            catch_nll = catch_nll + term;
            nll = nll + term;
          }
        }

        std::array<T, kAges> pred_age_comp{};
        T selected_numbers_sum = T(0.0);
        for (int a = 0; a < kAges; ++a)
        {
          const auto i = static_cast<std::size_t>(a);
          pred_age_comp[i] = n[i] * selectivity[i];
          selected_numbers_sum = selected_numbers_sum + pred_age_comp[i];
        }
        for (int a = 0; a < kAges; ++a)
        {
          const auto i = static_cast<std::size_t>(a);
          pred_age_comp[i] =
              pred_age_comp[i] / max_t(selected_numbers_sum, min_positive);
        }

        {
          T term = age_comp_nll(obs.age_comp, pred_age_comp,
                                age_comp_effective_n, min_positive);
          age_comp_nll_total = age_comp_nll_total + term;
          nll = nll + term;
        }

        std::array<T, kAges> next{};
        next[0] = r0 * exp_t(rec_dev);

        for (int a = 1; a < kAges; ++a)
        {
          const auto prev = static_cast<std::size_t>(a - 1);
          const T f_prev = fbar * selectivity[prev];
          const T z_prev = m + f_prev;
          next[static_cast<std::size_t>(a)] = n[prev] * exp_t(-z_prev);
        }

        const auto last = static_cast<std::size_t>(kAges - 1);
        const T f_last = fbar * selectivity[last];
        const T z_last = m + f_last;
        next[last] = next[last] + n[last] * exp_t(-z_last);

        n = next;
      }

      return nll;
    }

  private:
    std::vector<Observation> observations_;
  };

  void write_fit_summary(const std::string &path,
                         const quadra::OptResult &fit)
  {
    std::ofstream out(path);
    if (!out)
    {
      throw std::runtime_error("Could not open fit summary CSV: " + path);
    }

    out << "field,value\n";
    out << std::setprecision(12);
    out << "objective," << fit.value << "\n";
    out << "joint_objective," << fit.joint_objective << "\n";
    out << "laplace_logdet," << fit.laplace_logdet << "\n";
    out << "laplace_constant," << fit.laplace_constant << "\n";
    out << "grad_norm," << fit.grad_norm << "\n";
    out << "iterations," << fit.iterations << "\n";
    out << "converged," << (fit.converged ? "yes" : "no") << "\n";
    out << "message," << fit.message << "\n";
    out << "laplace,yes\n";
    out << "random_effects," << fit.u_hat.size() << "\n";

    if (fit.par.size() >= 3)
    {
      out << "log_r0," << fit.par[0] << "\n";
      out << "r0," << std::exp(fit.par[0]) << "\n";
      out << "log_fbar," << fit.par[1] << "\n";
      out << "fbar," << std::exp(fit.par[1]) << "\n";
      out << "log_q," << fit.par[2] << "\n";
      out << "q," << std::exp(fit.par[2]) << "\n";
      if (fit.par.size() >= 5)
      {
        const double sel_a50 = 1.0 + 9.0 / (1.0 + std::exp(-fit.par[3]));
        const double sel_slope = std::exp(fit.par[4]);
        out << "logit_sel_a50," << fit.par[3] << "\n";
        out << "sel_a50," << sel_a50 << "\n";
        out << "log_sel_slope," << fit.par[4] << "\n";
        out << "sel_slope," << sel_slope << "\n";
      }
    }
  }

} // namespace sefsc_red_snapper

void write_fitted_trajectory(
    const std::string &path,
    const std::vector<sefsc_red_snapper::Observation> &observations,
    const quadra::OptResult &fit)
{
  if (fit.par.size() < 3)
  {
    throw std::runtime_error("Cannot write fitted trajectory: expected at least 3 fixed parameters");
  }

  sefsc_red_snapper::AgeStructuredParams params;
  params.log_r0 = fit.par[0];
  params.log_fbar = fit.par[1];
  params.log_q = fit.par[2];
  if (fit.par.size() >= 5)
  {
    params.sel_a50 = 1.0 + 9.0 / (1.0 + std::exp(-fit.par[3]));
    params.sel_slope = std::exp(fit.par[4]);
  }

  const auto rows =
      sefsc_red_snapper::run_deterministic_age_structured_model(observations,
                                                                params);

  std::ofstream out(path);
  if (!out)
  {
    throw std::runtime_error("Could not open fitted trajectory CSV: " + path);
  }

  out << "year,recruitment,total_biomass,ssb_proxy,depletion,Fbar,"
      << "catch_obs,catch_hat,catch_log_residual,index_obs,index_hat,"
      << "index_log_residual\n";

  out << std::fixed << std::setprecision(6);

  for (const auto &row : rows)
  {
    const double catch_log_residual =
        std::log(std::max(row.catch_obs, 1.0e-12)) -
        std::log(std::max(row.catch_hat, 1.0e-12));
    const double index_log_residual =
        std::log(std::max(row.index_obs, 1.0e-12)) -
        std::log(std::max(row.index_hat, 1.0e-12));

    out << row.year << "," << row.recruitment << "," << row.total_biomass
        << "," << row.ssb_proxy << "," << row.depletion << ","
        << row.fbar << "," << row.catch_obs << "," << row.catch_hat
        << "," << catch_log_residual << "," << row.index_obs << ","
        << row.index_hat << "," << index_log_residual << "\n";
  }
}

struct ResidualDiagnostics
{
  int n = 0;
  double catch_rmse_log = 0.0;
  double index_rmse_log = 0.0;
  double catch_mean_log_residual = 0.0;
  double index_mean_log_residual = 0.0;
  double max_abs_catch_log_residual = 0.0;
  double max_abs_index_log_residual = 0.0;
};

void write_residual_diagnostics(
    const std::string &path,
    const std::vector<sefsc_red_snapper::Observation> &observations,
    const quadra::OptResult &fit)
{
  sefsc_red_snapper::AgeStructuredParams params;
  params.log_r0 = fit.par[0];
  params.log_fbar = fit.par[1];
  params.log_q = fit.par[2];
  if (fit.par.size() >= 5)
  {
    params.sel_a50 = 1.0 + 9.0 / (1.0 + std::exp(-fit.par[3]));
    params.sel_slope = std::exp(fit.par[4]);
  }

  const auto rows =
      sefsc_red_snapper::run_deterministic_age_structured_model(observations,
                                                                params);

  ResidualDiagnostics d;
  d.n = static_cast<int>(rows.size());

  double catch_sum = 0.0, catch_ss = 0.0;
  double index_sum = 0.0, index_ss = 0.0;

  for (const auto &row : rows)
  {
    const double cr = std::log(std::max(row.catch_obs, 1.0e-12)) -
                      std::log(std::max(row.catch_hat, 1.0e-12));
    const double ir = std::log(std::max(row.index_obs, 1.0e-12)) -
                      std::log(std::max(row.index_hat, 1.0e-12));

    catch_sum += cr;
    catch_ss += cr * cr;
    index_sum += ir;
    index_ss += ir * ir;

    d.max_abs_catch_log_residual =
        std::max(d.max_abs_catch_log_residual, std::abs(cr));
    d.max_abs_index_log_residual =
        std::max(d.max_abs_index_log_residual, std::abs(ir));
  }

  if (d.n > 0)
  {
    d.catch_mean_log_residual = catch_sum / d.n;
    d.index_mean_log_residual = index_sum / d.n;
    d.catch_rmse_log = std::sqrt(catch_ss / d.n);
    d.index_rmse_log = std::sqrt(index_ss / d.n);
  }

  std::ofstream out(path);
  out << "metric,value,note\n";
  out << std::setprecision(12);
  out << "n," << d.n << ",number of fitted years\n";
  out << "catch_rmse_log," << d.catch_rmse_log << ",root mean squared log catch residual\n";
  out << "index_rmse_log," << d.index_rmse_log << ",root mean squared log index residual\n";
  out << "catch_mean_log_residual," << d.catch_mean_log_residual << ",mean log observed minus predicted catch\n";
  out << "index_mean_log_residual," << d.index_mean_log_residual << ",mean log observed minus predicted index\n";
  out << "max_abs_catch_log_residual," << d.max_abs_catch_log_residual << ",maximum absolute log catch residual\n";
  out << "max_abs_index_log_residual," << d.max_abs_index_log_residual << ",maximum absolute log index residual\n";
}

void write_selectivity_at_age(const std::string &path,
                              const quadra::OptResult &fit)
{
  if (fit.par.size() < 5)
  {
    return;
  }

  const double a50 = 1.0 + 9.0 / (1.0 + std::exp(-fit.par[3]));
  const double slope = std::exp(fit.par[4]);

  std::ofstream out(path);
  out << "age,selectivity\n";

  for (int age = 1; age <= sefsc_red_snapper::kAges; ++age)
  {
    const double sel = 1.0 / (1.0 + std::exp(-slope * (age - a50)));
    out << age << "," << sel << "\n";
  }
}

void write_recruitment_deviations(const std::string &path,
                                  const quadra::OptResult &fit)
{
  std::ofstream out(path);
  out << "year,log_rec_dev,rec_multiplier\n";
  out << std::setprecision(12);

  for (std::size_t i = 0; i < fit.u_hat.size(); ++i)
  {
    const double u = fit.u_hat[i];
    out << (i + 1) << "," << u << "," << std::exp(u) << "\n";
  }
}

void write_objective_components(
    const std::string &path,
    const std::vector<sefsc_red_snapper::Observation> &observations,
    const quadra::OptResult &fit)
{
  if (fit.par.size() < 5 || fit.u_hat.size() < observations.size())
  {
    throw std::runtime_error("Cannot write objective components: missing fit values");
  }

  const double log_r0 = fit.par[0];
  const double log_fbar = fit.par[1];
  const double log_q = fit.par[2];
  const double logit_sel_a50 = fit.par[3];
  const double log_sel_slope = fit.par[4];

  const double r0 = std::exp(log_r0);
  const double m = 0.18;
  const double fbar = std::exp(log_fbar);
  const double q = std::exp(log_q);
  const double sel_a50 = 1.0 + 9.0 / (1.0 + std::exp(-logit_sel_a50));
  const double sel_slope = std::exp(log_sel_slope);

  const double sigma_log_index = 0.20;
  const double sigma_log_catch = 0.15;
  const double sigma_rec_dev = 0.35;
  const double age_comp_effective_n = 2.0;
  const double min_positive = 1.0e-12;

  const auto weight = sefsc_red_snapper::default_weight_at_age();

  std::array<double, sefsc_red_snapper::kAges> selectivity{};
  for (int a = 0; a < sefsc_red_snapper::kAges; ++a)
  {
    selectivity[static_cast<std::size_t>(a)] =
        sefsc_red_snapper::logistic_selectivity(
            static_cast<double>(a + 1), sel_a50, sel_slope);
  }

  std::array<double, sefsc_red_snapper::kAges> n{};
  n[0] = r0;
  for (int a = 1; a < sefsc_red_snapper::kAges; ++a)
  {
    n[static_cast<std::size_t>(a)] =
        n[static_cast<std::size_t>(a - 1)] * std::exp(-m);
  }
  n[static_cast<std::size_t>(sefsc_red_snapper::kAges - 1)] =
      n[static_cast<std::size_t>(sefsc_red_snapper::kAges - 1)] /
      (1.0 - std::exp(-m));

  auto normal_prior = [](double x, double mean, double sd)
  {
    const double z = (x - mean) / sd;
    return 0.5 * z * z;
  };

  double fixed_prior_nll = 0.0;
  double rec_prior_nll = 0.0;
  double index_nll = 0.0;
  double catch_nll = 0.0;
  double age_comp_nll = 0.0;

  fixed_prior_nll += normal_prior(log_r0, std::log(1200.0), 1.0);
  fixed_prior_nll += normal_prior(log_fbar, std::log(0.025), 0.75);
  fixed_prior_nll += normal_prior(log_q, std::log(0.00005), 1.0);
  fixed_prior_nll += normal_prior(sel_a50, 4.0, 0.75);
  fixed_prior_nll += normal_prior(log_sel_slope, std::log(1.2), 0.35);

  for (std::size_t t = 0; t < observations.size(); ++t)
  {
    const auto &obs = observations[t];
    const double rec_dev = fit.u_hat[t];

    rec_prior_nll += 0.5 * std::pow(rec_dev / sigma_rec_dev, 2.0);

    double biomass = 0.0;
    for (int a = 0; a < sefsc_red_snapper::kAges; ++a)
    {
      biomass += n[static_cast<std::size_t>(a)] *
                 weight[static_cast<std::size_t>(a)];
    }

    double catch_hat = 0.0;
    for (int a = 0; a < sefsc_red_snapper::kAges; ++a)
    {
      const auto i = static_cast<std::size_t>(a);
      const double f_a = fbar * selectivity[i];
      const double z_a = m + f_a;
      const double harvest_rate = (f_a / z_a) * (1.0 - std::exp(-z_a));
      catch_hat += n[i] * weight[i] * harvest_rate;
    }

    const double index_hat = q * biomass;

    if (obs.index > 0.0)
    {
      const double z =
          (std::log(obs.index) - std::log(std::max(index_hat, min_positive))) /
          sigma_log_index;
      index_nll += 0.5 * z * z;
    }

    if (obs.catch_mt > 0.0)
    {
      const double z =
          (std::log(obs.catch_mt) -
           std::log(std::max(catch_hat, min_positive))) /
          sigma_log_catch;
      catch_nll += 0.5 * z * z;
    }

    std::array<double, sefsc_red_snapper::kAges> pred_age_comp{};
    double selected_numbers_sum = 0.0;
    for (int a = 0; a < sefsc_red_snapper::kAges; ++a)
    {
      const auto i = static_cast<std::size_t>(a);
      pred_age_comp[i] = n[i] * selectivity[i];
      selected_numbers_sum += pred_age_comp[i];
    }

    for (int a = 0; a < sefsc_red_snapper::kAges; ++a)
    {
      const auto i = static_cast<std::size_t>(a);
      pred_age_comp[i] =
          pred_age_comp[i] / std::max(selected_numbers_sum, min_positive);

      const double obs_a = std::max(obs.age_comp[i], 0.0);
      if (obs_a > 0.0)
      {
        age_comp_nll -= age_comp_effective_n * obs_a *
                        std::log(std::max(pred_age_comp[i], min_positive));
      }
    }

    std::array<double, sefsc_red_snapper::kAges> next{};
    next[0] = r0 * std::exp(rec_dev);
    for (int a = 1; a < sefsc_red_snapper::kAges; ++a)
    {
      const auto prev = static_cast<std::size_t>(a - 1);
      const auto cur = static_cast<std::size_t>(a);
      const double f_prev = fbar * selectivity[prev];
      const double z_prev = m + f_prev;
      next[cur] = n[prev] * std::exp(-z_prev);
    }

    const int plus_group = sefsc_red_snapper::kAges - 1;
    const auto pg = static_cast<std::size_t>(plus_group);
    const double f_pg = fbar * selectivity[pg];
    const double z_pg = m + f_pg;
    next[pg] += n[pg] * std::exp(-z_pg);

    n = next;
  }

  std::ofstream out(path);
  if (!out)
  {
    throw std::runtime_error("Could not open component CSV: " + path);
  }

  out << "component,value\n";
  out << std::setprecision(12);
  out << "fixed_prior_nll," << fixed_prior_nll << "\n";
  out << "rec_prior_nll," << rec_prior_nll << "\n";
  out << "index_nll," << index_nll << "\n";
  out << "catch_nll," << catch_nll << "\n";
  out << "age_comp_nll," << age_comp_nll << "\n";
  out << "joint_total,"
      << fixed_prior_nll + rec_prior_nll + index_nll + catch_nll + age_comp_nll
      << "\n";
}

int main()
{
  const std::string input_path =
      "examples/NMFS/sefsc_red_snapper/data/synthetic_red_snapper_observations.csv";
  const std::string summary_path =
      "examples/NMFS/sefsc_red_snapper/outputs/quadra_fit_summary.csv";
  const std::string trajectory_path =
      "examples/NMFS/sefsc_red_snapper/outputs/quadra_fitted_trajectory.csv";
  const std::string residual_diagnostics_path =
      "examples/NMFS/sefsc_red_snapper/outputs/quadra_fit_residual_diagnostics.csv";
  const std::string selectivity_path =
      "examples/NMFS/sefsc_red_snapper/outputs/selectivity_at_age.csv";
  const std::string recruitment_deviations_path =
      "examples/NMFS/sefsc_red_snapper/outputs/recruitment_deviations.csv";
  const std::string objective_components_path =
      "examples/NMFS/sefsc_red_snapper/outputs/quadra_fit_objective_components.csv";
  const auto observations = sefsc_red_snapper::read_observations(input_path);

  sefsc_red_snapper::RedSnapperQuadraObjective objective(observations);

  quadra::ParameterVector params;
  params.add({"log_r0", std::log(1200.0), quadra::ParameterTransform::Identity, false});
  params.add({"log_fbar", std::log(0.025), quadra::ParameterTransform::Identity, false});
  params.add({"log_q", std::log(0.00005), quadra::ParameterTransform::Identity, false});
  params.add({"logit_sel_a50", 0.0, quadra::ParameterTransform::Identity, false});
  params.add({"log_sel_slope", std::log(1.2), quadra::ParameterTransform::Identity, false});

  for (std::size_t t = 0; t < observations.size(); ++t)
  {
    params.add({"log_rec_dev_" + std::to_string(t + 1),
                0.0,
                quadra::ParameterTransform::Identity,
                true});
  }

  quadra::LaplaceOptions opts;

  auto fit = quadra::optimize_lbfgs(objective, params, opts);

  sefsc_red_snapper::write_fit_summary(summary_path, fit);
  write_fitted_trajectory(trajectory_path, observations, fit);
  write_residual_diagnostics(residual_diagnostics_path, observations, fit);
  write_selectivity_at_age(selectivity_path, fit);
  write_recruitment_deviations(recruitment_deviations_path, fit);
  write_objective_components(objective_components_path, observations, fit);
  std::cout
      << "SEFSC red-snapper-style Quadra Laplace recruitment-deviation fit\n";
  std::cout << "objective:  " << fit.value << "\n";
  std::cout << "grad_norm:  " << fit.grad_norm << "\n";
  std::cout << "converged:  " << (fit.converged ? "yes" : "no") << "\n";
  std::cout << "message:    " << fit.message << "\n";
  std::cout << "wrote:      " << summary_path << "\n";
  std::cout << "wrote:      " << trajectory_path << "\n";
  std::cout << "wrote:      " << residual_diagnostics_path << "\n";
  std::cout << "wrote:      " << selectivity_path << "\n";
  std::cout << "wrote:      " << recruitment_deviations_path << "\n";
  std::cout << "wrote:      " << objective_components_path << "\n";
  return 0;
}
