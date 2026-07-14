#pragma once

#include "bigeye_level21_parameter_layout.hpp"

#include "../objective/bigeye_quadra_objective.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>

namespace pifsc_bigeye_tuna {

using namespace level21_layout;

inline double level19_safe_logit(double p) {
  p = std::min(0.999, std::max(0.001, p));
  return std::log(p / (1.0 - p));
}

inline void
write_level21_parameter_sanity_diagnostics(const std::string &txt_path,
                                           const std::string &csv_path,
                                           const quadra::OptResult &fit) {
  std::ofstream txt(txt_path);
  std::ofstream csv(csv_path);
  if (!txt || !csv) {
    throw std::runtime_error(
        "Cannot write Level 21 parameter sanity diagnostics");
  }

  txt << std::setprecision(15);
  csv << std::setprecision(15);

  constexpr int kBaseFixed = 3;
  constexpr int kMParamOffset = kBaseFixed;
  constexpr int kMParams = 2;
  constexpr int kLonglineSelOffset = kMParamOffset + kMParams;
  constexpr int kLonglineSelDevs = kAges;
  constexpr int kInitialDevOffset = kLonglineSelOffset + kLonglineSelDevs;
  constexpr int kInitialDevs = kAges;
  constexpr int kPurseSeineSelOffset = kInitialDevOffset + kInitialDevs;
  constexpr int kPurseSeineSelDevs = kAges;

  const double sigma_init_dev = 0.75;
#ifdef BIGEYE_LL_SEL_SIGMA
  const double sigma_ll_sel_dev = static_cast<double>(BIGEYE_LL_SEL_SIGMA);
#else
  const double sigma_ll_sel_dev = 1.0;
#endif
  const double sigma_ps_sel_dev = 1.0;
  const double sigma_rec_dev = 0.35;

  const std::array<double, kAges> ll_template = {
      0.001, 0.005, 0.03, 0.12, 0.35, 0.70, 0.95, 0.90, 0.60, 0.12};
  const std::array<double, kAges> ps_template = {
      0.20, 0.90, 1.00, 0.80, 0.45, 0.20, 0.08, 0.03, 0.015, 0.005};

  auto sq = [](double x) { return x * x; };
  auto invlogit = [](double x) { return 1.0 / (1.0 + std::exp(-x)); };

  double m_age_penalty = 0.0;
  double init_penalty = 0.0;
  double ll_sel_penalty = 0.0;
  double ps_sel_penalty = 0.0;
  double rec_penalty = 0.0;

  double max_abs_init = 0.0;
  int max_abs_init_age = -1;
  double max_init_multiplier = -std::numeric_limits<double>::infinity();
  int max_init_multiplier_age = -1;
  double min_init_multiplier = std::numeric_limits<double>::infinity();
  int min_init_multiplier_age = -1;

  const double sigma_log_m_age_offset = 0.35;
  const double log_m_young_offset = fit.par[kMParamOffset + 0];
  const double log_m_old_offset = fit.par[kMParamOffset + 1];
  m_age_penalty += 0.5 * sq(log_m_young_offset / sigma_log_m_age_offset);
  m_age_penalty += 0.5 * sq(log_m_old_offset / sigma_log_m_age_offset);

  for (int a = 0; a < kAges; ++a) {
    const double v = fit.par[kInitialDevOffset + a];
    init_penalty += 0.5 * sq(v / sigma_init_dev);
    const double mult = std::exp(v);
    if (std::abs(v) > max_abs_init) {
      max_abs_init = std::abs(v);
      max_abs_init_age = a + 1;
    }
    if (mult > max_init_multiplier) {
      max_init_multiplier = mult;
      max_init_multiplier_age = a + 1;
    }
    if (mult < min_init_multiplier) {
      min_init_multiplier = mult;
      min_init_multiplier_age = a + 1;
    }
  }

  for (int a = 0; a < kAges; ++a) {
    const double raw = fit.par[kLonglineSelOffset + a];
    const double prior =
        level19_safe_logit(ll_template[static_cast<std::size_t>(a)]);
    ll_sel_penalty += 0.5 * sq((raw - prior) / sigma_ll_sel_dev);
  }

  for (int a = 0; a < kAges; ++a) {
    const double raw = fit.par[kPurseSeineSelOffset + a];
    const double prior =
        level19_safe_logit(ps_template[static_cast<std::size_t>(a)]);
    ps_sel_penalty += 0.5 * sq((raw - prior) / sigma_ps_sel_dev);
  }

  for (Eigen::Index i = 0; i < fit.u_hat.size(); ++i) {
    const double v = fit.u_hat[i];
    rec_penalty += 0.5 * sq(v / sigma_rec_dev);
  }

  const double total_flexible_penalty = m_age_penalty + init_penalty +
                                        ll_sel_penalty + ps_sel_penalty +
                                        rec_penalty;

  txt << "Level 21 Parameter Sanity Diagnostics\n";
  txt << "=====================================\n\n";
  txt << "Purpose\n";
  txt << "-------\n";
  txt << "Check whether Level 21 improved the age-composition fit by "
         "creating\n";
  txt << "pathological tradeoffs among initial numbers, longline age "
         "selectivity,\n";
  txt << "purse-seine age selectivity, and recruitment deviations.\n\n";

  csv << "section,name,value,note\n";

  txt << "Prior penalty by block\n";
  txt << "----------------------\n";
  txt << "age_based_m_prior_nll," << m_age_penalty << "\n";
  txt << "m_young," << 0.45 * std::exp(log_m_young_offset) << "\n";
  txt << "m_adult," << 0.45 << "\n";
  txt << "m_old," << 0.45 * std::exp(log_m_old_offset) << "\n";
  txt << "initial_numbers_prior_nll," << init_penalty << "\n";
  txt << "longline_selectivity_prior_nll," << ll_sel_penalty << "\n";
  txt << "purse_seine_selectivity_prior_nll," << ps_sel_penalty << "\n";
  txt << "recruitment_prior_nll," << rec_penalty << "\n";
  txt << "total_flexible_prior_nll," << total_flexible_penalty << "\n\n";

  csv << "prior_block,age_based_m_prior_nll," << m_age_penalty
      << ",sigma=0.35 on young/old log offsets\n";
  csv << "m_at_age,m_young," << 0.45 * std::exp(log_m_young_offset)
      << ",ages 1-3\n";
  csv << "m_at_age,m_adult," << 0.45 << ",ages 4-7\n";
  csv << "m_at_age,m_old," << 0.45 * std::exp(log_m_old_offset)
      << ",ages 8-10\n";
  csv << "prior_block,initial_numbers_prior_nll," << init_penalty
      << ",sigma=0.75\n";
  csv << "prior_block,longline_selectivity_prior_nll," << ll_sel_penalty
      << ",sigma=1.0 around template logits\n";
  csv << "prior_block,purse_seine_selectivity_prior_nll," << ps_sel_penalty
      << ",sigma=1.0 around template logits\n";
  csv << "prior_block,recruitment_prior_nll," << rec_penalty << ",sigma=0.35\n";
  csv << "prior_block,total_flexible_prior_nll," << total_flexible_penalty
      << ",sum of diagnostic blocks\n";

  txt << "Initial-number stress summary\n";
  txt << "-----------------------------\n";
  txt << "max_abs_init_log_dev," << max_abs_init << "\n";
  txt << "max_abs_init_log_dev_age," << max_abs_init_age << "\n";
  txt << "min_init_multiplier," << min_init_multiplier << "\n";
  txt << "min_init_multiplier_age," << min_init_multiplier_age << "\n";
  txt << "max_init_multiplier," << max_init_multiplier << "\n";
  txt << "max_init_multiplier_age," << max_init_multiplier_age << "\n\n";

  csv << "initial_stress,max_abs_init_log_dev," << max_abs_init << ",age "
      << max_abs_init_age << "\n";
  csv << "initial_stress,min_init_multiplier," << min_init_multiplier << ",age "
      << min_init_multiplier_age << "\n";
  csv << "initial_stress,max_init_multiplier," << max_init_multiplier << ",age "
      << max_init_multiplier_age << "\n";

  txt << "Age-specific parameters\n";
  txt << "-----------------------\n";
  txt << "age,init_log_dev,init_multiplier,ll_logit,ll_selectivity,ll_"
         "template,";
  txt << "ll_logit_minus_template,ps_logit,ps_selectivity,ps_template,";
  txt << "ps_logit_minus_template\n";

  csv << "age_parameter,age,init_log_dev,init_multiplier,ll_logit,ll_"
         "selectivity,";
  csv << "ll_template,ll_logit_minus_template,ps_logit,ps_selectivity,ps_"
         "template,";
  csv << "ps_logit_minus_template\n";

  for (int a = 0; a < kAges; ++a) {
    const auto i = static_cast<std::size_t>(a);
    const double init = fit.par[kInitialDevOffset + a];
    const double ll_raw = fit.par[kLonglineSelOffset + a];
    const double ll_prior = level19_safe_logit(ll_template[i]);
    const double ps_raw = fit.par[kPurseSeineSelOffset + a];
    const double ps_prior = level19_safe_logit(ps_template[i]);

    txt << (a + 1) << "," << init << "," << std::exp(init) << "," << ll_raw
        << "," << invlogit(ll_raw) << "," << ll_template[i] << ","
        << (ll_raw - ll_prior) << "," << ps_raw << "," << invlogit(ps_raw)
        << "," << ps_template[i] << "," << (ps_raw - ps_prior) << "\n";

    csv << "age_parameter," << (a + 1) << "," << init << "," << std::exp(init)
        << "," << ll_raw << "," << invlogit(ll_raw) << "," << ll_template[i]
        << "," << (ll_raw - ll_prior) << "," << ps_raw << ","
        << invlogit(ps_raw) << "," << ps_template[i] << ","
        << (ps_raw - ps_prior) << "\n";
  }

  txt << "\nInterpretation\n";
  txt << "--------------\n";
  txt << "Large initial-number multipliers, especially paired with large "
         "selectivity\n";
  txt << "departures, indicate that Level 21 may be fitting composition data "
         "by\n";
  txt << "trading initial age structure against selectivity. If the "
         "initial-number\n";
  txt << "prior penalty dominates the improvement in objective, the next step "
         "should\n";
  txt << "be either tightening the initial-number prior or reducing "
         "selectivity freedom.\n";
}

} // namespace pifsc_bigeye_tuna

using pifsc_bigeye_tuna::write_level21_parameter_sanity_diagnostics;
