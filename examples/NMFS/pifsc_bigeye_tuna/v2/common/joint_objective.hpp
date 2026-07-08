#pragma once

#include "agecomp.hpp"
#include "catch.hpp"
#include "index.hpp"
#include "life_history.hpp"
#include "likelihood.hpp"
#include "mortality.hpp"
#include "population.hpp"
#include "recruitment.hpp"
#include "selectivity.hpp"

namespace bigeye_v2 {

struct BigeyeJointObjective {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const BigeyeModelParameters<T> &p,
                  BigeyeDerived<T> &d) const {
    d.total_nll = T(0.0);
    d.catch_nll = T(0.0);
    d.index_nll = T(0.0);
    d.agecomp_nll = T(0.0);

    BigeyeLifeHistory{}(data, p, d);
    FixedRecruitment{}(data, p, d);
    LogisticSelectivity{}(data, p, d);
    FishingMortality{}(data, p, d);
    UnfishedPopulation{}(data, p, d);
    BaranovCatch{}(data, p, d);
    BiomassIndexPrediction{}(data, p, d);
    CatchAgeCompositionPrediction{}(data, p, d);

    LognormalCatchLikelihood{}(data, p, d);
    LognormalIndexLikelihood{}(data, p, d);
    MultinomialAgeCompLikelihood{}(data, p, d);
  }
};

} // namespace bigeye_v2
