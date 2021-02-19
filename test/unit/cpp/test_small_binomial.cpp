#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include "steps/rng/small_binomial.hpp"

#include "gtest/gtest.h"

using namespace steps::rng;



// - - - One-Sample Chi-Squared Goodness of Fit - - -
// Measures distance from expected counts for each possible value
//Significance level for rejecting H0: alpha = 0.01
TEST(SmallBinomial, Density) {
    // sample size for test
    constexpr int n_samples = 1e3;

    // define probability p and number of trials t
    constexpr double p = 0.55;
    constexpr int t = 20;

    const int sol[t+1] = {0, 0, 0, 0, 1, 8, 12, 39, 75, 123, 157, 176, 147, 121, 83, 35, 16, 5, 2, 0, 0}; // done with seed = 3

    std::mt19937 gen(3);
    steps::rng::small_binomial_distribution<> X(t, p);
    int counts_obs[t+1] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    for (int i=0; i<n_samples; ++i)
        ++counts_obs[X(gen)];

    // compute err
    int err = 0;
    for (int i=0; i<t+1; ++i) {
        err += abs(sol[i]- counts_obs[i]);
    }
    ASSERT_EQ(err, 0);
}



TEST(smallBinomial, RNGlimits) {
  // Struct that fakes a "broken" rng generator that produces only its limits
  struct OnlyLimitsRNG {
    OnlyLimitsRNG() {}
    static constexpr uint min() { return 0; }
    static constexpr uint max() { return 0xffffffffu; }
    uint operator()() { parity = !parity; return parity ? max() : min(); }
  private:
    bool parity = 0;
  };

  constexpr long long int n_samples = 7;
  constexpr int t = 4;
  OnlyLimitsRNG olrng;

  {
    // with p = 1.0 every bernulli trial is a success
    constexpr double p = 1.0;
    steps::rng::small_binomial_distribution<> X(t, p);

    int counts_obs[t+1] = {0, 0, 0, 0, 0};
    for (int i=0; i<n_samples; ++i)
      ++counts_obs[X(olrng)];

    ASSERT_EQ(counts_obs[t], n_samples);
  }

  {
    // with p = 0.0 every bernulli trial is a failure
    constexpr double p = 0.0;
    steps::rng::small_binomial_distribution<> X(t, p);

    int counts_obs[t+1] = {0, 0, 0, 0, 0};
    for (int i=0; i<n_samples; ++i)
      ++counts_obs[X(olrng)];

    ASSERT_EQ(counts_obs[0], n_samples);
  }
}


int main(int argc, char **argv) {
    int r=0;
    ::testing::InitGoogleTest(&argc, argv);
    r=RUN_ALL_TESTS();
    return r;
}

