#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include "steps/rng/small_binomial.hpp"

#include "gtest/gtest.h"

using namespace steps::rng;


// - - - One-Sample Kolmogorov-Smirnov Goodness of Fit - - -
// Measures discrepancy of empirical cdf w.r.t. theoretical one
// Significance level for rejecting H0: alpha = 0.01
TEST(SmallBinomial, KolmogorovSmirnov) {
    // sample size for test
    constexpr int n_samples = 1e3;
    // critical value for discrepancy at level 0.01
    constexpr double D_critical = 0.051545;

    // define probability p and number of trials t
    constexpr double p = 0.2;
    constexpr int t = 3;
    // exact cumulative distribution function
    constexpr double cdf_exa[t+1] = {0.512, 0.896, 0.992, 1.};

    // compute empirical cdf
    std::random_device rd;
    std::mt19937 gen(rd());
    steps::rng::small_binomial_distribution<> X(t, p);
    double cdf_emp[t+1] = {0, 0, 0, 0};
    for (int i=0; i<n_samples; ++i) {
        int x = X(gen);
        for (int j=0; j<t+1; ++j)
            if (x<=j)
                cdf_emp[j] +=1;
    }
    for (auto& val : cdf_emp)
        val /= n_samples;

    // compute discrepancy of the cdf and compare with critical
    double D_sample =  std::abs(cdf_emp[0] - cdf_exa[0]);
    for (int i=0; i<t+1; ++i) {
        double Di = std::abs(cdf_emp[i] - cdf_exa[i]);
        if (Di > D_sample)
            D_sample = Di;
    }
    ASSERT_LT(D_sample, D_critical) << "Note: this test can fail 1% of the time for stochastic reasons";
}


// - - - One-Sample Chi-Squared Goodness of Fit - - -
// Measures distance from expected counts for each possible value
//Significance level for rejecting H0: alpha = 0.01
TEST(SmallBinomial, ChiSquared) {
    // sample size for test
    constexpr int n_samples = 1e3;

    // define probability p and number of trials t
    constexpr double p = 0.55;
    constexpr int t = 10;
    // exact probability denisity function
    constexpr double pdf_exa[t+1] = { 0.00034051, 0.00416174, 0.02288959,
                                      0.07460311, 0.15956775, 0.23403271, 0.23836665, 0.16647829, 0.07630255,
                                      0.02072415, 0.00253295};

    // compute observed counts

    // we need a static seed (anything that worked once works) otherwise
    // this test has 1% fail chance even if everything works as intended
    std::random_device rd;
    std::mt19937 gen(rd());
    steps::rng::small_binomial_distribution<> X(t, p);
    int counts_obs[t+1] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    for (int i=0; i<n_samples; ++i)
        ++counts_obs[X(gen)];

    // critical value for chi-squared distribution with t degrees of freedom, at level 0.01
    constexpr double chi2_critical = 23.209;

    // compute chi-squared test statistsics and compare with critical
    double chi2_sample = 0;
    for (int i=0; i<t+1; ++i) {
        double counts_exp = n_samples*pdf_exa[i];
        chi2_sample += (counts_obs[i]-counts_exp)*(counts_obs[i]-counts_exp)/counts_exp;
    }

    ASSERT_LT(chi2_sample, chi2_critical) << "Note: this test can fail 1% of the time for stochastic reasons";
}


int main(int argc, char **argv) {
    int r=0;
    ::testing::InitGoogleTest(&argc, argv);
    r=RUN_ALL_TESTS();
    return r;
}
