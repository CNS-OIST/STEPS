#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <map>
#include <utility>
#include <vector>

#include "math/tools.hpp"
#include "rng/create.hpp"
#include "util/error.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace steps::rng;
using namespace steps::math;


///--- sign function -----------------------------------------------------------
inline int sgn(double x) {
    return (x > 0) ? 1 : -1;
}

///--- CDF Normal( mu, sigma^2 ) -------------------------------------------------------
inline double cdf_normal_distribution(double x, double mean = 0., double stddev = 1.) {
    return 0.5 * (1. + std::erf((x - mean) / (stddev * sqrt(2.))));
}

///--- normalized incomplete gamma(a,x) ---------------------------------------
inline double norm_inc_gamma(double x, double a, double rtol = 1e-6, size_t iter_max = 100) {
    double term = 1;
    double sum = term;
    for (size_t k = 1; k < iter_max; ++k) {
        term = std::exp(std::log(term) + std::log(x) - std::log(a + k));
        sum += term;
        if (std::abs(term) < rtol * std::abs(sum))
            return std::exp(std::log(sum) + a * std::log(x) - x - lgamma(a + 1));
    }
    throw(std::runtime_error("Max iterations in incomplete gamma computation!"));
}

///--- CDF chi_sq( nu ) -------------------------------------------------------
inline double cdf_chi_squared_distribution(double x,
                                           double nu,
                                           double rtol = 1e-6,
                                           size_t iter_max = 1000) {
    return norm_inc_gamma(x / 2., nu / 2., rtol, iter_max);
}

/// Chi-squared goodness of fit test for binomial distribution
void binomial_check(const std::string& str,
                    uint t,
                    double p,
                    uint t_test,
                    double p_test,
                    const uint n_sample,
                    double level_confidence) {
    std::vector<uint> vec_rng(n_sample);

    auto rng = create(str, n_sample);
    rng->initialize(12345u);

    /// Generating samples from the supposed binomial
    for (uint i = 0; i < n_sample; ++i)
        vec_rng[i] = rng->getBinom(t, p);

    /// Calculating observed and expected counts from binomial
    std::vector<uint> vec_observed(t + 1, 0);
    std::vector<double> vec_expected(t + 1, 0);
    for (uint n = 0; n < n_sample; ++n)
        ++vec_observed[vec_rng[n]];
    for (uint k = 0; k <= t; ++k)
        vec_expected[k] = n_sample * binom_coeff(t_test, k) * pow(p_test, k) *
                          pow(1 - p_test, t_test - k);

    /// Compute chi-squared statistic
    double chi = 0.0;
    for (uint k = 0; k <= t; ++k)
        chi += (vec_expected[k] - vec_observed[k]) * (vec_expected[k] - vec_observed[k]) /
               vec_expected[k];

    double p_value = 1 - cdf_chi_squared_distribution(chi, n_sample - 1);

    if ((t == t_test) && (fabs(p - p_test) < 1e-14 * p)) {
        REQUIRE(p_value >= 1 - level_confidence);
        std::cerr << "FAILED Goodness of Fit test (chi-squared) with level of confidence "
                  << level_confidence << ": p_value is " << p_value
                  << " which is not big enough to support the hypothesis of binomial distribution"
                  << std::endl;
    } else {
        REQUIRE(p_value <= level_confidence);
        std::cerr
            << "This test was supposed to FAIL goodness of Fit test (chi-squared) with level of "
               "confidence "
            << level_confidence << ": p_value is " << p_value
            << " which is not small enough to reject the hypothesis that the sample follows the "
               "theoretical distribution"
            << std::endl;
    }
}

/// Kendall (type A) test to verify the independency of streams.
//  If p-value < (1-level_confidence) then the streams are correlated, and our test fails.
//  The p-value expression is that of a two-tailed statistical test.
double kendall_rank_correlation_check(const std::string& str,
                                      const uint n_sample,
                                      unsigned long seed1,
                                      unsigned long seed2) {
    std::vector<uint> samples1(n_sample);
    std::vector<uint> samples2(n_sample);

    auto rng1 = create(str, n_sample);
    rng1->initialize(seed1);
    auto rng2 = create(str, n_sample);
    rng2->initialize(seed2);

    for (uint i = 0; i < n_sample; ++i) {
        samples1[i] = rng1->get();
        samples2[i] = rng2->get();
    }

    long nc_nd = 0;
    for (uint i = 1; i < n_sample; ++i)
        for (uint j = 0; j < i; ++j)
            nc_nd += sgn(static_cast<long long>(samples1[i]) -
                         static_cast<long long>(samples1[j])) *
                     sgn(static_cast<long long>(samples2[i]) - static_cast<long long>(samples2[j]));

    const double tau = nc_nd / (n_sample * (n_sample - 1.) / 2.);
    const double sigma = std::sqrt(2. * (2. * n_sample + 5.) / (9. * n_sample * (n_sample - 1.)));
    const double p_value = 2. * (1 - cdf_normal_distribution(std::abs(tau / sigma)));

    return p_value;
}


void assert_pvalue(const double p_value, const double level_confidence) {
    REQUIRE(p_value >= 0.);
    REQUIRE(p_value <= 1.0);
    REQUIRE(p_value >= 1. - level_confidence);
    std::cerr << "Kendall correlation test has p_value = " << p_value
              << " below the threshold = " << (1. - level_confidence)
              << ". This indicates correlation in the streams" << std::endl;
}


TEST_CASE("rng_class_create_mt") {
    const uint n_sample = 100;

    REQUIRE_NOTHROW(create("mt19937", n_sample));
    REQUIRE_THROWS_AS(create("mt19937", 0), steps::ArgErr);
    REQUIRE_THROWS_AS(create("abc", n_sample), steps::ArgErr);
}

TEST_CASE("rng_class_create_r123") {
    const uint n_sample = 100;

    REQUIRE_NOTHROW(create("r123", n_sample));
    REQUIRE_THROWS_AS(create("r123", 0), steps::ArgErr);
    REQUIRE_THROWS_AS(create("abc", n_sample), steps::ArgErr);
}

TEST_CASE("rng_binomial_true_mt") {
    binomial_check("mt19937", 10, 0.3, 10, 0.3, 10000, 0.999);
    binomial_check("mt19937", 8, 0.4, 8, 0.4, 10000, 0.999);
}

TEST_CASE("rng_binomial_true_r123") {
    binomial_check("r123", 10, 0.3, 10, 0.3, 10000, 0.999);
    binomial_check("r123", 8, 0.4, 8, 0.4, 10000, 0.999);
}

TEST_CASE("rng_binomial_false_mt") {
    binomial_check("mt19937", 10, 0.3, 10, 0.45, 100, 0.9);  /// supposed to fail
    binomial_check("mt19937", 10, 0.3, 15, 0.4, 100, 0.9);   /// supposed to fail
}

TEST_CASE("rng_binomial_false_r123") {
    binomial_check("r123", 10, 0.3, 10, 0.45, 100, 0.9);  /// supposed to fail
    binomial_check("r123", 10, 0.3, 15, 0.4, 100, 0.9);   /// supposed to fail
}


TEST_CASE("rng_kendall_mt") {
    const uint n_sample = 10000;
    const double level_confidence = 0.95;
    const unsigned long seed1 = 1;
    const unsigned long seed2 = 2;
    const double p_value = kendall_rank_correlation_check("mt19937", n_sample, seed1, seed2);
    assert_pvalue(p_value, level_confidence);
}


TEST_CASE("rng_kendall_r123") {
    const uint n_sample = 10000;
    const double level_confidence = 0.95;
    const unsigned long seed1 = 1;
    const unsigned long seed2 = 2;
    const double p_value = kendall_rank_correlation_check("r123", n_sample, seed1, seed2);
    assert_pvalue(p_value, level_confidence);
}
