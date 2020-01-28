#include <cmath>
#include <vector>
#include <map>
#include <utility>
#include <iostream>
#include <iterator>
#include <algorithm>

#include "steps/rng/create.hpp"
#include "steps/error.hpp"
#include "steps/math/tools.hpp"

#include "gtest/gtest.h"

using namespace steps::rng;
using namespace steps::math;


///--- sign function -----------------------------------------------------------
inline int sgn(double x) {
   return (x>0)?1:-1;
}

///--- CDF Normal( mu, sigma^2 ) -------------------------------------------------------
inline double cdf_normal_distribution( double x,
                                       double mean=0.,
                                       double stddev=1.) {
    return 0.5 * ( 1. + std::erf((x-mean)/(stddev*sqrt(2.))) );
}

///--- normalized incomplete gamma(a,x) ---------------------------------------
inline double norm_inc_gamma(double x,
                             double a,
                             double rtol=1e-6,
                             size_t iter_max=100) {
    double term = 1;
    double sum  = term;
    for (size_t k=1; k<iter_max; ++k) {
       term = std::exp ( std::log(term) +  std::log(x) - std::log(a+k) );
       sum += term;
       if (std::abs(term) < rtol * std::abs(sum) )
           return std::exp( std::log(sum) + a * std::log(x) - x - lgamma(a+1) );
    }
    throw(std::runtime_error("Max iterations in incomplete gamma computation!"));
}

///--- CDF chi_sq( nu ) -------------------------------------------------------
inline double cdf_chi_squared_distribution(double x,
                                            double nu,
                                            double rtol=1e-6,
                                            size_t iter_max=1000) {
    return norm_inc_gamma(x/2.,nu/2.,rtol,iter_max);
}

/// Chi-squared goodness of fit test for binomial distribution
void binomial_check(const std::string &str, uint t, double p, uint t_test, double p_test, const uint n_sample, double level_confidence) {
    const uint n_dof = n_sample - 1;
    std::vector<uint> vec_rng(n_sample);

    auto rng = create(str, n_sample);
    rng->initialize(12345u);

    /// Generating samples from the supposed binomial
    for (uint i = 0; i < n_sample; ++i)
        vec_rng[i] = rng->getBinom(t, p);

    /// Calculating observed and expected counts from binomial
    std::vector<uint>   vec_observed(t + 1, 0);
    std::vector<double> vec_expected(t + 1, 0);
    for (uint n = 0; n < n_sample; ++n)
        ++vec_observed[vec_rng[n]];
    for (uint k = 0; k <= t; ++k)
        vec_expected[k] = n_sample * binom_coeff(t_test, k) * pow(p_test, k) * pow (1 - p_test, t_test - k);

    /// Compute chi-squared statistic
    double chi = 0.0;
    for (uint k = 0; k <= t; ++k)
        chi += (vec_expected[k]-vec_observed[k]) * (vec_expected[k]-vec_observed[k]) / static_cast<double>(vec_expected[k]);

    std::pair<double, uint> qpair = std::make_pair(level_confidence, n_sample-1);

    double p_value = cdf_chi_squared_distribution(chi, n_sample - 1);
    if ((t == t_test) && (fabs(p - p_test) < 1e-14*p))
    {
        ASSERT_LE(p_value, 1 - level_confidence) << "FAILED Goodness of Fit test (chi-squared) with level of confidence "
                                                 << level_confidence << ": p_value is " << p_value
                                                 << " which is not small enough to support the hypothesis of binomial distribution" << std::endl;
    }
    else
    {
        ASSERT_GT(p_value, level_confidence) << "This test was supposed to FAIL goodness of Fit test (chi-squared) with level of confidence "
                                             << level_confidence << ": p_value is " << p_value
                                             << " which is not large enough to reject the hypothesis of binomial distribution" << std::endl;
    }
}

/// Kendall rank correlation test for the independency of streams
void kendall_rank_correlation_check(const std::string &str, const uint n_sample, double level_confidence, ulong seed1, ulong seed2) {
    std::vector<uint>  samples1(n_sample);
    std::vector<uint>  samples2(n_sample);

    auto rng1 = create(str, n_sample);
    rng1->initialize(seed1);
    auto rng2 = create(str, n_sample);
    rng2->initialize(seed2);

    for (uint i=0; i < n_sample; ++i) {
        samples1[i] = rng1->get();
        samples2[i] = rng2->get();
    }

    long nc_nd = 0;
    for (uint i = 1; i < n_sample; ++i)
        for (uint j = 0; j < i; ++j)
            nc_nd +=  sgn((int)samples1[i] - (int)samples1[j]) * sgn((int)samples2[i] - (int)samples2[j]);

    const double tau    = nc_nd / (n_sample*(n_sample-1.)/2.);
    const double sigma  = std::sqrt( 2.*(2.*n_sample+5.) / (9.*n_sample*(n_sample-1.)) );

    const double p_value = 2. * cdf_normal_distribution(tau/sigma) - 1.;

    ASSERT_LE(p_value, 1. - level_confidence) << "Failed Kendall correlation test with the level of confidence " << level_confidence
                                              << ". p_value " << p_value << " is not large enough to reject the hypothesis of correlated streams" << std::endl;
}


TEST(rng, class_create_mt) {
    const uint n_sample = 100;

    ASSERT_NO_THROW(create("mt19937", n_sample));
    ASSERT_THROW(create("abc", n_sample), steps::ArgErr);
}

TEST(rng, class_create_r123) {
    const uint n_sample = 100;

    ASSERT_NO_THROW(create("r123", n_sample));
    ASSERT_THROW(create("abc", n_sample), steps::ArgErr);
}

TEST(rng, binomial_true_mt) {
    binomial_check("mt19937", 10, 0.3, 10, 0.3, 10000, 0.999); 
    binomial_check("mt19937", 8, 0.4, 8, 0.4, 10000, 0.999); 
}

TEST(rng, binomial_true_r123) {
    binomial_check("r123", 10, 0.3, 10, 0.3, 10000, 0.999); 
    binomial_check("r123", 8, 0.4, 8, 0.4, 10000, 0.999); 
}

TEST(rng, binomial_false_mt) {
    binomial_check("mt19937", 10, 0.3, 10, 0.45, 100, 0.9); /// supposed to fail
    binomial_check("mt19937", 10, 0.3, 15, 0.4, 100, 0.9); /// supposed to fail
}

TEST(rng, binomial_false_r123) {
    binomial_check("r123", 10, 0.3, 10, 0.45, 100, 0.9); /// supposed to fail
    binomial_check("r123", 10, 0.3, 15, 0.4, 100, 0.9); /// supposed to fail
}


TEST(rng, kendall_mt) {
    kendall_rank_correlation_check("mt19937", 1000, 0.95, 1, 2);
}


TEST(rng, kendall_r123) {
    kendall_rank_correlation_check("r123", 1000, 0.95, 1, 2);
}

