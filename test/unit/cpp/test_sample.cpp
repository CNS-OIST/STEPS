#include <random>
#include <iterator>

#include "steps/math/sample.hpp"
#include "steps/util/distribute.hpp"
#include "gtest/gtest.h"

using namespace steps::math;
using namespace steps::util;

TEST(Sample, adjpareto_param) {
    using namespace std;
    double weights[] = {0.1,0.2,0.3,0.5,0.8,0.7,0.7,0.8,0.9};

    adjusted_pareto_sampler<double> S1(5, begin(weights), end(weights));
    ASSERT_EQ(5, S1.size());

    adjusted_pareto_sampler<double> S2;
    ASSERT_EQ(0, S2.size());

    S2.param(S1.param());
    ASSERT_EQ(S1.size(), S2.size());
}

TEST(Sample, adjpareto_allpop) {
    using namespace std;
    minstd_rand R;
    double weights[] = {1,1,1,1};
    
    adjusted_pareto_sampler<double> S(4, begin(weights), end(weights));

    const char *pop = "abcd";
    char sample[4];
    size_t n = S(pop, pop+4, sample, R);

    ASSERT_EQ(4,n);
    int na = 0,nb = 0,nc = 0,nd = 0;

    for (auto c: sample) {
        switch (c) {
        case 'a': ++na; break;
        case 'b': ++nb; break;
        case 'c': ++nc; break;
        case 'd': ++nd; break;
        default:
            FAIL() << "Found sample item not from population";
        }
    }

    ASSERT_EQ(1, na);
    ASSERT_EQ(1, nb);
    ASSERT_EQ(1, nc);
    ASSERT_EQ(1, nd);
}

TEST(Sample, adjpareto_fair) {
    using namespace std;
    minstd_rand R;

    vector<int> population = {0,1,2,3,4,5,6,7,8,9};
    vector<double> weights = {0.2,0.2,0.2,0.2,0.2,
                             0.6,0.6,0.6,0.6,0.6};

    vector<int> sample(4);
    vector<double> histogram(population.size());

    adjusted_pareto_sampler<double> S(sample.size(),begin(weights),end(weights));

    int N=10000;
    for (int i=0; i<N; ++i) {
        S(population.begin(),population.end(),sample.begin(),R);
        for (auto i: sample) ++histogram.at(i);
    }

    for (auto &x: histogram) x/=N;

    // Check that each observed mean is within (say) 4 standard deviations of
    // inclusion probability pi. sd = sqrt(pi(1-pi)/N)
    
    for (size_t i=0; i<histogram.size(); ++i) {
        double z = std::abs(histogram[i]-weights[i])/std::sqrt(weights[i]*(1-weights[i])/N);
        EXPECT_LT(z,4);
    }
}


TEST(Sample,distribute_close) {
    using namespace std;
    minstd_rand R;

    struct vol {
        double volume;
        unsigned long count;

        vol(double volume_): volume(volume_) {}
    };

    vector<vol> volumes = {1.2, 1.9, 3, 2.3, 9.7, 2, 0.4, 5};

    double total_volume = 0;
    for (const auto &v: volumes) total_volume += v.volume;
 
    double test_x[] = { 3000, 3000.3, 7, 0 };

    for (auto x: test_x) {
        distribute_quantity(x, volumes.begin(), volumes.end(),
                [](vol &v) -> double  { return v.volume; }, // weight
                [](vol &v,unsigned c) { v.count=c; },  // setter
                [](vol &v,int c)      { v.count+=c; },  // incrementer
                R);

        unsigned total_count = 0;
        for (auto &v: volumes) {
            double min_count = floor(floor(x)*v.volume/total_volume);
            double max_count = ceil(ceil(x)*v.volume/total_volume);

            EXPECT_GE(v.count, min_count);
            EXPECT_LE(v.count, max_count);
            total_count += v.count;
        }

        EXPECT_GE(total_count, floor(x));
        EXPECT_LE(total_count, ceil(x));
    }
}
