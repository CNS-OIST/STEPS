#include "math/sample.hpp"
#include "util/distribute.hpp"

#include <iterator>
#include <random>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace steps::math;
using namespace steps::util;

TEST_CASE("Sample_adjpareto_param") {
    double weights[] = {0.1, 0.2, 0.3, 0.5, 0.8, 0.7, 0.7, 0.8, 0.9};

    adjusted_pareto_sampler<double> S1(5, std::begin(weights), std::end(weights));
    REQUIRE(5u == S1.size());

    adjusted_pareto_sampler<double> S2;
    REQUIRE(0u == S2.size());

    S2.param(S1.param());
    REQUIRE(S1.size() == S2.size());
}

static void test_adjpareto_weights(const std::vector<double>& weights, const std::string& pop) {
    uint nbSampled = std::count_if(weights.begin(), weights.end(), [](double w) { return w > 0; });

    std::minstd_rand R;

    adjusted_pareto_sampler<double> S(weights.size(), weights.begin(), weights.end());

    char* sample = new char[nbSampled];
    size_t n = S(pop.begin(), pop.end(), sample, R);

    REQUIRE(nbSampled == n);

    std::vector<uint> sampleCounts(weights.size(), 0);
    for (uint i = 0; i < nbSampled; ++i) {
        auto pos = pop.find(sample[i]);
        if (pos == std::string::npos) {
            FAIL("Found sample item not from population");
        } else {
            sampleCounts[pos]++;
        }
    }

    for (uint i = 0; i < weights.size(); ++i) {
        REQUIRE_THAT(weights[i] > 0 ? 1u : 0u,
                     Catch::Matchers::WithinULP(static_cast<double>(sampleCounts[i]), 4));
    }

    delete[] sample;
}

TEST_CASE("adjpareto") {
    test_adjpareto_weights({1, 1, 1, 1}, "abcd");
    test_adjpareto_weights({0, 0, 1, 1}, "abcd");
    test_adjpareto_weights({1, 1, 0, 0}, "abcd");
    test_adjpareto_weights({0, 0}, "ab");
    test_adjpareto_weights({1}, "a");
    test_adjpareto_weights({0}, "a");
    test_adjpareto_weights({}, "");
}

TEST_CASE("Sample_adjpareto_fair") {
    std::minstd_rand R;

    std::vector<int> population = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::vector<double> weights = {0.2, 0.2, 0.2, 0.2, 0.2, 0.6, 0.6, 0.6, 0.6, 0.6};

    std::vector<int> sample(4);
    std::vector<double> histogram(population.size());

    adjusted_pareto_sampler<double> S(sample.size(), std::begin(weights), std::end(weights));

    int N = 10000;
    for (int i = 0; i < N; ++i) {
        S(population.begin(), population.end(), sample.begin(), R);
        for (auto j: sample)
            ++histogram.at(j);
    }

    for (auto& x: histogram)
        x /= N;

    // Check that each observed mean is within (say) 4 standard deviations of
    // inclusion probability pi. sd = sqrt(pi(1-pi)/N)

    for (size_t i = 0; i < histogram.size(); ++i) {
        double z = std::abs(histogram[i] - weights[i]) /
                   std::sqrt(weights[i] * (1 - weights[i]) / N);
        REQUIRE(z < 4);
    }
}

template <typename T>
struct vol {
    double volume;
    T count;

    vol(double volume_)
        : volume(volume_) {}
};

template <typename Quantity>
void test_distribute_quantity_approx(const std::vector<double> quantities) {
    using volumes_type = std::vector<vol<Quantity>>;
    std::vector<volumes_type> volumes_cases{{0, 1.2, 1.9, 3, 2.3, 9.7, 2, 0.4, 5}, {}};

    for (auto& volumes: volumes_cases) {
        double total_volume = 0;
        for (const auto& v: volumes) {
            total_volume += v.volume;
        }

        std::minstd_rand R;
        auto weight = [](const typename volumes_type::iterator& v) -> double { return v->volume; };
        auto set_count = [](typename volumes_type::value_type& v, Quantity c) { v.count = c; };
        auto inc_count = [](typename volumes_type::value_type& v) { ++v.count; };

        for (auto x: quantities) {
            distribute_quantity<Quantity>(
                x, volumes.begin(), volumes.end(), weight, set_count, inc_count, R);

            Quantity total_count = 0;
            for (auto& v: volumes) {
                double min_count = floor(floor(x) * v.volume / total_volume);
                double max_count = ceil(ceil(x) * v.volume / total_volume);

                REQUIRE(v.count >= min_count);
                REQUIRE(v.count <= max_count);
                total_count += v.count;
            }
            if (volumes.size() > 0) {
                REQUIRE(total_count >= floor(x));
                REQUIRE(total_count <= ceil(x));
            }
        }
    }
}

TEST_CASE("Sample_distribute_quantity_approx") {
    test_distribute_quantity_approx<uint>({3000, 3000.3, 7, 0});
    test_distribute_quantity_approx<std::int64_t>(
        {3000, 3000.3, 7, 0, static_cast<double>(std::numeric_limits<uint>::max()) * 1e7});
}
