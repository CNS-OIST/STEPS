#pragma once
#include <algorithm>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>


namespace zee {

//-----------------------------------------------

template <typename Num, typename T>
struct MultinomialDistribution {
    MultinomialDistribution(Num n, const std::vector<T>& rates)
        : n_(n)
        , rates_(rates)
        , ret_(rates.size()) {
        if (!std::all_of(rates.begin(), rates.end(), [](auto v) { return v >= 0; })) {
            throw std::logic_error("Negative rates encountered.");
        }
    }


    template <typename RNG>
    const std::vector<Num>& operator()(RNG& rng) {
        size_t rem_n(n_);
        sum_rates_ = std::accumulate(rates_.begin(), rates_.end(), T{}, std::plus<T>());
        for (size_t i = 0; i < rates_.size(); ++i) {
            const auto prob = rates_[i] / sum_rates_;
            if (prob < 1.0) {
                std::binomial_distribution<Num> bin(rem_n, prob);
                ret_[i] = bin(rng);
                rem_n -= ret_[i];
                sum_rates_ -= rates_[i];
            } else {
                ret_[i] = rem_n;
            }
        }
        return ret_;
    }
    Num n_;
    const std::vector<T>& rates_;
    T sum_rates_;
    std::vector<Num> ret_;
};

};  // namespace zee
