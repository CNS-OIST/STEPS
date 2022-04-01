#pragma once

#include <algorithm>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

#include "rng/rng.hpp"
#include "tools.hpp"
#include "util/collections.hpp"
#include "util/distribute.hpp"

namespace steps {
namespace math {

/// TODO make this a scoped enum once STEPS requires cython 3x
enum DistributionMethod { DIST_UNIFORM, DIST_MULTINOMIAL };

//-----------------------------------------------
/**
 * Distributor of objects based on rates
 *
 * We need to distribute n objects (usually molecules) among rates_.size() options.
 *
 * Uniform_distribute distributes them in the most uniform way possible weighted with
 * rates_[i]/sum(rates_) with fair sampling. This means that, if rates are all equal, there is at
 * most 1 molecule difference between the least full bin and the fullest bin max(counts) -
 * min(counts) <= 1
 *
 * Multinomial_distribute distributes them in the bins
 * with a probability of each bin = rates_[i]/sum(rates_) in a multinomial way
 *
 * EXPLANATION
 * In Uniform_distribute the number of objects in each bins is proportional to the associated 
 * weight and is deterministic, such that final distribution among the bins is homogeneous.
 * In Multinomial_distribute the number of objects in each bins is randomly sampled 
 * from a multinomial distribution.
 *
 * distribute() produces the distribution
 * 
 * \tparam Num
 * \tparam Container
 */
template <typename Num, typename Container>
class Distribution {
  public:
    using value_type = typename util::sequence_container_traits<Container>::value_type;
    using size_type = typename util::sequence_container_traits<Container>::size_type;
    template <typename T>
    using write_type = typename util::sequence_container_traits<Container>::template write_type<T>;

    Distribution(Num n, const Container& rates) noexcept
        : n_(n)
        , rates_(rates)
        , sum_rates_(std::accumulate(rates_.begin(), rates_.end(), value_type{}, std::plus<>()))
        , ret_(rates.size(), 0) {
        assert(std::all_of(rates.begin(), rates.end(), [](auto v) { return v >= 0; }));
        // If all the rates are 0 it is possible that you are dealing with a mesh that is not
        // supported. For example, a patch in between compartments is not supported in steps 4
        assert(!rates.size() || sum_rates_ > 0);
    }

    Distribution(Num n, const Container& rates, write_type<Num>& output) noexcept
        : n_(n)
        , rates_(rates)
        , sum_rates_(std::accumulate(rates_.begin(), rates_.end(), value_type{}, std::plus<>()))
        , ret_(output) {
        assert(std::all_of(rates.begin(), rates.end(), [](auto v) { return v >= 0; }));
    }

    /// selector of distribution method
    template <typename RNG>
    const write_type<Num>& distribute(RNG& rng, const DistributionMethod distribution) {
        switch (distribution) {
        case DistributionMethod::DIST_UNIFORM:
            return distribute_uniform(rng);
        case DistributionMethod::DIST_MULTINOMIAL:
            return distribute_multinomial(rng);
        default:
            throw std::logic_error("Unknown distribution method: " +
                                   std::to_string(static_cast<int>(distribution)));
        }
    }

  private:
    /// Distribute objects in bins in multinomial way
    template <typename RNG>
    const write_type<Num>& distribute_multinomial(RNG& rng) {
        Num rem_n(this->n_);
        value_type sum_rates = this->sum_rates_;
        for (size_type i = 0; i < this->rates_.size(); ++i) {
            const auto prob = this->rates_[i] / sum_rates;
            if (prob < 1.0) {
                if constexpr (std::is_same_v<RNG, steps::rng::RNG>) {
                    this->ret_[i] = static_cast<Num>(
                        rng.getBinom(static_cast<uint>(rem_n), static_cast<double>(prob)));
                } else {
                    std::binomial_distribution<Num> bin(rem_n, prob);
                    this->ret_[i] = bin(rng);
                }
                rem_n -= this->ret_[i];
                sum_rates -= this->rates_[i];
            } else {
                this->ret_[i] = rem_n;
            }
        }
        return this->ret_;
    }

    /// Distribute objects in bins in uniform way
    template <typename RNG>
    const write_type<Num>& distribute_uniform(RNG& rng) {
        if (this->rates_.size() == 0) {
            return this->ret_;
        }

        auto set_count = [](Num& tet, uint c) { tet = c; };
        auto inc_count = [](Num& tet, int c) { tet += c; };

        auto weight = [this](const auto tet) {
            // we need to use this way for finding the distance because osh::Write
            // does not implement real iterators
            const size_t ii = static_cast<size_t>(&(*tet) - this->ret_.data());
            return this->rates_[ii];
        };

        steps::util::distribute_quantity(this->n_,
                                         this->ret_.begin(),
                                         this->ret_.end(),
                                         weight,
                                         set_count,
                                         inc_count,
                                         rng,
                                         this->sum_rates_);

        return this->ret_;
    }

  private:
    /// Number of objects (molecules) to be distributed
    const Num n_;
    /// Probabilities as: rates_[i]/sum_rates_
    const Container& rates_;
    /// Total sum of rates_
    const value_type sum_rates_;
    /// Distribution. Valid only after operator()
    write_type<Num> ret_;
};

template <typename Num, typename Container>
Distribution<Num, Container> make_dist(Num n, const Container& rates) noexcept {
    return {n, rates};
}

template <typename Num, typename Container, typename Output>
Distribution<Num, Container> make_dist(Num n, const Container& rates, Output& output) noexcept {
    return {n, rates, output};
}


} // namespace math
} // namespace steps
