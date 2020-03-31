#pragma once

#include <cmath>
#include <random>

namespace zee {

/**
 * Function that rounds a double value to one of the closest straddling integers with a probability
 * dependent on the proximity
 * \tparam T the type of the integer returned
 * \tparam RNG random number generator
 * \param value to round to one of the closest integer
 */

template <typename T, class RNG>
T stochastic_round(double value, RNG& rng) {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const auto floored_value = std::floor(value);
    return static_cast<T>(floored_value) + ((value - floored_value) > uniform(rng) ? 1 : 0);
}

/**
 * Function that rounds a double value to one of the closest straddling integers with a probability
 * dependent on the proximity
 * \tparam T the type of the integer returned
 * \tparam RNG random number generator
 * \param value to round to one of the closest integer
 * \param upper_bound is the maximum value allowed for the rounded value
 */

template <typename T, class RNG>
T stochastic_round(double value, RNG& rng, T upper_bound) {
    return std::min(upper_bound, stochastic_round<T>(value, rng));
}

}  // namespace zee
