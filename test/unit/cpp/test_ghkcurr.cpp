#include <cmath>
#include <iterator>
#include <numeric>

#include "math/constants.hpp"
#include "math/ghk.hpp"

#include "gtest/gtest.h"

using steps::math::FARADAY;
using steps::math::GAS_CONSTANT;
using steps::math::GHKcurrent;
using steps::math::permeability;

constexpr double T = 310.15;
constexpr double iconc = 1e-6;
constexpr double oconc = 1e-3;

TEST(GHKCurrent, permeability) {
    // 0 membrane potential
    ASSERT_EQ(permeability(1, 0, 1, T, iconc, oconc),
              2 / (iconc * 1000.0 + oconc * 1000.0) * GAS_CONSTANT * T / (FARADAY * FARADAY));
}

TEST(GHKCurrent, current) {
    // 0 Permeablity
    ASSERT_EQ(GHKcurrent(0, 1, 1, T, iconc, oconc), 0);
    // 0 membrane potential
    ASSERT_EQ(GHKcurrent(1, 0, 1, T, iconc, oconc), FARADAY * (iconc - oconc));
    // Equal concentrations on both sides
    double P = permeability(1, 1, -1, T, 0.001, 0.001);
    ASSERT_EQ(GHKcurrent(P, 1, -1, T, 1, 1), 1);
}
