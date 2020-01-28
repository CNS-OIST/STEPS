/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2020 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   

 */


#ifndef STEPS_UTIL_DISTRIBUTE_HPP
#define STEPS_UTIL_DISTRIBUTE_HPP

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <random>

// logging
#include <easylogging++.h>

#include "steps/error.hpp"
#include "steps/math/sample.hpp"

namespace steps {
namespace util {

/** Refence wrapper utility class for distribute_quantity implementation. */

template <typename T>
struct boxed_reference {
    T *p=nullptr;

    boxed_reference() {}
    boxed_reference(const boxed_reference &x): p(x.p) {}

    operator bool() const { return p; }
    T &get() const { return *p; }

    boxed_reference &operator=(const boxed_reference &x) {
        p=x.p;
        return *this;
    }

    boxed_reference &operator=(T &t) {
        p=&t;
        return *this;
    }

    void reset() { p=0; }
};


/** Fair stochastic distribution of quantity with weighted sampling.
 *
 * Takes a quantity, a collection to distribute over, a weight function,
 * count setter and count incrementer for items in the collection, a RNG,
 * and an optional total weight.
 *
 * Operates in two passes:
 * 1. Compute fractional allocations per item, and set the rounded-down
 *    value for each item.
 * 2. Allocate the remainder via fair sampling.
 *
 * The functional arguments should have signatures compatible with the following:
 *
 *     double weight(Item);
 *     void set_count(Item, uint);
 *     void inc_count(Item, int);
 *
 * where Item is the FwdIter reference type.
 */

template <typename FwdIter, typename Weight, typename SetCount, typename IncCount, typename Rng>
void distribute_quantity(double x, FwdIter b, FwdIter e, Weight weight, SetCount set_count, IncCount inc_count, Rng &g, double total_weight=0)
{
    static std::uniform_real_distribution<double> U;

    if (b==e) return;

    if (x<0) ArgErrLog("negative quantity to distribute");
    if (x==0) {
        // Everybody gets zero!
        for (auto i=b; i!=e; ++i) set_count(*i,0);
        return;
    }

    if (total_weight==0) 
        for (auto i=b; i!=e; ++i) total_weight += weight(*i);

    if (total_weight<=0)
        ArgErrLog("non-positive total weight for distribution");

    // Allocate rounded-down fractions and determine weights for sampling.
    if (U(g) < x-std::floor(x))
        x = 1+std::floor(x);
    else
        x = std::floor(x);

    double x_o_total = x/total_weight;
    uint allocated = 0;

    std::vector<double> pi;
    for (auto i=b; i!=e; ++i) {
        double xi = x_o_total*weight(*i);
        double xi_floor = std::floor(xi);
        pi.push_back(xi-xi_floor);

        if (xi_floor-1>std::numeric_limits<uint>::max())
            ArgErrLog("quantity too large to distribute (integer limit)");

        uint ni = static_cast<uint>(xi_floor);
        set_count(*i, ni);
        allocated += ni;
    }
    
    if (allocated>x)
        ProgErrLog("internal error in count rounding");
    auto remainder = static_cast<uint>(x-allocated);

    if (remainder == 0) return;

    // Use fractional parts as weights for sampling round.
    steps::math::adjusted_pareto_sampler<double> S(remainder,pi.begin(),pi.end());

    using item_type=typename std::remove_reference<decltype(*b)>::type;
    std::vector<boxed_reference<item_type>> extra(remainder);
    S(b,e,extra.begin(),g);

    for (auto c: extra) inc_count(c.get(),1);
}


}} // namespace steps::util


#endif // ndef STEPS_UTIL_DISTRIBUTE_HPP
