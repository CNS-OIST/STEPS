/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
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

#ifndef STEPS_MATH_SMALLSORT_HPP
#define STEPS_MATH_SMALLSORT_HPP 1

/// Inline sorting of small, fixed-size
/// random-access sequences.

#include <utility>
#include <algorithm>
#include <iostream>

namespace steps {
namespace math {

namespace impl {
    // compare-and-sort two elements of a:

    template <int m,int n>
    struct S {
        template <typename A>
        [[gnu::always_inline]] static void run(A &a) {
            if (a[m]>a[n]) std::swap(a[m],a[n]);
        }
    };

    // minimal networks up to 6, then recusive merge

    template <int n,typename A>
    struct small_sort_inplace {
        [[gnu::always_inline]] static void run(A &a) {
            std::sort(std::begin(a),std::end(a));
        }
    };

    template <typename A>
    struct small_sort_inplace<0,A> {
        static void run(A &a) {}
    };

    template <typename A>
    struct small_sort_inplace<1,A> {
        static void run(A &a) {}
    };

    template <typename A>
    struct small_sort_inplace<2,A> {
        [[gnu::always_inline]] static void run(A &a) { S<0,1>::run(a); }
    };

    template <typename A>
    struct small_sort_inplace<3,A> {
        [[gnu::always_inline]] static void run(A &a) {
            S<0,1>::run(a);
            S<1,2>::run(a);
            S<0,1>::run(a);
        }
    };

    template <typename A>
    struct small_sort_inplace<4,A> {
        [[gnu::always_inline]] static void run(A &a) {
            S<0,1>::run(a);
            S<2,3>::run(a);
            S<0,2>::run(a);
            S<1,3>::run(a);
            S<1,2>::run(a);
        }
    };

    template <typename A>
    struct small_sort_inplace<5,A> {
        [[gnu::always_inline]] static void run(A &a) {
            S<0,1>::run(a);
            S<2,4>::run(a);
            S<0,3>::run(a);
            S<1,4>::run(a);
            S<1,2>::run(a);
            S<3,4>::run(a);
            S<0,1>::run(a);
            S<2,3>::run(a);
            S<1,2>::run(a);
        }
    };

    template <typename A>
    struct small_sort_inplace<6,A> {
        [[gnu::always_inline]] static void run(A &a) {
            S<0,1>::run(a);
            S<2,3>::run(a);
            S<4,5>::run(a);
            S<0,2>::run(a);
            S<1,4>::run(a);
            S<3,5>::run(a);
            S<0,1>::run(a);
            S<2,3>::run(a);
            S<4,5>::run(a);
            S<1,2>::run(a);
            S<3,4>::run(a);
            S<2,3>::run(a);
        }
    };

} // namespace impl

template <int n,typename A>
[[gnu::always_inline]] inline void small_sort_inplace(A &a) {
    impl::small_sort_inplace<n,A>::run(a);
}

template <int n,typename A>
[[gnu::always_inline]] inline A small_sort(A a) {
    impl::small_sort_inplace<n,A>::run(a);
    return a;
}

}} // namespace steps::math

#endif // ndef STEPS_MATH_SMALLSORT_HPP


