/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
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

#pragma once

/** \file Extra type traits to ease generic programming.
 */

#include <cstddef>
#include <type_traits>

namespace steps::util {

/** Test if type is scalar or a contiguous fixed-size array of
 *  scalar or array types.
 */

template <typename T>
struct is_scalar_or_array: std::is_scalar<T>::type {};

template <typename T, size_t n>
struct is_scalar_or_array<T[n]>: is_scalar_or_array<T>::type {};

template <typename T, size_t n>
struct is_scalar_or_array<std::array<T, n>>: is_scalar_or_array<T>::type {};

/** Generic value_type and size() function for standard libaray
 * containers and array types.
 */

template <typename C>
struct container_traits {
    typedef typename C::value_type value_type;
    static size_t size(const C& c) {
        return c.size();
    }
};

template <typename T, size_t n>
struct container_traits<T[n]> {
    typedef T value_type;
    static constexpr size_t size(...) {
        return n;
    }
};

template <typename T, size_t n>
struct container_traits<std::array<T, n>> {
    typedef T value_type;
    static constexpr size_t size(...) {
        return n;
    }
};

// helper constant
template <class>
inline constexpr bool always_false_v = false;

}  // namespace steps::util
