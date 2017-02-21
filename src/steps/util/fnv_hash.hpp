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


#ifndef STEPS_UTIL_FNV_HASH_HPP
#define STEPS_UTIL_FNV_HASH_HPP

#include <cinttypes>

#include "steps/util/type_traits.hpp"

/** /file Hashing support for aggregates, containers, based
 * on the Fowler-Noll-Vo FNV-1a hash function.
 */

namespace steps {
namespace util {

typedef uint64_t hash_type;

template <typename T>
inline hash_type fnv1a_combine(hash_type h, const T &v) {
    const unsigned char *k=reinterpret_cast<const unsigned char *>(&v);
    for (size_t i=0; i<sizeof(T); ++i) {
        h ^= k[i];
        h *= 0x100000001b3ull;
    }
    return h;
}

template <typename T, typename... Rest>
inline hash_type fnv1a_combine(hash_type h, const T &v, const Rest &... vs) {
    return fnv1a_combine(fnv1a_combine(h,v),vs...);
}

template <typename... T>
inline hash_type fnv1a(const T &... vs) {
    return fnv1a_combine(0xcbf29ce484222325ull, vs...);
}


namespace impl {
    template <typename T, typename enable = void>
    struct fnv_hash;

    template <typename T>
    struct fnv_hash<T, typename std::enable_if<is_scalar_or_array<T>::value>::type> {
        size_t operator()(const T &v) const { return static_cast<size_t>(fnv1a(v)); }
    };
}

/* Define fnv::hash for scalar types and arrays of scalar types by default */

template <typename T>
struct fnv_hash: impl::fnv_hash<T> {};

}} // namespace steps::util


#endif // ndef STEPS_UTIL_FNV_HASH_HPP
