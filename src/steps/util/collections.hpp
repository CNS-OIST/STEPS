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


#ifndef STEPS_UTIL_COLLECTIONS_HPP
#define STEPS_UTIL_COLLECTIONS_HPP 1

/** \file Interface to common functionality across geom classes
 */

#include <functional>
#include <iterator>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "steps/geom/fwd.hpp"
#include "steps/util/fnv_hash.hpp"
#include "steps/util/type_traits.hpp"
#include "steps/util/strong_id.hpp"

namespace steps {
namespace util {

namespace impl {
    template <typename H>
    struct hash_ref {
        hash_ref(H hasher_): hasher(hasher_) {}
        H hasher;
        template <typename V> size_t operator()(std::reference_wrapper<V> v) const { return hasher(v.get()); }
    };

    struct equal_to_ref {
        template <typename V>
        size_t operator()(const std::reference_wrapper<V> &u,const std::reference_wrapper<V> &v) const {
            return u.get()==v.get();
        }
    };

    template <typename Table, typename C1, typename C2, typename H>
    std::vector<bool> inline map_membership(const C1 &items, const C2 &collection, const H &hasher) {
        Table hash_table(0,hasher);
        hash_table.insert(std::begin(collection),std::end(collection));

        std::vector<bool> r(container_traits<C1>::size(items));
        size_t i=0;
        for (const auto &x: items) r[i++] = hash_table.count(x)>0;

        return r;
    }
}

struct hash_references_tag {};

/** Return bit vector representing set membership.
 *
 * \param tag        (Optionally) specify hash table should hash references, not values,
 * \param items      Collection of items to test for membership.
 * \param collection Collection to test membership against.
 * \param hash       (Optionally) specify hash function for elements of collection.
 * \return           Vector v of size items.size(), where v[i] is
 *                   true if ith element of items exists in collection.
 *
 * When hash_references is specified, items and collection must have the same value_type,
 * that is, no conversions between their value_types will be deduced.
 */

template <typename C1, typename C2, typename H>
std::vector<bool> inline map_membership(hash_references_tag, const C1 &items, const C2 &collection, const H &hasher) {
    typedef typename container_traits<C2>::value_type V;
    typedef std::unordered_set<std::reference_wrapper<const V>, impl::hash_ref<H>, impl::equal_to_ref> hash_table_type;

    return impl::map_membership<hash_table_type>(items, collection, hasher);
}

template <typename C1, typename C2>
std::vector<bool> inline map_membership(hash_references_tag, const C1 &items, const C2 &collection) {
    typedef typename std::hash<typename container_traits<C2>::value_type> H;
    typedef typename container_traits<C2>::value_type V;
    typedef std::unordered_set<std::reference_wrapper<const V>, impl::hash_ref<H>, impl::equal_to_ref> hash_table_type;

    return impl::map_membership<hash_table_type>(items, collection, H());
}

template <typename C1, typename C2, typename H = std::hash<typename container_traits<C2>::value_type>>
std::vector<bool> inline map_membership(const C1 &items, const C2 &collection, const H &hasher) {
    typedef typename container_traits<C2>::value_type value_type;
    typedef std::unordered_set<value_type, H> hash_table_type;

    return impl::map_membership<hash_table_type>(items, collection, hasher);
}

template <typename C1, typename C2>
std::vector<bool> inline map_membership(const C1 &items, const C2 &collection) {
    typedef typename std::hash<typename container_traits<C2>::value_type> H;
    typedef typename container_traits<C2>::value_type value_type;
    typedef std::unordered_set<value_type, H> hash_table_type;

    return impl::map_membership<hash_table_type>(items, collection, H());
}

/** Build vector generically from any collection
 *
 * \param items   Collection of items
 * \return        std::vector containing copies of elements of items.
 */

template <typename C>
inline std::vector<typename container_traits<C>::value_type> as_vector(const C &items) {
    return std::vector<typename container_traits<C>::value_type>(std::begin(items),std::end(items));
}


/** Maintain a sequentially indexed collection of unique items.
 *
 * \tparam V      Type of items.
 * \tparam Hasher Type of hashing function.
 * \tparam Out    Type of output iterator.
 *
 * Items that are inserted are allocated an index sequentially
 * if they have not already been inserted. Uniquely inserted items
 * are appended to the supplied output iterator.
 */

template <typename V, typename Hasher, typename Out>
struct unique_indexer {
    std::unordered_map<V, index_t, Hasher> lookup;
    Out out;
    index_t count;

    explicit unique_indexer(Out out_, Hasher hash = Hasher()):
        lookup(0, hash), out(out_), count(0) {}

    unique_indexer(unique_indexer&&) = default;

    /** Return index of item x, allocating a new index if required. */
    index_t operator[](const V &x) {
        auto p = lookup.find(x);
        if (p != lookup.end()) return p->second;

        *out++ = deref_strongid(x);
        lookup.insert(std::pair<V,size_t>(x,count));
        return count++;
    }

    /** Allocate a new index for x if it is new. */
    void insert(const V &x) { (void)(*this)[x]; }

    /** Insert items given by iterator range [b,e).*/
    template <typename In>
    void insert(In b,In e) { while (b !=e) insert(*b++); }

    /** Return number of assigned indices. */
    index_t size() const { return count; }
};

/** Construct unique_indexer given value type and output iterator */

template <typename V, typename Out, typename H = steps::util::fnv_hash<V>>
static unique_indexer<V, H, Out> make_unique_indexer(Out out, H hasher = H()) {
    return unique_indexer<V,H,Out>(out, hasher);
}

}} // namespace steps::util

#endif // ndef STEPS_UTIL_COLLECTIONS_HPP


