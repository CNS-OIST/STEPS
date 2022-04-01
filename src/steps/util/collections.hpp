/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
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

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <type_traits>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#include <boost/bind/bind.hpp>
#include <boost/functional/value_factory.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <gsl>

#ifdef STEPS_USE_DIST_MESH
#include <Omega_h_array.hpp>
#endif // STEPS_USE_DIST_MESH

#include "fnv_hash.hpp"
#include "strong_id.hpp"
#include "type_traits.hpp"
#include "util/vocabulary.hpp"

using namespace boost::placeholders;

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
} // namespace impl

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

/** Check if two given floating point values almost equal.
 * This function checks if the given floating point values x and y are almost equal, with 4 ulps as default precision
 * The choice of 4 ulps is according to googletest, but may need adjustment for special cases
 * https://github.com/google/googletest/blob/472cd8fd8b1c665bddfd021ad0f62d6747fe8e72/googletest/include/gtest/internal/gtest-internal.h#L290
 *
 * \param x       the first value for comaprison
 * \param y       the second value for comaprison
 * \param ulp     the desired precision in ULPs (units in the last place)
 * \return        true of x and y are almost equal, false if they are not.
 */

template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
    almost_equal(T x, T y, int ulp = 4)
{
    return std::fabs(x-y) <= std::numeric_limits<T>::epsilon() * std::fabs(x+y) * ulp
        || std::fabs(x-y) < std::numeric_limits<T>::min();
}

namespace impl {

// sink to consume expanded arguments
struct sink {
    template <typename T>
    sink(T const& /*v*/) {}
};

template <int... Values>
struct integer_sequence {};

template <std::size_t Size, int... Accu>
struct ones_traits {
    using type = typename ones_traits<Size - 1, 1, Accu...>::type;
};

template <int... Accu>
struct ones_traits<0, Accu...> {
    using type = integer_sequence<Accu...>;
};

template <typename T, int... Ones>
std::initializer_list<T> initializer_list_impl(T value, const integer_sequence<Ones...>&) {
    static auto array = {(sink{Ones}, value)...};
    return array;
}

}  // namespace impl

template <std::size_t Size, typename T>
std::initializer_list<T> initializer_list(T value) {
    return impl::initializer_list_impl<T>(value, typename impl::ones_traits<Size>::type());
}

/**
 * Generate a class providing \c begin() and \c end() methods
 * to iterate over instances of type \c Entity accepting an integer in
 * constructor
 * @tparam Entity Type manipulated by the iterators
 */
template <typename Entity, typename Int = int> struct EntityIterator {
  /// Construct a range from 0 to \a num_entities
  explicit EntityIterator(Int num_entities) noexcept
      : num_entities_(num_entities) {}

  /// get an iterator at the beginning of the range providing \c Entity(Int{})
  auto begin() const noexcept {
    const auto ctor = boost::bind(boost::value_factory<Entity>(), _1);
    return boost::make_transform_iterator(boost::make_counting_iterator(Int{}),
                                          ctor);
  }

  /// get an iterator at the end of the range providing \c Entity(num_entities)
  auto end() const noexcept {
    const auto ctor = boost::bind(boost::value_factory<Entity>(), _1);
    return boost::make_transform_iterator(
        boost::make_counting_iterator(num_entities_), ctor);
  }

private:
  /// Right bound of the range
  const Int num_entities_;
};

/**
 * A type traits useful to write generic code using either
 * a std::vector or Omega_h::Read<T> arrays.
 */
template <typename Container> struct sequence_container_traits {};

template <typename T> struct sequence_container_traits<std::vector<T>> {
  using container_type = std::vector<T>;
  using value_type = typename container_type::value_type;
  using size_type = typename container_type::size_type;
  template <typename U> using write_type = std::vector<U>;
};

template <typename T> struct sequence_container_traits<gsl::span<T>> {
  using container_type = gsl::span<T>;
  using value_type = T;
  using size_type = typename container_type::index_type;
  template <typename U> using write_type = gsl::span<U>;
};

#ifdef STEPS_USE_DIST_MESH

template <> struct sequence_container_traits<Omega_h::Reals> {
  using container_type = Omega_h::Reals;
  using value_type = Omega_h::Real;
  using size_type = Omega_h::LO;
  template <typename T> using write_type = Omega_h::Write<T>;
};

template <> struct sequence_container_traits<Omega_h::LOs> {
  using container_type = Omega_h::LOs;
  using value_type = Omega_h::LO;
  using size_type = Omega_h::LO;
  template <typename T> using write_type = Omega_h::Write<T>;
};

template <> struct sequence_container_traits<Omega_h::GOs> {
  using container_type = Omega_h::GOs;
  using value_type = Omega_h::GO;
  using size_type = Omega_h::LO;
  template <typename T> using write_type = Omega_h::Write<T>;
};

template <typename T> Omega_h::Read<T> createRead(Omega_h::Write<T> array) {
  return {array};
}
#endif // STEPS_USE_DIST_MESH

} // namespace util
} // namespace steps

#endif // ndef STEPS_UTIL_COLLECTIONS_HPP
