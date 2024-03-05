#pragma once

#include <unordered_map>

#include "strong_id.hpp"

namespace steps::util {

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

    explicit unique_indexer(Out out_, Hasher hash = Hasher())
        : lookup(0, hash)
        , out(out_)
        , count(0) {}

    unique_indexer(unique_indexer&&) = default;

    /** Return index of item x, allocating a ne~w index if required. */
    index_t operator[](const V& x) {
        auto p = lookup.find(x);
        if (p != lookup.end())
            return p->second;

        *out++ = deref_strongid(x);
        lookup.emplace(x, count);
        return count++;
    }

    /** Allocate a new index for x if it is new. */
    void insert(const V& x) {
        (void) (*this)[x];
    }

    /** Insert items given by iterator range [b,e).*/
    template <typename In>
    void insert(In b, In e) {
        while (b != e)
            insert(*b++);
    }

    /** Return number of assigned indices. */
    index_t size() const {
        return count;
    }
};

/** Construct unique_indexer given value type and output iterator */

template <typename V, typename Out, typename H = util::fnv_hash<V>>
static unique_indexer<V, H, Out> make_unique_indexer(Out out, H hasher = H()) {
    return unique_indexer<V, H, Out>(out, hasher);
}

}  // namespace steps::util