#pragma once

#include <iterator>
#include <stdexcept>
#include <type_traits>

#include <Omega_h_array.hpp>
#include <Omega_h_defines.hpp>

namespace zee {

/**
 * Iterator class for \a strong_ids
 * \tparam StrongIds \a strong_ids container type
 */
template <typename StrongIds>
struct strong_ids_iterator {
    using size_type = typename StrongIds::size_type;
    using value_type = typename StrongIds::value_type;
    using difference_type = size_type;
    using pointer = value_type*;
    using reference = std::false_type;
    using iterator_category = std::random_access_iterator_tag;

    inline strong_ids_iterator(const StrongIds& strong_ids, size_type pos) noexcept;
    inline strong_ids_iterator(const strong_ids_iterator& other) noexcept;

    inline strong_ids_iterator& operator++() noexcept;
    inline strong_ids_iterator operator++(int) noexcept;

    inline value_type operator*() const noexcept;

    inline bool operator==(const strong_ids_iterator& other) const noexcept;
    inline bool operator!=(const strong_ids_iterator& other) const noexcept;

  private:
    const StrongIds& strong_ids_;
    size_type pos_;
};

/**
 * \brief Read only random-access container over strong identifiers
 *
 * Rely on Omega_h::Read container internally so only the following integral types are supported:
 * I8, I32, I64, Real
 * \tparam StrongId Concrete zee:strong_id<...> type
 */
template <typename StrongId>
struct strong_ids {
    using value_type = StrongId;
    using size_type = Omega_h::LO;
    using data_type = typename StrongId::value_type;

    using this_type = strong_ids<value_type>;
    using iterator_type = strong_ids_iterator<this_type>;

    /**
     * \name Ctors & dtors
     * \{
     */

    inline strong_ids();
    inline strong_ids(const Omega_h::Read<data_type>& ids);
    inline strong_ids(const Omega_h::Write<data_type>& ids);
    inline strong_ids(std::initializer_list<data_type> ids, std::string const& name = "");

    /** \} */

    /**
     * \name Element access
     * \{
     */

    /**
     * access the element at specified location pos, with bounds checking
     * \param pos position of the element to return
     * \throw std::out_of_range if pos is not within the range of the container
     * \return element at requested location
     */
    inline value_type at(size_type pos) const;

    /**
     * access the element at specified location pos
     * \param pos position of the element to return
     * \return element at requested location
     */
    inline value_type operator[](size_type pos) const noexcept;

    /// access the first element
    inline value_type front() const noexcept;

    /// access the last element
    inline value_type back() const noexcept;

    /// direct access to the underlying array
    inline const Omega_h::Read<data_type>& data() const noexcept;

    /** \} */

    /**
     * \name Iterators
     * \{
     */

    /// \return an iterator to the beginning
    inline iterator_type begin() const noexcept;

    /// \return an iterator to the end
    inline iterator_type end() const noexcept;
    /** \} */

    /**
     * \name Capacity
     */

    /// \return the number of elements
    inline size_type size() const noexcept;

    /// \return true if the container is empty, false otherwise
    inline bool empty() const noexcept;

    /** \} */

  private:
    Omega_h::Read<data_type> ids_;
};

}  // namespace zee

/// inline class definitions

namespace zee {

// class strong_ids

template <typename StrongId>
inline strong_ids<StrongId>::strong_ids()
    : strong_ids(Omega_h::Write<data_type>(0)) {}

template <typename StrongId>
inline strong_ids<StrongId>::strong_ids(const Omega_h::Read<data_type>& ids)
    : ids_(ids) {}

template <typename StrongId>
inline strong_ids<StrongId>::strong_ids(const Omega_h::Write<data_type>& ids)
    : ids_(ids) {}

template <typename StrongId>
inline strong_ids<StrongId>::strong_ids(std::initializer_list<data_type> ids,
                                        const std::string& name)
    : ids_(ids, name) {}

template <typename StrongId>
inline typename strong_ids<StrongId>::value_type strong_ids<StrongId>::at(size_type pos) const {
    if (pos < 0 || pos >= this->ids_.size()) {
        throw std::out_of_range("Invalid access");
    }
    return value_type(this->ids_[pos]);
}

template <typename StrongId>
inline
    typename strong_ids<StrongId>::value_type strong_ids<StrongId>::operator[](size_type pos) const
    noexcept {
    return value_type(ids_[pos]);
}

template <typename StrongId>
inline typename strong_ids<StrongId>::value_type strong_ids<StrongId>::front() const noexcept {
    return value_type(ids_.first());
}

template <typename StrongId>
inline typename strong_ids<StrongId>::value_type strong_ids<StrongId>::back() const noexcept {
    return value_type(ids_.last());
}

template <typename StrongId>
inline const typename Omega_h::Read<typename StrongId::value_type>& strong_ids<StrongId>::data()
    const noexcept {
    return ids_;
}

template <typename StrongId>
inline typename strong_ids<StrongId>::iterator_type strong_ids<StrongId>::begin() const noexcept {
    return iterator_type(*this, 0);
}

template <typename StrongId>
inline typename strong_ids<StrongId>::iterator_type strong_ids<StrongId>::end() const noexcept {
    return iterator_type(*this, size());
}

template <typename StrongId>
inline typename strong_ids<StrongId>::size_type strong_ids<StrongId>::size() const noexcept {
    return ids_.size();
}

template <typename StrongId>
inline bool strong_ids<StrongId>::empty() const noexcept {
    return size() == 0;
}

// class strong_ids_iterator

template <typename StrongIds>
inline strong_ids_iterator<StrongIds>::strong_ids_iterator(const StrongIds& strong_ids,
                                                           size_type pos) noexcept
    : strong_ids_(strong_ids)
    , pos_(pos) {}

template <typename StrongIds>
inline strong_ids_iterator<StrongIds>::strong_ids_iterator(
    const strong_ids_iterator& other) noexcept
    : strong_ids_(other.strong_ids_)
    , pos_(other.pos_) {}

template <typename StrongIds>
inline strong_ids_iterator<StrongIds>& strong_ids_iterator<StrongIds>::operator++() noexcept {
    ++pos_;
    return *this;
}

template <typename StrongIds>
inline strong_ids_iterator<StrongIds> strong_ids_iterator<StrongIds>::operator++(int) noexcept {
    strong_ids_iterator copy(*this);
    ++pos_;
    return copy;
}

template <typename StrongIds>
inline
    typename strong_ids_iterator<StrongIds>::value_type strong_ids_iterator<StrongIds>::operator*()
        const noexcept {
    return strong_ids_[pos_];
}

template <typename StrongIds>
inline bool strong_ids_iterator<StrongIds>::operator==(const strong_ids_iterator& other) const
    noexcept {
    return &strong_ids_ == &other.strong_ids_ && pos_ == other.pos_;
}

template <typename StrongIds>
inline bool strong_ids_iterator<StrongIds>::operator!=(const strong_ids_iterator& other) const
    noexcept {
    return pos_ != other.pos_ || &strong_ids_ != &other.strong_ids_;
}

}  // namespace zee
