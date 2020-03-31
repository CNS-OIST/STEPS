/**
 * \file Provides implementation of flat_multimap data structure.
 */
#pragma once

#include <iterator>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include <Omega_h_array.hpp>
#include <gsl/gsl-lite.hpp>

namespace zee {

/// tag representing the Standard Template Library
class STL {};
/// tag representing Omega_h
class OSH {};

namespace traits {

/**
 * Data trait providing types used to manipulate a dynamic random-access container
 * \tparam T the type of the elements
 * \tparam tag Data structure family (STL, OSH, ...)
 * \tparam Allocator Allocator class (only used by std::vector)
 */
template <class T, class tag, class Allocator = std::allocator<T>>
struct dynamic_vector {
    static_assert(true, "missing implementation for this tag");
};

template <class T, class Allocator>
struct dynamic_vector<T, STL, Allocator> {
    /// random access container for read operations only.
    using ro_type = std::vector<T, Allocator>;
    /// random access container for read and write operations.
    using rw_type = std::vector<T, Allocator>;
    /// the type of the data indices.
    using size_type = std::size_t;
};

template <class T, class Allocator>
struct dynamic_vector<T, OSH, Allocator> {
    /// random access container for read operations only.
    using ro_type = Omega_h::Read<T>;
    /// random access container for read and write operations.
    using rw_type = Omega_h::Write<T>;
    /// the type of the data indices.
    using size_type = Omega_h::LO;
};

}  // namespace traits


template <typename T, int Size, class Tag>
class flat_multimap;

namespace {

/**
 * Iterator over the values of a specific element data.
 * \tparam T stored data type
 * \tparam Size number of values per data
 * \tparam Tag data structure implementation
 */
template <typename T, int Size, class Tag>
class flat_multimap_data_iterator: public std::iterator<std::input_iterator_tag, gsl::span<T>> {
  public:
    using super_type = std::iterator<std::input_iterator_tag, gsl::span<T>>;
    using value_type = typename super_type::value_type;
    using fmm_type = typename std::conditional<
        std::is_const<T>::value,
        const flat_multimap<typename std::remove_const<T>::type, Size, Tag>,
        flat_multimap<T, Size, Tag>>::type;
    using size_type = typename fmm_type::size_type;

    inline flat_multimap_data_iterator(fmm_type& fmm, size_type element_idx) noexcept
        : data_(fmm.ab2c_.data() + fmm.a2ab_[element_idx], Size) {}

    inline flat_multimap_data_iterator(const flat_multimap_data_iterator& other) noexcept
        : data_(other.data_) {}

    inline flat_multimap_data_iterator& operator++() noexcept {
        data_ = value_type(data_.data() + Size, Size);
        return *this;
    }
    inline flat_multimap_data_iterator operator++(int) noexcept {
        flat_multimap_data_iterator copy(*this);
        data_ = value_type(data_.data() + Size, Size);
        return copy;
    }

    inline bool operator==(const flat_multimap_data_iterator& other) const noexcept {
        return data_.begin() == other.data_.begin();
    }

    inline bool operator!=(const flat_multimap_data_iterator& other) const noexcept {
        return data_.begin() != other.data_.begin();
    }

    typename super_type::value_type& operator*() noexcept {
        return data_;
    }

    const typename super_type::value_type& operator*() const noexcept {
        return data_;
    }

  private:
    value_type data_;
};

template <typename T, class Tag>
class flat_multimap_data_iterator<T, 1, Tag>: public std::iterator<std::input_iterator_tag, T> {
  public:
    using super_type = std::iterator<std::input_iterator_tag, T>;
    using value_type = typename super_type::value_type;
    using fmm_type =
        typename std::conditional<std::is_const<T>::value,
                                  const flat_multimap<typename std::remove_const<T>::type, 1, Tag>,
                                  flat_multimap<T, 1, Tag>>::type;
    using size_type = typename fmm_type::size_type;


    inline flat_multimap_data_iterator(fmm_type& fmm, size_type element_idx)
        : data_(fmm.ab2c_.data() + fmm.a2ab_[element_idx]) {}

    inline flat_multimap_data_iterator(const flat_multimap_data_iterator& other)
        : data_(other.data_) {}

    inline flat_multimap_data_iterator& operator++() noexcept {
        ++data_;
        return *this;
    }
    inline flat_multimap_data_iterator operator++(int) noexcept {
        flat_multimap_data_iterator copy(*this);
        ++data_;
        return copy;
    }

    inline bool operator==(const flat_multimap_data_iterator& other) const noexcept {
        return data_ == other.data_;
    }

    inline bool operator!=(const flat_multimap_data_iterator& other) const noexcept {
        return data_ != other.data_;
    }

    inline typename super_type::value_type& operator*() const noexcept {
        return *data_;
    }

  private:
    T* data_;
};


template <typename T, int Size, class Tag>
class flat_multimap_element {
  public:
    using iterator = flat_multimap_data_iterator<T, Size, Tag>;
    using const_iterator = flat_multimap_data_iterator<const T, Size, Tag>;
    using fmm_type = typename std::conditional<
        std::is_const<T>::value,
        const flat_multimap<typename std::remove_const<T>::type, Size, Tag>,
        flat_multimap<T, Size, Tag>>::type;
    using size_type = typename fmm_type::size_type;

    inline flat_multimap_element(fmm_type& fmm, size_type element_idx) noexcept
        : fmm_(fmm)
        , element_idx_(element_idx) {}

    /**
     * \name Iterators
     * \{
     */
    /// returns an iterator to the beginning
    inline iterator begin() noexcept {
        return {fmm_, element_idx_};
    }
    /// returns an iterator to the beginning
    inline const_iterator begin() const noexcept {
        return {fmm_, element_idx_};
    }
    /// returns a constant iterator to the beginning
    inline const_iterator cbegin() const noexcept {
        return {fmm_, element_idx_};
    }
    /// returns an iterator to the end
    inline iterator end() noexcept {
        return {fmm_, element_idx_ + 1};
    }
    /// returns an iterator to the end
    inline const_iterator end() const noexcept {
        return {fmm_, element_idx_ + 1};
    }
    /// returns a constant iterator to the end
    inline const_iterator cend() const noexcept {
        return {fmm_, element_idx_ + 1};
    }
    /** \} */

    /**
     * \name Capacity
     * \{
     */

    /**
     *
     * \return number of data for the current element
     */
    inline size_type size() const noexcept {
        return fmm_.size(element_idx_);
    }

    inline bool empty() const noexcept {
        return size() == 0;
    }

    /**
     * \}
     */

    inline flat_multimap_element& operator++() noexcept {
        ++element_idx_;
        return *this;
    }

    inline flat_multimap_element operator++(int) noexcept {
        flat_multimap_element copy(*this);
        ++element_idx_;
        return copy;
    }

    inline bool operator==(const flat_multimap_element& other) const noexcept {
        return &fmm_ == &other.fmm_ && element_idx_ == other.element_idx_;
    }

    inline bool operator!=(const flat_multimap_element& other) const noexcept {
        return &fmm_ != &other.fmm_ || element_idx_ != other.element_idx_;
    }

  private:
    fmm_type& fmm_;
    size_type element_idx_;
};

template <typename T, int Size, class Tag>
class flat_multimap_element_iterator
    : public std::iterator<std::input_iterator_tag, flat_multimap_element<T, Size, Tag>> {
  public:
    using super_type = std::iterator<std::input_iterator_tag, flat_multimap_element<T, Size, Tag>>;
    using value_type = typename super_type::value_type;
    using fmm_type = typename std::conditional<
        std::is_const<T>::value,
        const flat_multimap<typename std::remove_const<T>::type, Size, Tag>,
        flat_multimap<T, Size, Tag>>::type;
    using size_type = typename fmm_type::size_type;

    inline flat_multimap_element_iterator(fmm_type& fmm, size_type element_idx) noexcept
        : element_(fmm, element_idx) {}

    inline flat_multimap_element_iterator(const flat_multimap_element_iterator& other) noexcept
        : element_(other.element_) {}

    inline flat_multimap_element_iterator& operator++() noexcept {
        ++element_;
        return *this;
    }

    inline flat_multimap_element_iterator operator++(int) noexcept {
        const flat_multimap_element_iterator copy(*this);
        ++element_;
        return copy;
    }

    inline bool operator==(const flat_multimap_element_iterator& other) const noexcept {
        return element_ == other.element_;
    }

    inline bool operator!=(const flat_multimap_element_iterator& other) const noexcept {
        return element_ != other.element_;
    }

    inline const value_type& operator*() const noexcept {
        return element_;
    }

    inline value_type& operator*() noexcept {
        return element_;
    }

  private:
    /// iterator instance
    value_type element_;
};

template <typename T, int Size, class Tag>
struct flat_multimap_item_accessor {
    using fmm_type = typename std::conditional<
        std::is_const<T>::value,
        const flat_multimap<typename std::remove_const<T>::type, Size, Tag>,
        flat_multimap<T, Size, Tag>>::type;
    using size_type = typename fmm_type::size_type;
    using value_type = typename flat_multimap_data_iterator<T, Size, Tag>::value_type;
    using const_value_type = typename flat_multimap_data_iterator<const T, Size, Tag>::value_type;

    inline value_type operator()(fmm_type& fmm, size_type element_idx, size_type data_idx) {
        const auto index = fmm.ab(element_idx, data_idx);
        return {fmm.ab2c_.data() + index, Size};
    }

    inline const_value_type operator()(const fmm_type& fmm,
                                       size_type element_idx,
                                       size_type data_idx) const {
        const auto index = fmm.ab(element_idx, data_idx);
        return {fmm.ab2c_.data() + index, Size};
    }
};

template <typename T, class Tag>
struct flat_multimap_item_accessor<T, 1, Tag> {
    using fmm_type =
        typename std::conditional<std::is_const<T>::value,
                                  const flat_multimap<typename std::remove_const<T>::type, 1, Tag>,
                                  flat_multimap<T, 1, Tag>>::type;
    using size_type = typename fmm_type::size_type;
    using const_value_type = typename flat_multimap_data_iterator<T, 1, Tag>::value_type;
    using value_type = typename flat_multimap_data_iterator<T, 1, Tag>::value_type&;

    inline value_type operator()(fmm_type& fmm,
                                 size_type element_idx,
                                 size_type data_idx) noexcept {
        const auto index = fmm.ab(element_idx, data_idx);
        return fmm.ab2c_[index];
    }

    inline const_value_type operator()(const fmm_type& fmm,
                                       size_type element_idx,
                                       size_type data_idx) const noexcept {
        const auto index = fmm.ab(element_idx, data_idx);
        return fmm.ab2c_[index];
    }
};

}  // namespace

/**
 * random-access data structure where elements can have 0 or more associated values.
 * Number of elements and the size of their associated data is defined in the constructor
 * and cannot be modified afterward.
 *
 * \tparam T stored data type
 * \tparam Size number of values per data
 * \tparam Tag data structure implementation
 */
template <typename T, int Size, class Tag>
class flat_multimap {
  public:
    /// index used to reference elements
    using size_type = typename traits::dynamic_vector<T, Tag>::size_type;
    /// array mapping element index a to the shift ab for the index pair (element, data)
    using a2ab_type = typename traits::dynamic_vector<size_type, Tag>::ro_type;
    /// array mapping from ab to the value c of the data b at element a
    using ab2c_type = typename traits::dynamic_vector<T, Tag>::rw_type;
    using iterator = flat_multimap_element_iterator<T, Size, Tag>;
    using const_iterator = flat_multimap_element_iterator<const T, Size, Tag>;
    using const_element_type = flat_multimap_element<const T, Size, Tag>;
    using element_type = flat_multimap_element<T, Size, Tag>;

    friend struct flat_multimap_item_accessor<T, Size, Tag>;
    friend class flat_multimap_data_iterator<T, Size, Tag>;
    friend class flat_multimap_data_iterator<const T, Size, Tag>;

    /**
     * \brief Constructs the container with value \a value
     * \param num_data contains the number of value stored for each element
     * \param value data default to 0
     */
    inline explicit flat_multimap(const a2ab_type& num_data, T value = {})
        : flat_multimap(flat_multimap::build_a2ab(num_data, value)) {}

    /**
     * \brief Construct an empty container
     */
    inline flat_multimap()
        : a2ab_(1) {}

    /**
     * \name Iterators
     * \{
     */

    /// returns an iterator to the beginning
    inline iterator begin() noexcept {
        return {*this, 0};
    }

    /// returns an iterator to the beginning
    inline const_iterator begin() const noexcept {
        return {*this, 0};
    }

    /// returns a constant iterator to the beginning
    inline const_iterator cbegin() const noexcept {
        return begin();
    }

    /// returns an iterator to the end
    inline iterator end() noexcept {
        return {*this, size()};
    }

    /// returns an iterator to the end
    inline const_iterator end() const noexcept {
        return {*this, size()};
    }

    /// returns a constant iterator to the end
    inline const_iterator cend() const noexcept {
        return end();
    }

    /** \} */

    /**
     * \name Element access
     * \{
     *
     */

    /**
     * \return data reference at the given \a element_idx and \a data_idx
     */
    inline typename flat_multimap_item_accessor<T, Size, Tag>::const_value_type operator()(
        size_type element_idx,
        size_type data_idx) const noexcept {
        return flat_multimap_item_accessor<T, Size, Tag>()(*this, element_idx, data_idx);
    }

    /**
     * \return data at the given \a element_idx and \a data_idx
     */
    inline typename flat_multimap_item_accessor<T, Size, Tag>::value_type operator()(
        size_type element_idx,
        size_type data_idx) {
        return flat_multimap_item_accessor<T, Size, Tag>()(*this, element_idx, data_idx);
    }

    /**
     * \return data collection at specified location \a element_idx.
     */
    inline const_element_type operator[](size_type element_idx) const noexcept {
        return {*this, element_idx};
    }

    /**
     * \return data collection at specified location \a element_idx.
     */
    inline element_type operator[](size_type element_idx) noexcept {
        return {*this, element_idx};
    }

    /**
     * \return data collection at specified location \a element_idx, with bounds checking.
     */
    inline const_element_type at(size_type element_idx) const {
        if (element_idx > size()) {
            throw std::out_of_range("Invalid element");
        }
        return (*this)[element_idx];
    }

    /**
     * \return data collection at specified location \a element_idx, with bounds checking.
     */
    inline element_type at(size_type element_idx) {
        if (element_idx >= size()) {
            throw std::out_of_range("Invalid index " + std::to_string(element_idx));
        }
        return (*this)[element_idx];
    }

    inline void assign(const T value) noexcept {
        std::fill(ab2c_.begin(), ab2c_.end(), value);
    }

    /** \} */

    /**
     * \name Capacity
     * \{
     */

    /**
     * \return number of map elements
     */
    inline size_type size() const noexcept {
        return static_cast<size_type>(a2ab_.size()) - 1;
    }

    /**
     * \return data collection size at specified location \a element_idx
     */
    inline size_type size(size_type element_idx) const noexcept {
        return (a2ab_[element_idx + 1] - a2ab_[element_idx]) / Size;
    }

    /**
     * \return Check whether the container is empty
     */
    inline bool empty() const noexcept {
        return size() == 0;
    }

    /**
     * \return true if there is no data at location \a element_idx
     */
    inline bool empty(size_type element_idx) const noexcept {
        return size(element_idx) == 0;
    }

    /**
     * Update organization of data
     * \param num_data contains the number of value stored for each element
     * \param keep_data reset the data buffer
     * (undefined behavior if the sum of data in \a num_data if different than the previous one)
     * \param value data default to 0
     */
    void reshape(const a2ab_type& num_data, bool keep_data = false, T value = {}) {
        const auto& construct = build_a2ab(num_data, value);
        a2ab_ = std::get<0>(construct);
        if (!keep_data) {
            ab2c_ = ab2c_type(std::get<1>(construct), value);
        }
    }

    /** \} */

  protected:
    /**
     * \return index in ab2c array of the data at the given \a element_idx and \a data_idx
     */
    inline size_type ab(size_type element_idx, size_type data_idx) const noexcept {
        const auto shift = a2ab_[element_idx];
        return shift + data_idx * Size;
    }

    inline const a2ab_type& a2ab() const noexcept {
        return a2ab_;
    }
    inline const ab2c_type& ab2c() const noexcept {
        return ab2c_;
    }

  private:
    /// mapping from element index a to the shift ab for the index pair (element, data)
    a2ab_type a2ab_;
    /// mapping from ab to the value c of the data b at element a
    ab2c_type ab2c_;

    inline explicit flat_multimap(const std::tuple<a2ab_type, size_type, T>& fields)
        : a2ab_(std::get<0>(fields))
        , ab2c_(std::get<1>(fields), std::get<2>(fields)) {}

    inline static std::tuple<a2ab_type, size_type, T> build_a2ab(const a2ab_type& num_data,
                                                                 T default_value) {
        typename traits::dynamic_vector<size_type, Tag>::rw_type a2ab(num_data.size() + 1);
        size_type ab2c_size{};
        a2ab[0] = 0;
        for (size_type element_idx = {}; element_idx < num_data.size(); ++element_idx) {
            ab2c_size += num_data[element_idx] * Size;
            a2ab[element_idx + 1] = ab2c_size;
        }
        return std::make_tuple(std::move(a2ab), ab2c_size, default_value);
    }
};

}  // namespace zee
