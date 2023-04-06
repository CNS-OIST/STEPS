/**
 * \file Provides implementation of flat_multimap data structure.
 */
#pragma once

#include <iterator>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#ifdef STEPS_USE_DIST_MESH
#include <Omega_h_array.hpp>
#include <Omega_h_for.hpp>
#endif  // STEPS_USE_DIST_MESH
#include <gsl>

namespace steps {

/// to use the Standard Template Library containers
constexpr int fmm_stl = 0b1;
#ifdef STEPS_USE_DIST_MESH
/// to use Omega_h containers
constexpr int fmm_osh = 0b10;
constexpr int fmm_default_policy = fmm_osh;
#else   // STEPS_USE_DIST_MESH
constexpr int fmm_default_policy = fmm_stl;
#endif  // STEPS_USE_DIST_MESH
/// to enable data padding in the ab2c array
constexpr int fmm_ab2c_padding = 0b100;
constexpr int fmm_backend_mask = 0b11;

namespace util {
namespace traits {

/**
 * Data trait providing types used to manipulate a dynamic random-access
 * container
 * \tparam T the type of the elements
 * \tparam Policy Data structure family (STL, OSH, ...)
 * \tparam Allocator Allocator class (only used by std::vector)
 */
template <class T, int Policy, class Allocator = std::allocator<T>>
struct dynamic_vector {
    static_assert(true, "missing implementation for this policy");
};

template <class T, class Allocator>
struct dynamic_vector<T, fmm_stl, Allocator> {
    /// random access container for read operations only.
    using ro_type = std::vector<T, Allocator>;
    /// random access container for read and write operations.
    using rw_type = std::vector<T, Allocator>;
    /// the type of the data indices.
    using size_type = std::size_t;
    /// assigns the given value to the elements of the vector.
    static inline void fill(rw_type& vector, const T& value) {
        std::fill(vector.begin(), vector.end(), value);
    }
};

#ifdef STEPS_USE_DIST_MESH

template <class T, class Allocator>
struct dynamic_vector<T, fmm_osh, Allocator> {
    /// random access container for read operations only.
    using ro_type = Omega_h::Read<T>;
    /// random access container for read and write operations.
    using rw_type = Omega_h::Write<T>;
    /// the type of the data indices.
    using size_type = Omega_h::LO;
    /// assigns the given value to the elements of the vector.
    static inline void fill(rw_type& array, const T& value) {
        Omega_h::fill(array, value);
    }
};

#endif  // STEPS_USE_DIST_MESH

}  // namespace traits

template <typename T, int Size, int Policy>
class flat_multimap;

/**
 * Iterator over the values of a specific element data.
 * \tparam T stored data type
 * \tparam Size number of values per data
 * \tparam Policy data structure implementation
 */
template <typename T, int Size, int Policy>
class flat_multimap_data_iterator {
  public:
    using iterator_category = std::input_iterator_tag;
    using value_type = gsl::span<T>;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type&;
    using fmm_type = typename std::conditional<
        std::is_const<T>::value,
        const flat_multimap<typename std::remove_const<T>::type, Size, Policy>,
        flat_multimap<T, Size, Policy>>::type;
    using size_type = typename fmm_type::size_type;

    inline flat_multimap_data_iterator(fmm_type& fmm,
                                       size_type element_idx,
                                       size_type data_idx = 0) noexcept {
        data_ = {fmm.ab2c_.data() + fmm.ab(element_idx, data_idx), Size};
    }

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

    value_type& operator*() noexcept {
        return data_;
    }

    const value_type& operator*() const noexcept {
        return data_;
    }

  private:
    value_type data_;
};

template <typename T, int Policy>
class flat_multimap_data_iterator<T, 1, Policy> {
  public:
    using iterator_category = std::input_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = T*;
    using reference = T&;
    using fmm_type = typename std::conditional<
        std::is_const<T>::value,
        const flat_multimap<typename std::remove_const<T>::type, 1, Policy>,
        flat_multimap<T, 1, Policy>>::type;
    using size_type = typename fmm_type::size_type;

    inline flat_multimap_data_iterator(fmm_type& fmm,
                                       size_type element_idx,
                                       size_type data_idx = 0) {
        data_ = fmm.ab2c_.data() + fmm.ab(element_idx, data_idx);
    }

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

    inline const value_type& operator*() const noexcept {
        return *data_;
    }

    inline value_type& operator*() noexcept {
        return *data_;
    }

  private:
    T* data_;
};

template <typename T, int Size, int Policy>
class flat_multimap_element {
  public:
    using iterator = flat_multimap_data_iterator<T, Size, Policy>;
    using const_iterator = flat_multimap_data_iterator<const T, Size, Policy>;
    using fmm_type = typename std::conditional<
        std::is_const<T>::value,
        const flat_multimap<typename std::remove_const<T>::type, Size, Policy>,
        flat_multimap<T, Size, Policy>>::type;
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
        return {fmm_, element_idx_, fmm_.size(element_idx_)};
    }
    /// returns an iterator to the end
    inline const_iterator end() const noexcept {
        return {fmm_, element_idx_, fmm_.size(element_idx_)};
    }
    /// returns a constant iterator to the end
    inline const_iterator cend() const noexcept {
        return {fmm_, element_idx_, fmm_.size(element_idx_)};
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

    /**
     * \return the element index
     */
    inline size_type element_idx() const noexcept {
        return element_idx_;
    }

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

template <typename T, int Size, int Policy>
class flat_multimap_element_iterator {
  public:
    using iterator_category = std::input_iterator_tag;
    using value_type = flat_multimap_element<T, Size, Policy>;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type&;
    using fmm_type = typename std::conditional<
        std::is_const<T>::value,
        const flat_multimap<typename std::remove_const<T>::type, Size, Policy>,
        flat_multimap<T, Size, Policy>>::type;
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

template <typename T, int Size, int Policy>
struct flat_multimap_item_accessor {
    using fmm_type = typename std::conditional<
        std::is_const<T>::value,
        const flat_multimap<typename std::remove_const<T>::type, Size, Policy>,
        flat_multimap<T, Size, Policy>>::type;
    using size_type = typename fmm_type::size_type;
    using value_type = typename flat_multimap_data_iterator<T, Size, Policy>::value_type;
    using const_value_type =
        typename flat_multimap_data_iterator<const T, Size, Policy>::value_type;

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

template <typename T, int Policy>
struct flat_multimap_item_accessor<T, 1, Policy> {
    using fmm_type = typename std::conditional<
        std::is_const<T>::value,
        const flat_multimap<typename std::remove_const<T>::type, 1, Policy>,
        flat_multimap<T, 1, Policy>>::type;
    using size_type = typename fmm_type::size_type;
    using const_value_type = typename flat_multimap_data_iterator<T, 1, Policy>::value_type;
    using value_type = typename flat_multimap_data_iterator<T, 1, Policy>::value_type&;

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


/**
 * random-access data structure where elements can have 0 or more associated
 * data. Each data is made of a fixed number of values.
 * Number of elements and the number of their associated data is defined
 * in the constructor and cannot be modified afterward.
 *
 * \tparam T stored data type
 * \tparam Size number of values per data
 * \tparam Policy bit-mask describing policies:
 *   - fmm_osh to use Omega_h array, fmm_stl to use std::vector
 *   - fmm_ab2c_padding to add padding between element data for the ab2c array
 */
template <typename T, int Size, int Policy = fmm_default_policy>
class flat_multimap {
  public:
    static constexpr int backend() {
        return Policy & fmm_backend_mask;
    }
    static_assert(backend() != 0, "backend has to be specified in Policy template parameter");
    static constexpr bool ab2c_padding() {
        return Policy & fmm_ab2c_padding;
    }
    /// index used to reference elements
    using size_type = typename traits::dynamic_vector<T, backend()>::size_type;
    /// array mapping element index a to the shift ab for the index pair (element,
    /// data)
    using a2ab_type = typename traits::dynamic_vector<size_type, backend()>::ro_type;
    /// array mapping from ab to the value c of the data b at element a
    using ab2c_type = typename traits::dynamic_vector<T, backend()>::rw_type;
    using iterator = flat_multimap_element_iterator<T, Size, Policy>;
    using const_iterator = flat_multimap_element_iterator<const T, Size, Policy>;
    using const_element_type = flat_multimap_element<const T, Size, Policy>;
    using element_type = flat_multimap_element<T, Size, Policy>;

    friend struct flat_multimap_item_accessor<T, Size, Policy>;
    friend class flat_multimap_data_iterator<T, Size, Policy>;
    friend class flat_multimap_data_iterator<const T, Size, Policy>;

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
        : a2ab_(1)
        , occupancy_rate_(0)
        , max_num_data_per_element_(0)
        , num_data_(0) {}

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
    inline typename flat_multimap_item_accessor<T, Size, Policy>::const_value_type operator()(
        size_type element_idx,
        size_type data_idx) const noexcept {
        return flat_multimap_item_accessor<T, Size, Policy>()(*this, element_idx, data_idx);
    }

    /**
     * \return data at the given \a element_idx and \a data_idx
     */
    inline typename flat_multimap_item_accessor<T, Size, Policy>::value_type operator()(
        size_type element_idx,
        size_type data_idx) {
        return flat_multimap_item_accessor<T, Size, Policy>()(*this, element_idx, data_idx);
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
     * \return data collection at specified location \a element_idx, with bounds
     * checking.
     */
    inline const_element_type at(size_type element_idx) const {
        if (element_idx > size()) {
            throw std::out_of_range("Invalid element");
        }
        return (*this)[element_idx];
    }

    /**
     * \return data collection at specified location \a element_idx, with bounds
     * checking.
     */
    inline element_type at(size_type element_idx) {
        if (element_idx >= size()) {
            throw std::out_of_range("Invalid index " + std::to_string(element_idx));
        }
        return (*this)[element_idx];
    }

    /// Assign value to all the elements in the flat_multimap
    inline void assign(const T value) noexcept {
        traits::dynamic_vector<T, backend()>::fill(ab2c_, value);
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
     * \return total number of data.
     * The sum of number of data per element
     */
    inline size_type num_data() const noexcept {
        return num_data_;
    }

    /**
     * \return total number of values.
     * The total number of data multiplied by the number of values per data.
     */
    inline size_type num_values() const noexcept {
        return num_data() * Size;
    }

    /**
     * \return data collection size at specified location \a element_idx
     */
    inline size_type size(size_type element_idx) const noexcept {
        if constexpr (ab2c_padding()) {
            return a2ab_[element_idx];
        } else {
            return (a2ab_[element_idx + 1] - a2ab_[element_idx]) / Size;
        }
    }

    /**
     * \return number of values per data
     */
    constexpr auto num_values_per_data() const noexcept -> int {
        return Size;
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
     * \return the maximum number of data among the elements
     */
    inline auto max_num_data_per_element() const noexcept -> size_type {
        return max_num_data_per_element_;
    }

    /// \return the occupancy rate of ab2c array
    inline float occupancy_rate() const noexcept {
        return occupancy_rate_;
    }

    /**
     * Update organization of data
     * \param num_data contains the number of value stored for each element
     * \param keep_data reset the data array
     * (undefined behavior if the sum of data in \a num_data if different than the
     * previous one).
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

    /**
     * \return index in ab2c array of the data at the given \a element_idx and \a
     * data_idx
     */
    inline size_type ab(size_type element_idx, size_type data_idx) const noexcept {
        if constexpr (ab2c_padding()) {
            return (max_num_data_per_element() * element_idx + data_idx) * Size;
        } else {
            const auto shift = a2ab_[element_idx];
            return shift + data_idx * Size;
        }
    }


  protected:
    inline const a2ab_type& a2ab() const noexcept {
        return a2ab_;
    }

    inline const ab2c_type& ab2c() const noexcept {
        return ab2c_;
    }

  private:
    /// mapping from element index a to the shift ab for the index pair (element,
    /// data)
    a2ab_type a2ab_;
    /// mapping from ab to the value c of the data b at element a
    ab2c_type ab2c_;
    /// occupancy rate of ab2c array
    float occupancy_rate_;
    /// maximum number of data in an element
    size_type max_num_data_per_element_;
    /// total number of data. The sum of number of data per element
    size_type num_data_;

    inline explicit flat_multimap(
        const std::tuple<a2ab_type, size_type, float, size_type, size_type, T>& fields)
        : a2ab_(std::get<0>(fields))
        , ab2c_(std::get<1>(fields), std::get<5>(fields))
        , occupancy_rate_(std::get<2>(fields))
        , max_num_data_per_element_(std::get<3>(fields))
        , num_data_(std::get<4>(fields)) {}

    /**
     * \return a tuple with the following fields:
     *   - a2ab array properly filled
     *   - ab2c size to allocate
     *   - ab2c occupancy_rate
     *   - maximum number of data per element
     *   - ab2c default filled value
     */
    inline static std::tuple<a2ab_type, size_type, float, size_type, size_type, T> build_a2ab(
        const a2ab_type& num_data,
        T default_value) {
        typename traits::dynamic_vector<size_type, backend()>::rw_type a2ab(num_data.size() + 1);
        size_type ab2c_size{};
        float occupancy_rate{1};
        size_type max_num_data_per_element{};
        size_type total_num_data{};

        if constexpr (ab2c_padding()) {
            for (size_type element_idx = {}; element_idx < num_data.size(); ++element_idx) {
                total_num_data += num_data[element_idx];
                a2ab[element_idx] = num_data[element_idx];
                max_num_data_per_element = std::max(max_num_data_per_element,
                                                    num_data[element_idx]);
            }
            a2ab[a2ab.size() - 1] = max_num_data_per_element;
            ab2c_size = num_data.size() * max_num_data_per_element * Size;
            occupancy_rate = static_cast<float>(max_num_data_per_element / ab2c_size);
        } else {
            a2ab[0] = 0;
            for (size_type element_idx = {}; element_idx < num_data.size(); ++element_idx) {
                total_num_data += num_data[element_idx];
                ab2c_size += num_data[element_idx] * Size;
                a2ab[element_idx + 1] = ab2c_size;
                max_num_data_per_element = std::max(max_num_data_per_element,
                                                    num_data[element_idx]);
            }
        }
        return std::make_tuple(std::move(a2ab),
                               ab2c_size,
                               occupancy_rate,
                               max_num_data_per_element,
                               total_num_data,
                               default_value);
    }
};

}  // namespace util
}  // namespace steps
