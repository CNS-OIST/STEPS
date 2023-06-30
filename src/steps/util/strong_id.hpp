#pragma once

#include <algorithm>
#include <functional>
#include <optional>
#include <string>
#include <typeinfo>
#include <vector>

#include <boost/functional/hash.hpp>

#include "util/collections.hpp"

namespace steps::util {

// Generic template handles all other types
template <typename, typename, typename = void>
struct is_promotion: std::false_type {};

// Specialization recognizes enabled types
template <typename TFrom, typename TTo>
struct is_promotion<
    TFrom,
    TTo,
    typename std::enable_if<std::is_arithmetic<TFrom>::value && std::is_arithmetic<TTo>::value &&
                            std::is_integral<TFrom>::value == std::is_integral<TTo>::value &&
                            std::numeric_limits<TFrom>::is_specialized &&
                            std::numeric_limits<TTo>::is_specialized>::type>
    : std::integral_constant<bool,
                             std::numeric_limits<TFrom>::max() <= std::numeric_limits<TTo>::max()> {
};

template <typename TFrom, typename TTo>
struct is_conversion: std::false_type {};

template <>
struct is_conversion<int, signed char>: std::true_type {};
template <>
struct is_conversion<unsigned int, unsigned char>: std::true_type {};
template <>
struct is_conversion<int, short>: std::true_type {};
template <>
struct is_conversion<unsigned int, unsigned short>: std::true_type {};

template <typename T, typename Parameter, typename = std::enable_if<std::is_arithmetic<T>::value>>
class strong_id final {
    T m_value;

  public:
    using value_type = T;

    constexpr strong_id() noexcept
        : m_value(unknown_value()) {}
    constexpr strong_id(std::nullopt_t /*unused*/) noexcept
        : m_value(unknown_value()) {}

    template <typename U,
              typename = typename std::enable_if<std::is_same<T, U>::value ||
                                                 is_promotion<U, T>::value>::type>
    constexpr explicit strong_id(U const& value) noexcept
        : m_value(value) {}

    template <typename U, typename = typename std::enable_if<is_conversion<U, T>::value>::type>
    constexpr static strong_id from(U const& other) noexcept {
        return strong_id(T(other));
    }

    constexpr static auto range(value_type stop) {
        return EntityIterator<strong_id>(stop);
    }

    constexpr static auto range(strong_id stop) {
        return range(stop.get());
    }

    constexpr auto range() const noexcept {
        return range(m_value);
    }

    strong_id(strong_id const&) = default;
    strong_id(strong_id&&) noexcept = default;

    strong_id& operator=(strong_id const&) = default;
    strong_id& operator=(strong_id&&) = default;

    /**
     * \return special constant meaning that the identifier is unknown
     */
    static constexpr T unknown_value() noexcept {
        return std::numeric_limits<T>::max();
    }

    /**
     * \return true if the index is set to "no value"
     */
    [[nodiscard]] constexpr bool unknown() const noexcept {
        return m_value == unknown_value();
    }

    /**
     * \return true if the index is set to a value within the normal range
     */
    [[nodiscard]] bool valid() const noexcept {
        return !unknown();
    }

    /// \return actual identifier
    constexpr T const& get() const noexcept {
        return m_value;
    }

    constexpr void set(T value) noexcept {
        m_value = value;
    }

    constexpr T& reference() noexcept {
        return m_value;
    }

    template <typename U = T, typename = std::enable_if<!std::is_same<U, bool>::value>>
    constexpr strong_id const& operator+() const noexcept {
        return *this;
    }
    template <typename U = T, typename = std::enable_if<!std::is_same<U, bool>::value>>
    constexpr strong_id operator-() const noexcept {
        return strong_id(-m_value);
    }

    template <typename U = T, typename = std::enable_if<std::is_same<U, bool>::value>>
    constexpr bool operator!() const noexcept {
        return !m_value;
    }

    strong_id& operator++() noexcept {
        ++m_value;
        return *this;
    }
    strong_id operator++(int) noexcept {
        return strong_id(m_value++);
    }

    strong_id& operator--() noexcept {
        --m_value;
        return *this;
    }
    strong_id operator--(int) noexcept {
        return strong_id(m_value--);
    }

    template <typename U>
    strong_id& operator+=(U const& other) noexcept {
        m_value += other;
        return *this;
    }
    template <typename U>
    strong_id& operator+=(strong_id<U, Parameter> const& other) noexcept {
        m_value += other.get();
        return *this;
    }

    template <typename U>
    strong_id& operator-=(U const& other) noexcept {
        m_value -= other;
        return *this;
    }
    template <typename U>
    strong_id& operator-=(strong_id<U, Parameter> const& other) noexcept {
        m_value -= other.get();
        return *this;
    }
};

template <typename T, typename P>
constexpr bool operator==(strong_id<T, P> const& lhs, T const& rhs) noexcept {
    return lhs.get() == rhs;
}
template <typename T, typename P>
constexpr bool operator==(T const& lhs, strong_id<T, P> const& rhs) noexcept {
    return lhs == rhs.get();
}
template <typename T, typename P>
constexpr bool operator==(strong_id<T, P> const& lhs, strong_id<T, P> const& rhs) noexcept {
    return lhs.get() == rhs.get();
}

template <typename T, typename P>
constexpr bool operator!=(strong_id<T, P> const& lhs, T const& rhs) noexcept {
    return lhs.get() != rhs;
}
template <typename T, typename P>
constexpr bool operator!=(T const& lhs, strong_id<T, P> const& rhs) noexcept {
    return lhs != rhs.get();
}
template <typename T, typename P>
constexpr bool operator!=(strong_id<T, P> const& lhs, strong_id<T, P> const& rhs) noexcept {
    return lhs.get() != rhs.get();
}

template <typename T, typename P>
constexpr bool operator>(strong_id<T, P> const& lhs, T const& rhs) noexcept {
    return lhs.get() > rhs;
}
template <typename T, typename P>
constexpr bool operator>(T const& lhs, strong_id<T, P> const& rhs) noexcept {
    return lhs > rhs.get();
}
template <typename T, typename P>
constexpr bool operator>(strong_id<T, P> const& lhs, strong_id<T, P> const& rhs) noexcept {
    return lhs.get() > rhs.get();
}

template <typename T, typename P>
constexpr bool operator>=(strong_id<T, P> const& lhs, T const& rhs) noexcept {
    return lhs.get() >= rhs;
}
template <typename T, typename P>
constexpr bool operator>=(T const& lhs, strong_id<T, P> const& rhs) noexcept {
    return lhs >= rhs.get();
}
template <typename T, typename P>
constexpr bool operator>=(strong_id<T, P> const& lhs, strong_id<T, P> const& rhs) noexcept {
    return lhs.get() >= rhs.get();
}

template <typename T, typename P>
constexpr bool operator<(strong_id<T, P> const& lhs, T const& rhs) noexcept {
    return lhs.get() < rhs;
}
template <typename T, typename P>
constexpr bool operator<(T const& lhs, strong_id<T, P> const& rhs) noexcept {
    return lhs < rhs.get();
}
template <typename T, typename P>
constexpr bool operator<(strong_id<T, P> const& lhs, strong_id<T, P> const& rhs) noexcept {
    return lhs.get() < rhs.get();
}

template <typename T, typename P>
constexpr bool operator<=(strong_id<T, P> const& lhs, T const& rhs) noexcept {
    return lhs.get() <= rhs;
}
template <typename T, typename P>
constexpr bool operator<=(T const& lhs, strong_id<T, P> const& rhs) noexcept {
    return lhs <= rhs.get();
}
template <typename T, typename P>
constexpr bool operator<=(strong_id<T, P> const& lhs, strong_id<T, P> const& rhs) noexcept {
    return lhs.get() <= rhs.get();
}

/// trait
template <typename>
struct is_strong_id: public std::false_type {};

template <typename T, typename P>
struct is_strong_id<strong_id<T, P>>: public std::true_type {};

/// generic deference function
template <typename T>
inline typename std::enable_if<is_strong_id<T>::value, typename T::value_type>::type deref_strongid(
    T id) noexcept {
    return id.get();
}

template <typename T>
inline typename std::enable_if<!is_strong_id<T>::value, T>::type deref_strongid(
    const T& id) noexcept {
    return id;
}

struct deref_strong_id_func {
    template <typename T>
    inline typename std::enable_if<is_strong_id<T>::value, typename T::value_type>::type operator()(
        T id) const noexcept {
        return id.get();
    }

    template <typename T>
    inline typename std::enable_if<!is_strong_id<T>::value, T>::type operator()(
        const T& id) const noexcept {
        return id;
    }
};

}  // namespace steps::util

namespace std {

/**
 * \brief Declare that all primitive types are scalars
 */
template <typename T, typename P>
struct is_scalar<::steps::util::strong_id<T, P>>: public is_scalar<T>::type {};

//
template <typename T, typename P>
class numeric_limits<steps::util::strong_id<T, P>> {
  public:
    static constexpr T max() noexcept {
        return numeric_limits<T>::max();
    };
};


/**
 * \brief Provide std::hash implementation for primitive types
 */
template <typename T, typename P>
struct hash<steps::util::strong_id<T, P>> {
    // typedef ::util::strong_id<T, P> argument_type;
    // typedef std::size_t result_type;
    constexpr size_t operator()(const steps::util::strong_id<T, P>& p) const {
        size_t seed = 0;
        boost::hash_combine(seed, typeid(P).hash_code());
        boost::hash_combine(seed, p.get());
        return seed;
    }
};

/**
 * \brief Provide implementation for boost::hash
 */
template <typename T, typename P>
size_t hash_value(const steps::util::strong_id<T, P>& p) {
    return std::hash<steps::util::strong_id<T, P>>()(p);
}

template <typename T, typename P>
string to_string(const steps::util::strong_id<T, P>& p) {
    return to_string(p.get());
}

template <typename T, typename P>
ostream& operator<<(ostream& ostr, const steps::util::strong_id<T, P>& p) {
    return ostr << p.get();
}

}  // namespace std

template <class Iterator, class T = typename std::iterator_traits<Iterator>::value_type>
std::vector<typename T::value_type> strong_type_to_value_type(Iterator begin, Iterator end) {
    std::vector<typename T::value_type> eax;
    eax.reserve(std::distance(begin, end));
    std::transform(begin, end, std::back_inserter(eax), [](const T& e) { return e.get(); });
    return eax;
}

template <typename T, typename P>
std::vector<T> strong_type_to_value_type(const std::vector<steps::util::strong_id<T, P>>& p) {
    return strong_type_to_value_type(p.begin(), p.end());
}
