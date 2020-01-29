#ifndef STEPS_STRONG_TYPE_HPP
#define STEPS_STRONG_TYPE_HPP

#include <algorithm>
#include <functional>
#include <string>
#include <vector>

namespace steps {

template<typename TFrom, typename TTo> struct is_promotion : std::false_type {};
template<typename TFrom, typename TTo> struct is_conversion : std::false_type {};

template<> struct is_conversion<int, signed char> : std::true_type {};
template<> struct is_conversion<unsigned int, unsigned char> : std::true_type {};
template<> struct is_conversion<int, short> : std::true_type {};
template<> struct is_conversion<unsigned int, unsigned short> : std::true_type {};

template<> struct is_promotion<signed char, short> : std::true_type {};

template<> struct is_promotion<unsigned char, unsigned short> : std::true_type {};

template<> struct is_promotion<signed char, int> : std::true_type {};
template<> struct is_promotion<short, int> : std::true_type {};

template<> struct is_promotion<unsigned char, unsigned int> : std::true_type {};
template<> struct is_promotion<unsigned short, unsigned int> : std::true_type {};

template<> struct is_promotion<signed char, long> : std::true_type {};
template<> struct is_promotion<short, long> : std::true_type {};
template<> struct is_promotion<int, long> : std::true_type {};

template<> struct is_promotion<unsigned char, unsigned long> : std::true_type {};
template<> struct is_promotion<unsigned short, unsigned long> : std::true_type {};
template<> struct is_promotion<unsigned int, unsigned long> : std::true_type {};

template<> struct is_promotion<signed char, long long> : std::true_type {};
template<> struct is_promotion<short, long long> : std::true_type {};
template<> struct is_promotion<int, long long> : std::true_type {};
template<> struct is_promotion<long, long long> : std::true_type {};

template<> struct is_promotion<unsigned char, unsigned long long> : std::true_type {};
template<> struct is_promotion<unsigned short, unsigned long long> : std::true_type {};
template<> struct is_promotion<unsigned int, unsigned long long> : std::true_type {};
template<> struct is_promotion<unsigned long, unsigned long long> : std::true_type {};

template<> struct is_promotion<signed char, double> : std::true_type {};
template<> struct is_promotion<unsigned char, double> : std::true_type {};
template<> struct is_promotion<short, double> : std::true_type {};
template<> struct is_promotion<unsigned short, double> : std::true_type {};
template<> struct is_promotion<int, double> : std::true_type {};
template<> struct is_promotion<unsigned int, double> : std::true_type {};
template<> struct is_promotion<long, double> : std::true_type {};
template<> struct is_promotion<unsigned long, double> : std::true_type {};
template<> struct is_promotion<float, double> : std::true_type {};

template<> struct is_promotion<signed char, long double> : std::true_type {};
template<> struct is_promotion<unsigned char, long double> : std::true_type {};
template<> struct is_promotion<short, long double> : std::true_type {};
template<> struct is_promotion<unsigned short, long double> : std::true_type {};
template<> struct is_promotion<int, long double> : std::true_type {};
template<> struct is_promotion<unsigned int, long double> : std::true_type {};
template<> struct is_promotion<long, long double> : std::true_type {};
template<> struct is_promotion<unsigned long, long double> : std::true_type {};
template<> struct is_promotion<long long, long double> : std::true_type {};
template<> struct is_promotion<unsigned long long, long double> : std::true_type {};
template<> struct is_promotion<float, long double> : std::true_type {};
template<> struct is_promotion<double, long double> : std::true_type {};

template<typename T, typename Parameter, typename = std::enable_if< std::is_arithmetic<T>::value >>
class strong_id final {
  T m_value;

public:
  using value_type = T;

  constexpr strong_id() noexcept: m_value() {}

  template<typename U, typename = typename std::enable_if<
    std::is_same<T, U>::value || is_promotion<U, T>::value
  >::type>
  constexpr strong_id(U const& value) noexcept : m_value(value) {}

  template<typename U, typename = typename std::enable_if< is_promotion<U, T>::value >::type>
  constexpr strong_id(strong_id<U, Parameter> const& other) noexcept : m_value(other.get()) {}

  template<typename U, typename = typename std::enable_if< is_conversion<U, T>::value >::type>
  constexpr static strong_id from(U const& other) noexcept {
    return strong_id(T(other));
  }

  strong_id(strong_id const&) = default;
  strong_id(strong_id &&) = default;

  strong_id& operator=(strong_id const&) = default;
  strong_id& operator=(strong_id &&) = default;

  /// \return actual identifier
  constexpr T const& get() const noexcept { return m_value; }

  template<typename U = T, typename = std::enable_if< !std::is_same<U, bool>::value  >>
  constexpr strong_id const& operator+() const noexcept {
    return *this;
  }
  template<typename U = T, typename = std::enable_if< !std::is_same<U, bool>::value  >>
  constexpr strong_id operator-() const noexcept {
    return strong_id(-m_value);
  }

  template<typename U = T, typename = std::enable_if< std::is_same<U, bool>::value >>
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

  template<typename U>
  strong_id& operator+=(U const& other) noexcept {
    m_value += other;
    return *this;
  }
  template<typename U>
  strong_id& operator+=(strong_id<U, Parameter> const& other) noexcept {
    m_value += other.get();
    return *this;
  }

  template<typename U>
  strong_id& operator-=(U const& other) noexcept {
    m_value -= other;
    return *this;
  }
  template<typename U>
  strong_id& operator-=(strong_id<U, Parameter> const& other) noexcept {
    m_value -= other.get();
    return *this;
  }
};

template<typename T, typename P>
constexpr bool operator==(strong_id<T, P> const& lhs, T const& rhs) noexcept {
  return lhs.get() == rhs;
}
template<typename T, typename P>
constexpr bool operator==(T const& lhs, strong_id<T, P> const& rhs) noexcept {
  return lhs == rhs.get();
}
template<typename T, typename P>
constexpr bool operator==(strong_id<T, P> const& lhs, strong_id<T, P> const& rhs) noexcept {
  return lhs.get() == rhs.get();
}

template<typename T, typename P>
constexpr bool operator!=(strong_id<T, P> const& lhs, T const& rhs) noexcept {
  return lhs.get() != rhs;
}
template<typename T, typename P>
constexpr bool operator!=(T const& lhs, strong_id<T, P> const& rhs) noexcept {
  return lhs != rhs.get();
}
template<typename T, typename P>
constexpr bool operator!=(strong_id<T, P> const& lhs, strong_id<T, P> const& rhs) noexcept {
  return lhs.get() != rhs.get();
}

template<typename T, typename P>
constexpr bool operator>(strong_id<T, P> const& lhs, T const& rhs) noexcept {
  return lhs.get() > rhs;
}
template<typename T, typename P>
constexpr bool operator>(T const& lhs, strong_id<T, P> const& rhs) noexcept {
  return lhs > rhs.get();
}
template<typename T, typename P>
constexpr bool operator>(strong_id<T, P> const& lhs, strong_id<T, P> const& rhs) noexcept {
  return lhs.get() > rhs.get();
}

template<typename T, typename P>
constexpr bool operator>=(strong_id<T, P> const& lhs, T const& rhs) noexcept {
  return lhs.get() >= rhs;
}
template<typename T, typename P>
constexpr bool operator>=(T const& lhs, strong_id<T, P> const& rhs) noexcept {
  return lhs >= rhs.get();
}
template<typename T, typename P>
constexpr bool operator>=(strong_id<T, P> const& lhs, strong_id<T, P> const& rhs) noexcept {
  return lhs.get() >= rhs.get();
}

template<typename T, typename P>
constexpr bool operator<(strong_id<T, P> const& lhs, T const& rhs) noexcept {
  return lhs.get() < rhs;
}
template<typename T, typename P>
constexpr bool operator<(T const& lhs, strong_id<T, P> const& rhs) noexcept {
  return lhs < rhs.get();
}
template<typename T, typename P>
constexpr bool operator<(strong_id<T, P> const& lhs, strong_id<T, P> const& rhs) noexcept {
  return lhs.get() < rhs.get();
}

template<typename T, typename P>
constexpr bool operator<=(strong_id<T, P> const& lhs, T const& rhs) noexcept {
  return lhs.get() <= rhs;
}
template<typename T, typename P>
constexpr bool operator<=(T const& lhs, strong_id<T, P> const& rhs) noexcept {
  return lhs <= rhs.get();
}
template<typename T, typename P>
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
inline
typename std::enable_if<is_strong_id<T>::value, typename T::value_type>::type
deref_strongid(T id) noexcept { return id.get(); }

template <typename T>
inline
typename std::enable_if<!is_strong_id<T>::value, T>::type
deref_strongid(const T& id) noexcept { return id; }

struct deref_strong_id_func {
  template <typename T>
  inline
  typename std::enable_if<is_strong_id<T>::value, typename T::value_type>::type
  operator()(T id) const noexcept { return id.get(); }

  template <typename T>
  inline
  typename std::enable_if<!is_strong_id<T>::value, T>::type
  operator()(const T& id) const noexcept { return id; }

};

} // namespace steps

namespace std {

/**
 * \brief Declare that all primitive types are scalars
 */
template<typename T, typename P>
struct is_scalar<::steps::strong_id<T, P>> : public is_scalar<T>::type {};

/**
 * \brief Provide std::hash implementation for primitive types
 */
template <typename T, typename P>
struct hash<::steps::strong_id<T, P>> {
public:
  constexpr size_t operator()(const ::steps::strong_id<T, P>& p) const {
    return std::hash<T>()(p.get());
  }
};

template <typename T, typename P>
string to_string(const ::steps::strong_id<T, P>& p) {
  return to_string(p.get());
}

template <typename T, typename P>
ostream& operator<<(ostream& ostr, const ::steps::strong_id<T, P>& p) {
  return ostr << p.get();
}

} // namespace std

template <class Iterator, class T = typename std::iterator_traits<Iterator>::value_type>
std::vector<typename T::value_type> strong_type_to_value_type(Iterator begin, Iterator end) {
  std::vector<typename T::value_type> eax;
  eax.reserve(std::distance(begin, end));
  std::transform(
    begin, end, std::back_inserter(eax),
    [](const T& e) { return e.get(); }
  );
  return eax;
}

template <typename T, typename P>
std::vector<T> strong_type_to_value_type(const std::vector<::steps::strong_id<T, P>>& p) {
  return strong_type_to_value_type(p.begin(), p.end());
}

#endif //!STEPS_STRONG_TYPE_HPP
