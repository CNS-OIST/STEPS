#pragma once

#include <cassert>
#include <limits>
#include <stdexcept>
#include <utility>

namespace steps::util {

struct NoneType {};
static constexpr NoneType nothing{};

/**
 * The class template \a std::util::OptionalNum manages an optional contained
 * value \tparam T the type of the value to manage \tparam NoneValue the value
 * indicating there is no contained value
 */
template <typename T, T NoneValue = std::numeric_limits<T>::max()>
struct OptionalNum {
  public:
    using value_type = T;

    constexpr OptionalNum() noexcept
        : value_(none_value()){};

    constexpr OptionalNum(NoneType) noexcept
        : value_(none_value()){};

    constexpr OptionalNum(T value) noexcept {
        assert(is_value(value));
        value_ = value;
    }

    /// destroys any contained value
    constexpr OptionalNum& operator=(NoneType) noexcept {
        value_ = none_value();
        return *this;
    }

    /// assign contents
    constexpr OptionalNum& operator=(const T value) noexcept {
        assert(is_value(value));
        value_ = value;
        return *this;
    }

    /**
     * \name Observers
     * \{
     */

    constexpr T operator*() const noexcept {
        assert(has_value());
        return value_;
    }

    constexpr explicit operator bool() const noexcept {
        return has_value();
    }

    constexpr bool has_value() const noexcept {
        return value_ != none_value();
    }

    /// \return the contained value if any, throws \a std::invalid_argument
    /// otherwise
    T value() const {
        if (!has_value()) {
            throw std::invalid_argument("Unexpected none value");
        }
        return **this;
    }

    /// \return the contained value if any, \a default_value otherwise
    template <typename U>
    constexpr T value_or(U&& default_value) const noexcept {
        if (has_value()) {
            return std::move(**this);
        } else {
            return static_cast<T>(std::forward<U>(default_value));
        }
    }

    /** \} */

    /// \return the value indicating there is no contained value
    static constexpr T none_value() noexcept {
        return NoneValue;
    }

  private:
    static constexpr bool is_value(const T num) noexcept {
        return num != none_value();
    }

    T value_;
};

}  // namespace steps::util
