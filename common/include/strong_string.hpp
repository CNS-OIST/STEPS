#pragma once

#include <string>

namespace zee {

/**
 * \brief Quick and dirty strong type wrapper on top of std::string
 */
template <typename>
class strong_string: public std::string {
  public:
    using std::string::string;
};

}  // namespace zee

namespace std {
/**
 * \brief Provide std::hash implementation for primitive types
 */
template <typename T>
struct hash<::zee::strong_string<T>> {
  public:
    constexpr size_t operator()(const ::zee::strong_string<T>& s) const {
        return std::hash<std::string>()(s);
    }
};

}  // namespace std
