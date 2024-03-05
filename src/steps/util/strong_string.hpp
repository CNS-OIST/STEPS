#pragma once

#include <string>

namespace steps::util {

/**
 * \brief Quick and dirty strong type wrapper on top of std::string
 */
template <typename>
class strong_string: public std::string {
  public:
    using std::string::string;
    explicit strong_string(const std::string& s)
        : std::string(s.c_str()){};
};

}  // namespace steps::util

namespace std {
/**
 * \brief Provide std::hash implementation for primitive types
 */
template <typename T>
struct hash<steps::util::strong_string<T>> {
  public:
    constexpr size_t operator()(const steps::util::strong_string<T>& s) const {
        return std::hash<std::string>()(s);
    }
};

}  // namespace std
