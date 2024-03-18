/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
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

#pragma once

#include <array>
#include <istream>
#include <map>
#include <memory>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "math/point.hpp"
#include "util/common.hpp"
#include "util/error.hpp"
#include "util/strong_id.hpp"
#include "util/strong_ra.hpp"

namespace steps::util {

// Generic checkpoint/restore methods

// Only allow trivially copyable types to avoid potential errors like we had with std::string
template <typename T,
          typename std::enable_if<std::is_trivially_copyable<T>::value, bool>::type = true>
void checkpoint(std::ostream& ostr, const T& v) {
    ostr.write(reinterpret_cast<const char*>(&v), sizeof(T));
}

template <typename T,
          typename std::enable_if<std::is_trivially_copyable<T>::value, bool>::type = true>
void restore(std::istream& istr, T& v) {
    istr.read(reinterpret_cast<char*>(&v), sizeof(T));
}

// Classes that have checkpoint and restore methods

template <typename T,
          typename std::enable_if<std::is_member_function_pointer_v<decltype(&T::checkpoint)>,
                                  bool>::type = true>
void checkpoint(std::ostream& ostr, const T& v) {
    v.checkpoint(ostr);
}

template <typename T,
          typename std::enable_if<std::is_member_function_pointer_v<decltype(&T::restore)>,
                                  bool>::type = true>
void restore(std::istream& istr, T& v) {
    v.restore(istr);
}

// Forward declarations for all overloads

void checkpoint(std::ostream& ostr, const std::string& v);
void restore(std::istream& istr, std::string& v);
template <typename Index, typename T>
void checkpoint(std::ostream& ostr, strong_id<Index, T> v);
template <typename T>
void restore(std::istream& istr, strong_id<index_t, T>& v);
template <typename T>
void checkpoint(std::ostream& ostr, const std::vector<T>& v, bool with_size = true);
template <typename T>
void restore(std::istream& istr, size_t nelems, std::vector<T>& v);
template <typename T, typename V>
void checkpoint(std::ostream& ostr, const std::pair<T, V>& v);
template <typename T, typename V>
void restore(std::istream& istr, std::pair<T, V>& v);
template <typename T, typename P>
void checkpoint(std::ostream& ostr, const util::strongid_vector<T, P>& v);
template <typename T, typename P>
void restore(std::istream& istr, strongid_vector<T, P>& v);
template <typename T>
void checkpoint(std::ostream& ostr, const std::set<T>& v);
template <typename T>
void restore(std::istream& istr, std::set<T>& v);
template <typename K, typename V>
void checkpoint(std::ostream& ostr, const std::map<K, V>& v);
template <typename K, typename V>
void restore(std::istream& istr, std::map<K, V>& v);
template <typename K, typename V>
void checkpoint(std::ostream& ostr, const std::unordered_map<K, V>& v);
template <typename K, typename V>
void restore(std::istream& istr, std::unordered_map<K, V>& v);
template <typename T>
void checkpoint(std::ostream& ostr, const std::unique_ptr<T>& ptr, size_t N);
template <typename T>
void restore(std::istream& istr, std::unique_ptr<T>& ptr, size_t N);
template <typename T>
void checkpoint(std::ostream& ostr, const T* c_array, size_t N);
template <typename T>
void restore(std::istream& istr, T* c_array, size_t N);
template <typename T, size_t N>
void checkpoint(std::ostream& ostr, const std::array<T, N>& arr);
template <typename T, size_t N>
void restore(std::istream& istr, std::array<T, N>& arr);

// std::string

inline void checkpoint(std::ostream& ostr, const std::string& v) {
    size_t sz = v.size();
    checkpoint(ostr, sz);
    ostr.write(v.data(), sizeof(char) * sz);
}

inline void restore(std::istream& istr, std::string& v) {
    size_t sz;
    restore(istr, sz);
    v.resize(sz);
    istr.read(v.data(), sizeof(char) * sz);
}

// strong_id

template <typename Index, typename T>
void checkpoint(std::ostream& ostr, strong_id<Index, T> v) {
    Index value = v.get();
    ostr.write(reinterpret_cast<char*>(&value), sizeof(Index));
}

template <typename T>
void restore(std::istream& istr, strong_id<index_t, T>& v) {
    index_t value;
    istr.read(reinterpret_cast<char*>(&value), sizeof(index_t));
    v.set(value);
}

// std::vector

template <typename T>
void checkpoint(std::ostream& ostr, const std::vector<T>& v, bool with_size) {
    if (with_size) {
        checkpoint(ostr, v.size());
    }
    if constexpr (std::is_trivially_copyable<T>::value) {
        ostr.write(reinterpret_cast<const char*>(v.data()), sizeof(T) * v.size());
    } else {
        for (const auto& e: v) {
            checkpoint(ostr, e);
        }
    }
}

template <typename T>
void restore(std::istream& istr, size_t nelems, std::vector<T>& v) {
    v.assign(nelems, {});
    if constexpr (std::is_trivially_copyable<T>::value) {
        istr.read(reinterpret_cast<char*>(v.data()), sizeof(T) * nelems);
    } else {
        for (auto& e: v) {
            restore(istr, e);
        }
    }
}

template <typename T>
void restore(std::istream& istr, std::vector<T>& v) {
    size_t nelems;
    restore(istr, nelems);
    restore(istr, nelems, v);
}

// std::pair

template <typename T, typename V>
void checkpoint(std::ostream& ostr, const std::pair<T, V>& v) {
    checkpoint(ostr, v.first);
    checkpoint(ostr, v.second);
}

template <typename T, typename V>
void restore(std::istream& istr, std::pair<T, V>& v) {
    restore(istr, v.first);
    restore(istr, v.second);
}

// std::vector<math::point3d>

inline void compare(std::istream& istr,
                    const std::vector<math::point3d>& v,
                    std::string err_msg = "") {
    std::vector<math::point3d> temp_v;
    restore(istr, temp_v);
    if (v != temp_v) {
        std::stringstream s;
        if (err_msg.empty()) {
            s << "Mismatched value from data restore comparison. Previous: ";
        } else {
            s << err_msg << " Previous: ";
        }

        for (auto const& ele: temp_v) {
            s << ele[0] << "," << ele[1] << "," << ele[2] << " ";
        }

        s << " Current: ";

        for (auto const& ele: v) {
            s << ele[0] << "," << ele[1] << "," << ele[2] << " ";
        }
        CheckpointErrLog(s.str());
    }
}

// strong_id_vector

template <typename T, typename P>
void checkpoint(std::ostream& ostr, const util::strongid_vector<T, P>& v) {
    checkpoint(ostr, v.container());
}

template <typename T, typename P>
void restore(std::istream& istr, strongid_vector<T, P>& v) {
    restore(istr, v.container());
}

// std::map

template <typename K, typename V>
void checkpoint(std::ostream& ostr, const std::map<K, V>& v) {
    checkpoint(ostr, v.size());
    for (auto const& [key, value]: v) {
        checkpoint(ostr, key);
        checkpoint(ostr, value);
    }
}

template <typename K, typename V>
void restore(std::istream& istr, std::map<K, V>& v) {
    v.clear();
    size_t nelems;
    restore(istr, nelems);
    K key;
    for (size_t e = 0; e < nelems; e++) {
        restore(istr, key);
        restore(istr, v[key]);
    }
}

// std::unordered_map

template <typename K, typename V>
void checkpoint(std::ostream& ostr, const std::unordered_map<K, V>& v) {
    checkpoint(ostr, v.size());
    for (auto const& [key, value]: v) {
        checkpoint(ostr, key);
        checkpoint(ostr, value);
    }
}

template <typename K, typename V>
void restore(std::istream& istr, std::unordered_map<K, V>& v) {
    v.clear();
    size_t nelems;
    restore(istr, nelems);
    K key;
    for (size_t e = 0; e < nelems; e++) {
        restore(istr, key);
        restore(istr, v[key]);
    }
}

// std::set

template <typename T>
void checkpoint(std::ostream& ostr, const std::set<T>& v) {
    checkpoint(ostr, v.size());
    for (auto const& d: v) {
        checkpoint(ostr, d);
    }
}

template <typename T>
void restore(std::istream& istr, std::set<T>& v) {
    v.clear();
    size_t nelems;
    restore(istr, nelems);
    for (size_t e = 0; e < nelems; e++) {
        T data;
        restore(istr, data);
        v.insert(v.end(), data);
    }
}

// Check method

template <typename T>
void compare(std::istream& istr, const T& v, std::string err_msg = "") {
    T temp_v;
    restore(istr, temp_v);
    // First check that the stream is in a good state
    if (not istr.good()) {
        CheckpointErrLog("Could not read data from checkpoint.");
    }
    if (v != temp_v) {
        std::stringstream s;
        if (err_msg.empty()) {
            s << "Mismatched value from data restore comparison. Previous: ";
        } else {
            s << err_msg << " Previous: ";
        }
        s << temp_v << " Current: " << v;
        CheckpointErrLog(s.str());
    }
}

// std::unique_ptr

template <typename T>
void checkpoint(std::ostream& ostr, const std::unique_ptr<T>& ptr, size_t N) {
    checkpoint(ostr, ptr.get(), N);
}

template <typename T>
void restore(std::istream& istr, std::unique_ptr<T>& ptr, size_t N) {
    restore(istr, ptr.get(), N);
}

// C pointer

template <typename T>
void checkpoint(std::ostream& ostr, const T* c_array, size_t N) {
    checkpoint(ostr, sizeof(T));
    checkpoint(ostr, N);
    ostr.write(reinterpret_cast<const char*>(c_array), sizeof(T) * N);
}

template <typename T>
void restore(std::istream& istr, T* c_array, size_t N) {
    compare(istr, sizeof(T));
    compare(istr, N);
    istr.read(reinterpret_cast<char*>(c_array), sizeof(T) * N);
}

// std::array

template <typename T, size_t N>
void checkpoint(std::ostream& ostr, const std::array<T, N>& arr) {
    checkpoint(ostr, sizeof(T));
    checkpoint(ostr, N);
    ostr.write(reinterpret_cast<const char*>(arr.data()), sizeof(T) * N);
}

template <typename T, size_t N>
void restore(std::istream& istr, std::array<T, N>& arr) {
    compare(istr, sizeof(T));
    compare(istr, N);
    istr.read(reinterpret_cast<char*>(arr.data()), sizeof(T) * N);
}

template <typename T, size_t N>
void compare(std::istream& istr, std::array<T, N>& arr, std::string err_msg = "") {
    std::array<T, N> temp_arr;
    restore(istr, temp_arr);
    if (arr != temp_arr) {
        std::stringstream s;
        if (err_msg.empty()) {
            s << "Mismatched value from data restore comparison. Previous: ";
        } else {
            s << err_msg << " Previous: ";
        }

        for (auto const& ele: temp_arr) {
            s << ele << " ";
        }

        s << " Current: ";

        for (auto const& ele: arr) {
            s << ele << " ";
        }
        CheckpointErrLog(s.str());
    }
}

}  // namespace steps::util
