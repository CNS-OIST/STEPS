#pragma once

#include <array>
#include <iosfwd>
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "flat_multimap.hpp"

#if USE_PETSC
#include <petscmat.h>
#include <petscvec.h>
#endif // USE_PETSC

#if STEPS_USE_DIST_MESH
#include <Omega_h_array.hpp>
#include <Omega_h_matrix.hpp>
#include <Omega_h_vector.hpp>

#include "strong_ids.hpp"
#endif // STEPS_USE_DIST_MESH

namespace steps {
namespace util {

void wait_for_gdb();

} // namespace util
} // namespace steps

#if USE_PETSC
/// It must be called after assemblyEnd and it could present problems if run with multiple nodes
std::ostream& operator<<(std::ostream& os, const Vec& v);
/// It must be called after assemblyEnd and it could present problems if run with multiple nodes
std::ostream& operator<<(std::ostream& os, const Mat& v);
#endif // USE_PETSC


#if STEPS_USE_DIST_MESH

namespace Omega_h {

/// Pretty print of Read
template <typename T>
std::ostream& operator<<(std::ostream& os, const Read<T>& v) {
    os << '(' << v.size() << "): [";
    for (const auto i : v) {
        os << i << ", ";
    }
    return os << ']';
}
/// Pretty print of Write
template <typename T>
std::ostream& operator<<(std::ostream& os, const Write<T>& v) {
    return os << Read(v);
}
/// Pretty print of Vector
template <Int N>
std::ostream& operator<<(std::ostream& os, const Vector<N>& v) {
    os << '(' << N << "): \n[";
    for (Int i = 0; i < N; ++i) {
        os << v[i] << ", ";
    }
    return os << ']';
}
/// Pretty print of Matrix
template <Int M, Int N>
std::ostream& operator<<(std::ostream& os, const Matrix<M, N>& mat) {
    os << '(' << M << "x" << N << "): \n[";
    for (Int i = 0; i < M; ++i) {
        os << '[';
        for (Int j = 0; j < N; ++j) {
            os << mat[j][i] << ", ";
        }
        os << "],\n";
    }
    return os << ']';
}

}  // namespace Omega_h

#endif // STEPS_USE_DIST_MESH

namespace std {

/// Pretty print of std::array
template <typename T, size_t n>
ostream& operator<<(ostream& os, const array<T, n>& v) {
    os << '(' << n << "): [";
    for (const auto& i: v) {
        os << i << ", ";
    }
    return os << "]";
}

/// Pretty print of std::vector
template <typename T>
ostream& operator<<(ostream& os, const vector<T>& v) {
    os << '(' << v.size() << "): [";
    for (const auto& i: v) {
        os << i << ", ";
    }
    return os << ']';
}

/// Pretty print of std::map
template <typename K, typename V>
ostream& operator<<(ostream& os, const map<K, V>& m) {
    os << '(' << m.size() << "): {";
    for (const auto& i: m) {
        os << i.first << ": " << i.second << ", ";
    }
    return os << '}';
}

/// Pretty print of std::multimap
template <typename K, typename V>
ostream& operator<<(ostream& os, const multimap<K, V>& m) {
    os << '(' << m.size() << "): {";
    for (const auto& i: m) {
        os << i.first << ": " << i.second << ", ";
    }
    return os << '}';
}

/// Pretty print of std::unordered_map
template <typename K, typename V>
ostream& operator<<(ostream& os, const unordered_map<K, V>& m) {
    os << '(' << m.size() << "): {";
    for (const auto& i: m) {
        os << i.first << ": " << i.second << ", ";
    }
    return os << '}';
}

/// Pretty print of std::set
template <typename T>
ostream& operator<<(ostream& os, const set<T>& s) {
    os << '(' << s.size() << "): {";
    for (const auto& i: s) {
        os << i << ", ";
    }
    return os << '}';
}

/// Pretty print of std::unordered_set
template <typename T>
ostream& operator<<(ostream& os, const unordered_set<T>& s) {
    os << '(' << s.size() << "): {";
    for (const auto& i: s) {
        os << i << ", ";
    }
    return os << '}';
}

/// Pretty print of gsl::span
template <typename T>
ostream& operator<<(ostream& os, const gsl::span<T>& v) {
    os << '(' << v.size() << "): [";
    for (const auto& i: v) {
        os << i << ", ";
    }
    return os << ']';
}

/// Pretty print of flat_multimap
template <typename T, int Size, int Policy>
ostream& operator<<(ostream& os, const steps::util::flat_multimap<T, Size, Policy>& fm) {
    os << "a2ab: (" << fm.size() << "): [";
    for (const auto& i: fm) {
        os << i.size() << ", ";
    }
    os << "] ab2c: (" << fm.num_values() << "): [";
    for (const auto& i: fm) {
        for (const auto& j: i) {
            os << j << ", ";
        }
    }
    return os << ']';
}

#ifdef STEPS_USE_DIST_MESH
/// Pretty print of strong_ids
template <typename T>
ostream& operator<<(ostream& os, const steps::util::strong_ids<T>& v) {
    os << '(' << v.size() << "): [";
    for (const auto& i: v) {
        os << i << ", ";
    }
    return os << "]";
}
#endif // STEPS_USE_DIST_MESH

}  // namespace std

