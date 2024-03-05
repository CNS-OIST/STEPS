#pragma once

#include <cmath>

namespace steps::util {

template <typename T>
using Vec1 = std::vector<T>;

template <typename T>
using Vec2 = std::vector<Vec1<T>>;

template <typename T>
using Vec3 = std::vector<Vec2<T>>;

template <typename T>
std::vector<T> arange(T start, T stop, T step = 1) {
    std::vector<T> values;
    for (T value = start; value < stop; value += step)
        values.push_back(value);
    return values;
}

template <typename T>
void iterMeanVec3(const Vec3<T>& data, Vec2<T>& result) {
    auto niters = data.size();
    auto ntpnts = data[0].size();
    auto nelems = data[0][0].size();
    result.assign(ntpnts, Vec1<T>(nelems, 0.0));
    for (auto t = 0u; t < ntpnts; t++) {
        for (auto e = 0u; e < nelems; e++) {
            for (auto iter = 0u; iter < niters; iter++) {
                result[t][e] += data[iter][t][e];
            }
            result[t][e] /= static_cast<T>(niters);
        }
    }
}

template <typename T>
void iterSTDVec3(const Vec3<T>& data, Vec2<T>& result) {
    auto niters = data.size();
    auto ntpnts = data[0].size();
    auto nelems = data[0][0].size();
    result.assign(ntpnts, Vec1<T>(nelems, 0.0));
    for (auto t = 0u; t < ntpnts; t++) {
        for (auto e = 0u; e < nelems; e++) {
            T sum = 0.0;
            for (auto iter = 0u; iter < niters; iter++) {
                sum += data[iter][t][e];
            }
            T mean = sum / static_cast<T>(niters);
            T standard_deviation = 0.0;
            for (auto iter = 0u; iter < niters; iter++) {
                standard_deviation += std::pow(data[iter][t][e] - mean, 2);
            }
            result[t][e] = std::sqrt(standard_deviation / static_cast<T>(niters));
        }
    }
}

}  // namespace steps::util
