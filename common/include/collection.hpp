#pragma once

#include <numeric>
#include <vector>

#include <petscsystypes.h>

/**
 * \file collection.hpp
 * Provide miscellaneous collection utilities
 */

namespace zee {

/**
 * \brief Internal helper class to update quickly the sum of an array
 */
class VectorSumUpdater {
  public:
    VectorSumUpdater()
        : sum_(0.0)
        , v_() {}

    void init(const std::vector<PetscScalar>& v) {
        v_.resize(v.size());
        std::copy(v.begin(), v.end(), v_.begin());
        sum_ = std::accumulate(v.begin(), v.end(), 0.0);
    }

    inline PetscScalar sum() const noexcept {
        return sum_;
    }

    inline PetscScalar operator[](const size_t idx) const noexcept {
        return v_[idx];
    }

    void update(size_t idx, PetscScalar val) {
        sum_ += val - v_[idx];
        v_[idx] = val;
    }

    inline size_t size() const noexcept {
        return v_.size();
    }

  private:
    PetscScalar sum_;
    std::vector<PetscScalar> v_;
};

}  // namespace zee
