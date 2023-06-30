#pragma once

#include <array>
#include <type_traits>

#include <petscsystypes.h>

namespace steps::util::petsc {

/**
 * This utility class allows passing from a random-access container of double to
 * a pointer of PetscScalar, given that PetscScalar can be either
 * a double or a std::complex<double> depending of whether PETSc has been compiled with or without
 * the configure option "--with-scalar-type=complex".
 * \tparam Array can be a std::array, std::vector, Omega_h::read or Omega_h::write
 * \tparam needs_copy whether a copy of the input container is required or not. It is
 * an internal template parameter that should not be set explicitly
 */

template <typename Array,
          bool needs_copy = !std::is_same_v<typename Array::value_type, PetscScalar>>
class scalars {
  public:
    scalars(const Array& array);
    const PetscScalar* data() const;
};

template <typename Array>
class scalars<Array, false> {
  public:
    inline scalars(const Array& array)
        : ref_(array) {}
    inline const PetscScalar* data() const noexcept {
        return ref_.data();
    }

  private:
    const Array& ref_;
};

template <typename T, size_t N>
class scalars<std::array<T, N>, true> {
  public:
    inline scalars(const std::array<T, N>& array) {
        std::copy(array.begin(), array.end(), copy_.begin());
    }
    inline const PetscScalar* data() const noexcept {
        return copy_.data();
    }

  private:
    std::array<PetscScalar, N> copy_;
};

template <typename Array>
class scalars<Array, true> {
  public:
    inline scalars(const Array& array) {
        copy_.reserve(array.size());
        std::copy(array.begin(), array.end(), std::back_inserter(copy_));
    }
    inline const PetscScalar* data() const noexcept {
        return copy_.data();
    }

  private:
    std::vector<PetscScalar> copy_;
};

[[maybe_unused]] constexpr PetscReal to_real(const PetscComplex& complex) {
    return complex.real();
}

[[maybe_unused]] constexpr PetscReal to_real(const PetscReal& real) {
    return real;
}

#if PETSC_USE_COMPLEX
#define MPI_PETSC_SCALAR MPI_DOUBLE_COMPLEX
#else
#define MPI_PETSC_SCALAR MPI_DOUBLE
#endif

}  // namespace steps::util::petsc
