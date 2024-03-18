#pragma once

#include <numeric>
#include <stdexcept>
#include <type_traits>

#include <Omega_h_adj.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_bbox.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_mesh.hpp>

#if USE_PETSC
#include <petscsys.h>
#include <petscvec.h>
#endif  // USE_PETSC

#include "util/collections.hpp"
#include "util/vocabulary.hpp"

namespace steps {
/**
 * \return Number of dimensions of the meshes supported by STEPS
 */
constexpr int mesh_dimensions() {
    return 3;
}

/**
 * Compute barycenter of an element
 * \tparam Dim Mesh dimension
 * \param points coordinates of the element points
 */
template <osh::Int Dim>
osh::Vector<Dim> barycenter(const osh::Matrix<Dim, Dim + 1>& points) {
    const auto& point = std::accumulate(points.begin(),
                                        points.end(),
                                        osh::Vector<Dim>(),
                                        std::plus<osh::Vector<Dim>>());
    return point / points.size();
}

/**
 * \brief Get neighbour elements of a given element
 * \tparam Dim Dimension of the mesh
 * \param boundary2elems get triangles/tetrahedrons of an edge/triangle
 * \param elems2boundary get edges/triangles of a triangle/tetrahedron
 * \param elem a given element identifier (triangle/tetrahedron)
 * \return neighbours of an element
 */
template <osh::Int Dim>
static osh::Few<osh::LO, Dim + 1> gather_neighbours(osh::Adj const& boundary2elems,
                                                    osh::LOs const& elems2boundary,
                                                    osh::LO elem) {
    osh::Few<osh::LO, Dim + 1> neighbours_indexes(util::initializer_list<Dim + 1>(-1));
    const auto elem2boundaries = osh::gather_down<Dim + 1>(elems2boundary, elem);
    const auto& b2be = boundary2elems.a2ab;  // boundary to (boundary,element)
    const auto& be2e = boundary2elems.ab2b;  // (boundary,element) to element
    for (auto boundary = 0; boundary < elem2boundaries.size(); ++boundary) {
        const auto boundary_index = elem2boundaries[boundary];
        const auto offset = b2be[boundary_index];
        if (b2be[boundary_index + 1] - offset > 1) {
            // ef2f[offset + 0 or 1] is either elem or the index of the neighbour
            neighbours_indexes[boundary] = be2e[offset + (be2e[offset] == elem ? 1 : 0)];
        }
    }
    return neighbours_indexes;
}

template <osh::Int N>
bool operator<(osh::Vector<N> const& a, osh::Vector<N> const& b) OMEGA_H_NOEXCEPT {
    for (auto i = 0; i < N; ++i) {
        if (a[i] >= b[i]) {
            return false;
        }
    }
    return true;
}

/**
 *
 * \tparam Dim mesh dimension
 * \param mesh Omega_h mesh instance
 * \param point relative position of the point in [0, 1]
 * \return mesh element index that contains the point, -1 if not found.
 */
template <osh::Int Dim>
osh::LO point_to_elem(osh::Mesh& mesh, osh::Vector<Dim> point) {
    const osh::BBox<Dim> bbox = osh::get_bounding_box<Dim>(&mesh);
    for (int i = 0; i < Dim; ++i) {
        point[i] = point[i] * (bbox.max[i] - bbox.min[i]) + bbox.min[0];
    }
    if (point < bbox.min || bbox.max < point) {
        return -1;
    }
    const auto& elem2verts = mesh.ask_elem_verts();
    const auto& coords = mesh.coords();
    for (auto elem = 0; elem < mesh.nelems(); ++elem) {
        const auto& elem_j2verts = osh::gather_verts<Dim + 1>(elem2verts, elem);
        const auto& elem_j2x = osh::gather_vectors<Dim + 1, Dim>(coords, elem_j2verts);

        osh::Matrix<Dim, Dim> M;
        osh::Vector<Dim> b;
        for (auto ir = 0; ir < Dim; ++ir) {
            for (auto ic = 0; ic < Dim; ++ic) {
                M[ic][ir] = elem_j2x[ic + 1][ir] - elem_j2x[0][ir];
            }
            b[ir] = point[ir] - elem_j2x[0][ir];
        }
        const auto coeff = invert(M) * b;
        const auto coeff_sum =
            std::accumulate(coeff.begin(), coeff.end(), osh::Real{0}, std::plus<osh::Real>());
        const auto all_coeff_positive = std::accumulate(coeff.begin(),
                                                        coeff.end(),
                                                        true,
                                                        [](bool all_positive, osh::Real c) -> bool {
                                                            return all_positive && c >= 0;
                                                        });
        if (all_coeff_positive && coeff_sum <= 1) {
            return elem;
        }
    }
    return -1;
}

struct invalid_mesh_dim: public std::logic_error {
    explicit invalid_mesh_dim(int dim)
        : std::logic_error("Invalid mesh dimension: " + std::to_string(dim)) {}
    virtual ~invalid_mesh_dim() noexcept;
};

namespace dist::mesh {

#if USE_PETSC

/// \brief true if PETSc uses 32bits for integers, which
/// can lead to an overflow when assigned a Omega_h global element.
static constexpr const bool petsc_can_overflow = std::numeric_limits<PetscInt>::max() <
                                                 std::numeric_limits<osh::GO>::max();

/// \brief Omega_h array used to store PETSc integers
using petsc_osh_element_array =
    typename std::conditional<petsc_can_overflow, osh::LOs, osh::GOs>::type;

constexpr tetrahedron_global_id_t petsc_cast(PetscInt element) {
    return tetrahedron_global_id_t(element);
}

namespace detail {

template <bool Overflow>
struct PETScCast;

template <bool Overflow>
struct PETScCast {
    inline osh::LOs operator()(const osh::GOs& array) const {
        osh::Write<osh::LO> copy(array.size());
        auto f = OMEGA_H_LAMBDA(osh::LO i) {
            copy[i] = static_cast<osh::LO>(array[i]);
        };
        osh::parallel_for(array.size(), f, "petsc_array_cast");
        return {copy};
    }

    inline constexpr PetscInt operator()(const tetrahedron_global_id_t& element) const noexcept {
        return static_cast<PetscInt>(element.get());
    }
};

template <>
struct PETScCast<false> {
    inline osh::GOs operator()(const osh::GOs& array) const noexcept {
        return array;
    }

    inline constexpr PetscInt operator()(const tetrahedron_global_id_t& element) const noexcept {
        return element.get();
    }
};

template <typename PETSC_INT, class Enable = void>
struct PETScGetValues {
    PetscErrorCode operator()(Vec x, PETSC_INT ni, const osh::GOs& ix, PetscScalar y[]) {
        osh::Write<PetscInt> ix_cp(ni);
        auto f = OMEGA_H_LAMBDA(osh::LO i) {
            ix_cp[i] = static_cast<PETSC_INT>(ix[i]);
        };
        osh::parallel_for(ni, f, "petsc_get_values");

        return VecGetValues(x, ni, ix_cp.data(), y);
    }

    PetscErrorCode operator()(Vec x, osh::GO ix, PetscScalar& y) {
        auto ix_cp = static_cast<PETSC_INT>(ix);
        return VecGetValues(x, 1, &ix_cp, &y);
    }

    const PETSC_INT* operator()(const petsc_osh_element_array& array) const noexcept {
        return array.data();
    }

    constexpr PETSC_INT operator()(int value) const noexcept {
        return value;
    }
};

// PETSC_INT == long int, osh::LO == long int
template <typename PETSC_INT>
struct PETScGetValues<
    PETSC_INT,
    typename std::enable_if<std::is_same<typename std::add_pointer<PETSC_INT>::type,
                                         std::add_pointer<osh::GO>::type>::value>::type> {
    PetscErrorCode operator()(Vec x, PETSC_INT ni, const osh::GOs& ix, PetscScalar y[]) const {
        return VecGetValues(x, ni, ix.data(), y);
    }

    PetscErrorCode operator()(Vec x, osh::GO ix, PetscScalar& y) const {
        return VecGetValues(x, PETSC_INT(1), &ix, &y);
    }

    const PETSC_INT* operator()(const petsc_osh_element_array& array) const noexcept {
        return array.data();
    }

    constexpr PETSC_INT operator()(int value) const noexcept {
        return static_cast<PETSC_INT>(value);
    }
};

// PETSC_INT == long long int, osh::GO == long int
template <typename PETSC_INT>
struct PETScGetValues<
    PETSC_INT,
    typename std::enable_if<std::numeric_limits<PETSC_INT>::max() ==
                                std::numeric_limits<osh::GO>::max() &&
                            !std::is_same<typename std::add_pointer<PETSC_INT>::type,
                                          std::add_pointer<osh::GO>::type>::value>::type> {
    PetscErrorCode operator()(Vec x, PETSC_INT ni, const osh::GOs& ix, PetscScalar y[]) const {
        return VecGetValues(x, ni, reinterpret_cast<const PETSC_INT*>(ix.data()), y);
    }

    PetscErrorCode operator()(Vec x, osh::GO ix, PetscScalar& y) const {
        return VecGetValues(x, PETSC_INT(1), reinterpret_cast<const PETSC_INT*>(&ix), &y);
    }

    const PETSC_INT* operator()(const petsc_osh_element_array& array) const noexcept {
        return reinterpret_cast<const PETSC_INT*>(array.data());
    }

    constexpr PETSC_INT operator()(int value) const noexcept {
        return static_cast<PETSC_INT>(value);
    }
};

}  // namespace detail

/// \brief Convert a tetrahedron identifier to a PETSc
constexpr PetscInt petsc_cast(tetrahedron_global_id_t element) noexcept {
    return detail::PETScCast<petsc_can_overflow>()(element);
}

/// \brief Convert an Omega_h array of global identifiers to an
/// Omega_h array of PETSc integers.
inline petsc_osh_element_array petsc_cast(const osh::GOs& elements) {
    return detail::PETScCast<petsc_can_overflow>()(elements);
}

/// \brief Convert an integer to a PetscInt
inline constexpr PetscInt petsc_int_cast(int i) noexcept {
    return detail::PETScGetValues<PetscInt>()(i);
}

/// \brief Wrapper of PETScGetValues function
inline PetscErrorCode petsc_get_values(Vec x, PetscInt ni, const osh::GOs& ix, PetscScalar y[]) {
    return detail::PETScGetValues<PetscInt>()(x, ni, ix, y);
}

/// \brief Wrapper of PETScGetValues function with size 1
inline PetscErrorCode petsc_get_values(Vec x, osh::GO ix, PetscScalar& y) {
    return detail::PETScGetValues<PetscInt>()(x, ix, y);
}

inline const PetscInt* petsc_pointer(const petsc_osh_element_array& array) noexcept {
    return detail::PETScGetValues<PetscInt>()(array);
}

#endif  // USE_PETSC

}  // namespace dist::mesh

}  // namespace steps
