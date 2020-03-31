#pragma once

#include <numeric>
#include <stdexcept>

#include <Omega_h_adj.hpp>
#include <Omega_h_bbox.hpp>
#include <Omega_h_mesh.hpp>

#include <hadoken/format/format.hpp>

#include "common.hpp"

namespace zee {
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
static inline osh::Few<osh::LO, Dim + 1> gather_neighbours(osh::Adj const& boundary2elems,
                                                           osh::LOs const& elems2boundary,
                                                           osh::LO elem) {
    osh::Few<osh::LO, Dim + 1> neighbours_indexes(zee::initializer_list<Dim + 1>(-1));
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
mesh::element_id point_to_elem(osh::Mesh& mesh, osh::Vector<Dim> point) {
    const osh::BBox<Dim> bbox = osh::get_bounding_box<Dim>(&mesh);
    for (int i = 0; i < Dim; ++i) {
        point[i] = point[i] * (bbox.max[i] - bbox.min[i]) + bbox.min[0];
    }
    if (point < bbox.min || bbox.max < point) {
        return mesh::element_id(-1);
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
            return mesh::element_id(elem);
        }
    }
    return mesh::element_id(-1);
}

struct invalid_mesh_dim: public std::logic_error {
    explicit invalid_mesh_dim(int dim)
        : std::logic_error(hadoken::scat("Invalid mesh dimension: ", dim)) {}
    virtual ~invalid_mesh_dim() noexcept;
};
}  // namespace zee
