#include <sstream>

#include "distcomp.hpp"
#include "easylogging++.h"
#include "util/error.hpp"
#include "util/mpitools.hpp"

namespace steps::dist {

DistComp::DistComp(const mesh::compartment_name& compartment, DistMesh& mesh, double cond)
    : DistComp(compartment, mesh, compartment, cond) {
    mesh.addComp(model::compartment_id(compartment), model::compartment_label(100), this);
}

////////////////////////////////////////////////////////////////////////////////

DistComp::DistComp(const mesh::compartment_name& compartment,
                   DistMesh& mesh,
                   mesh::compartment_physical_tag physical_tag,
                   double cond)
    : DistComp(compartment, mesh, std::to_string(physical_tag), cond) {
    mesh.addComp(model::compartment_id(compartment),
                 model::compartment_label(physical_tag.get()),
                 this);
}

////////////////////////////////////////////////////////////////////////////////

DistComp::DistComp(const mesh::compartment_name& compartment,
                   DistMesh& mesh,
                   std::string tag,
                   double cond)
    : wm::Comp(compartment, &mesh, 0.0)
    , meshRef(mesh)
    , ownedVol(0.0)
    , pConductivity(cond) {
    for (auto tet_local_index: meshRef.getEntities(model::compartment_id(tag))) {
        _addTet(tet_local_index);
    }

    _computeTotalVol();
    _computeBBox();
}

////////////////////////////////////////////////////////////////////////////////

DistComp::DistComp(const mesh::compartment_name& compartment,
                   DistMesh& mesh,
                   const std::vector<mesh::tetrahedron_global_id_t>& global_indices,
                   double cond)
    : wm::Comp(compartment, &mesh, 0.0)
    , meshRef(mesh)
    , ownedVol(0.0)
    , pConductivity(cond) {
    for (const auto& tet_global_index: global_indices) {
        const auto tet_local_index = mesh.getLocalIndex(tet_global_index, false);
        if (tet_local_index.valid()) {
            _addTet(tet_local_index);
        }
    }
    _computeTotalVol();
    _computeBBox();
    mesh.addComp(model::compartment_id(compartment), global_indices, this);
}

////////////////////////////////////////////////////////////////////////////////

DistComp::DistComp(const mesh::compartment_name& compartment,
                   DistMesh& mesh,
                   const std::vector<mesh::tetrahedron_local_id_t>& local_indices,
                   double cond)
    : wm::Comp(compartment, &mesh, 0.0)
    , meshRef(mesh)
    , ownedVol(0.0)
    , pConductivity(cond) {
    for (const auto& tet_local_index: local_indices) {
        if (mesh.isOwned(tet_local_index)) {
            _addTet(tet_local_index);
        }
    }
    _computeTotalVol();
    _computeBBox();
    mesh.addComp(model::compartment_id(compartment), local_indices, this);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<mesh::tetrahedron_global_id_t> DistComp::getAllTetIndices() const {
    std::vector<mesh::tetrahedron_global_id_t> tets_global;
    tets_global.reserve(ownedTetLocalIndices.size());
    for (const auto& tet: ownedTetLocalIndices) {
        tets_global.push_back(meshRef.getGlobalIndex(tet));
    }
    return meshRef.allGatherEntities(tets_global, MPI_INT64_T);
}

////////////////////////////////////////////////////////////////////////////////

const std::vector<mesh::tetrahedron_local_id_t>& DistComp::getLocalTetIndices(bool owned) const {
    if (owned) {
        return ownedTetLocalIndices;
    } else {
        return tetLocalIndices;
    }
}

////////////////////////////////////////////////////////////////////////////////

std::vector<mesh::triangle_global_id_t> DistComp::getSurfTris() {
    auto surface_local = getSurfLocalTris();
    std::vector<mesh::triangle_global_id_t> surface_global;
    surface_global.reserve(surface_local.size());
    for (const auto& tri: surface_local) {
        surface_global.push_back(meshRef.getGlobalIndex(tri));
    }
    return meshRef.allGatherEntities(surface_global, MPI_INT64_T);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<mesh::triangle_local_id_t> DistComp::getSurfLocalTris() {
    std::unordered_set<osh::LO> entss;
    for (auto e: tetLocalIndices) {
        entss.insert(e.get());
    }
    const auto& bound_owned = meshRef.owned_bounds_mask();
    const auto& bound2elems_a2ab = meshRef.bounds2elems_a2ab(meshRef.dim() - 1, meshRef.dim());
    const auto& bound2elems_ab2b = meshRef.bounds2elems_ab2b(meshRef.dim() - 1, meshRef.dim());
    std::vector<mesh::triangle_local_id_t> surface;
    for (osh::LO boundary = 0; boundary < bound_owned.size(); boundary++) {
        if (!bound_owned[boundary]) {
            continue;
        }
        const auto num_elems = bound2elems_a2ab[boundary + 1] - bound2elems_a2ab[boundary];
        const auto el1 = bound2elems_ab2b[bound2elems_a2ab[boundary]];
        if (num_elems == 1) {
            if (entss.find(el1) != entss.end()) {
                surface.emplace_back(boundary);
            }
        } else {
            const auto el2 = bound2elems_ab2b[bound2elems_a2ab[boundary] + 1];
            if ((entss.find(el1) != entss.end()) xor (entss.find(el2) != entss.end())) {
                surface.emplace_back(boundary);
            }
        }
    }
    return surface;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> DistComp::getBoundMin(bool local) const {
    if (local) {
        return {ownedBBoxMin.begin(), ownedBBoxMin.end()};
    } else {
        std::vector<double> minBound(mesh_dimensions());
        MPI_Allreduce(ownedBBoxMin.data(),
                      minBound.data(),
                      mesh_dimensions(),
                      MPI_DOUBLE,
                      MPI_MIN,
                      meshRef.comm_impl());
        return minBound;
    }
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> DistComp::getBoundMax(bool local) const {
    if (local) {
        return {ownedBBoxMax.begin(), ownedBBoxMax.end()};
    } else {
        std::vector<double> maxBound(mesh_dimensions());
        MPI_Allreduce(ownedBBoxMax.data(),
                      maxBound.data(),
                      mesh_dimensions(),
                      MPI_DOUBLE,
                      MPI_MAX,
                      meshRef.comm_impl());
        return maxBound;
    }
}

////////////////////////////////////////////////////////////////////////////////

void DistComp::_addTet(mesh::tetrahedron_local_id_t local_index) {
    if (meshRef.getTetComp(local_index) != nullptr) {
        ArgErrLog("Tetrahedron with local index " + std::to_string(local_index) +
                  " already belongs to a compartment.");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    meshRef.setTetComp(local_index, this);
    tetLocalIndices.push_back(local_index);
    if (meshRef.isOwned(local_index)) {
        ownedVol += meshRef.getTetInfo()[static_cast<size_t>(local_index.get())].vol;
        ownedTetLocalIndices.push_back(local_index);
    }
}

////////////////////////////////////////////////////////////////////////////////

void DistComp::_computeTotalVol() {
    MPI_Allreduce(&ownedVol, &pVol, 1, MPI_DOUBLE, MPI_SUM, meshRef.comm_impl());
}

////////////////////////////////////////////////////////////////////////////////

void DistComp::_computeBBox() {
    ownedBBoxMin.fill(std::numeric_limits<osh::Real>::max());
    ownedBBoxMax.fill(std::numeric_limits<osh::Real>::lowest());
    const auto& tets2verts = meshRef.ask_elem_verts();
    for (const auto elem: ownedTetLocalIndices) {
        const auto tet2verts = osh::gather_verts<4>(tets2verts, elem.get());
        const auto tet2x = osh::gather_vectors<4, mesh_dimensions()>(meshRef.coords(), tet2verts);
        for (const auto& p: tet2x) {
            for (int i = 0; i < mesh_dimensions(); ++i) {
                ownedBBoxMin[i] = std::min(ownedBBoxMin[i], p[i]);
                ownedBBoxMax[i] = std::max(ownedBBoxMax[i], p[i]);
            }
        }
    }
}

}  // namespace steps::dist
