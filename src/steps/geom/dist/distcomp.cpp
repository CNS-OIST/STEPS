#include <sstream>

#include "distcomp.hpp"
#include "easylogging++.h"
#include "util/error.hpp"
#include "util/mpitools.hpp"

namespace steps {
namespace dist {

DistComp::DistComp(const mesh::compartment_name &compartment, DistMesh &mesh,
                   double cond)
    : DistComp(compartment, mesh, compartment, cond) {
    mesh.addComp(model::compartment_id(compartment), 100, this);
}

////////////////////////////////////////////////////////////////////////////////

DistComp::DistComp(const mesh::compartment_name &compartment, DistMesh &mesh,
                   mesh::compartment_physical_tag physical_tag, double cond)
    : DistComp(compartment, mesh, std::to_string(physical_tag), cond) {
    mesh.addComp(model::compartment_id(compartment),
                 model::compartment_label(physical_tag.get()), this);
}

////////////////////////////////////////////////////////////////////////////////

DistComp::DistComp(const mesh::compartment_name &compartment, DistMesh &mesh,
                   std::string tag, double cond)
    : steps::wm::Comp(compartment, &mesh, 0.0), meshRef(mesh), ownedVol(0.0),
      pConductivity(cond) {

    for (auto tet_local_index :
         meshRef.getEntities(model::compartment_id(tag))) {
        _addTet(tet_local_index);
    }

    _computeTotalVol();
    _computeBBox();
}

////////////////////////////////////////////////////////////////////////////////

DistComp::DistComp(
    const mesh::compartment_name &compartment, DistMesh &mesh,
    const std::vector<mesh::tetrahedron_global_id_t> &global_indices,
    double cond)
    : steps::wm::Comp(compartment, &mesh, 0.0), meshRef(mesh), ownedVol(0.0),
      pConductivity(cond) {
    for (const auto &tet_global_index : global_indices) {
        const auto tet_local_index = mesh.getLocalIndex(tet_global_index);
        if (tet_local_index.valid()) {
            _addTet(tet_local_index);
        }
    }
    _computeTotalVol();
    _computeBBox();
    mesh.addComp(model::compartment_id(compartment), global_indices, this);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<mesh::tetrahedron_global_id_t> DistComp::getAllTetIndices() const {
    int local_size = ownedTetLocalIndices.size();
    std::vector<int> sizes(steps::util::mpi_comm_size(meshRef.comm_impl()));

    auto err = MPI_Allgather(&local_size, 1, MPI_INT, sizes.data(), 1, MPI_INT,
                             meshRef.comm_impl());
    if (err != MPI_SUCCESS) {
        MPI_Abort(meshRef.comm_impl(), err);
    }

    std::vector<int> offsets(sizes.size() + 1);
    std::partial_sum(sizes.begin(), sizes.end(), offsets.begin() + 1);

    std::vector<osh::GO> global_indices;
    global_indices.reserve(ownedTetLocalIndices.size());
    for (auto &ind : ownedTetLocalIndices) {
        global_indices.emplace_back(meshRef.getGlobalIndex(ind));
    }

    std::vector<osh::GO> all_indices(offsets.back());

    err = MPI_Allgatherv(global_indices.data(), global_indices.size(),
                         MPI_INT64_T, all_indices.data(), sizes.data(),
                         offsets.data(), MPI_INT64_T, meshRef.comm_impl());
    if (err != MPI_SUCCESS) {
        MPI_Abort(meshRef.comm_impl(), err);
    }

    return {all_indices.begin(), all_indices.end()};
}

////////////////////////////////////////////////////////////////////////////////

const std::vector<mesh::tetrahedron_local_id_t> &
DistComp::getLocalTetIndices(bool owned) const {
    if (owned) {
        return ownedTetLocalIndices;
    } else {
        return tetLocalIndices;
    }
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> DistComp::getBoundMin(bool local) const {
    if (local) {
        return {ownedBBoxMin.begin(), ownedBBoxMin.end()};
    } else {
        std::vector<double> minBound(mesh_dimensions());
        MPI_Allreduce(ownedBBoxMin.data(), minBound.data(), mesh_dimensions(),
                      MPI_DOUBLE, MPI_MIN, meshRef.comm_impl());
        return minBound;
    }
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> DistComp::getBoundMax(bool local) const {
    if (local) {
        return {ownedBBoxMax.begin(), ownedBBoxMax.end()};
    } else {
        std::vector<double> maxBound(mesh_dimensions());
        MPI_Allreduce(ownedBBoxMax.data(), maxBound.data(), mesh_dimensions(),
                      MPI_DOUBLE, MPI_MAX, meshRef.comm_impl());
        return maxBound;
    }
}

////////////////////////////////////////////////////////////////////////////////

void DistComp::_addTet(mesh::tetrahedron_local_id_t local_index) {
    if (meshRef.getTetComp(local_index) != nullptr) {
        ArgErrLog("Tetrahedron with local index " +
                  std::to_string(local_index) +
                  " already belongs to a compartment.");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    meshRef.setTetComp(local_index, this);
    tetLocalIndices.push_back(local_index);
    if (meshRef.isOwned(local_index)) {
        ownedVol +=
            meshRef.getTetInfo()[static_cast<size_t>(local_index.get())].vol;
        ownedTetLocalIndices.push_back(local_index.get());
    }
}

////////////////////////////////////////////////////////////////////////////////

void DistComp::_computeTotalVol() {
    MPI_Allreduce(&ownedVol, &pVol, 1, MPI_DOUBLE, MPI_SUM,
                  meshRef.comm_impl());
}

////////////////////////////////////////////////////////////////////////////////

void DistComp::_computeBBox() {
    ownedBBoxMin.fill(std::numeric_limits<osh::Real>::max());
    ownedBBoxMax.fill(std::numeric_limits<osh::Real>::min());
    const auto &tets2verts = meshRef.ask_elem_verts();
    for (const auto elem : ownedTetLocalIndices) {
        const auto tet2verts = osh::gather_verts<4>(tets2verts, elem.get());
        const auto tet2x = osh::gather_vectors<4, mesh_dimensions()>(
            meshRef.coords(), tet2verts);
        for (const auto &p : tet2x) {
            for (int i = 0; i < mesh_dimensions(); ++i) {
                ownedBBoxMin[i] = std::min(ownedBBoxMin[i], p[i]);
                ownedBBoxMax[i] = std::max(ownedBBoxMax[i], p[i]);
            }
        }
    }
}

} // namespace dist
} // namespace steps
