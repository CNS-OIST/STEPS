#include "distpatch.hpp"

#include <sstream>

#include <easylogging++.h>

#include "util/error.hpp"
#include "util/mpitools.hpp"

namespace steps {
namespace dist {

DistPatch::DistPatch(const mesh::patch_name &patch, DistMesh &mesh,
                     DistComp *icomp, DistComp *ocomp)
    : DistPatch(patch, mesh, patch, icomp, ocomp) {
}

////////////////////////////////////////////////////////////////////////////////

DistPatch::DistPatch(const mesh::patch_name &patch, DistMesh &mesh,
                     mesh::patch_physical_tag physical_tag, DistComp *icomp,
                     DistComp *ocomp)
    : DistPatch(patch, mesh, std::to_string(physical_tag), icomp, ocomp) {
}

////////////////////////////////////////////////////////////////////////////////

DistPatch::DistPatch(const mesh::patch_name &patch, DistMesh &mesh,
                     std::string tag, DistComp *icomp, DistComp *ocomp)
    : steps::wm::Patch(patch, &mesh, icomp, ocomp, 0.0), meshRef(mesh) {

    for (auto tri_local_index : meshRef.getEntities(model::patch_id(tag))) {
        _addTri(tri_local_index);
    }

    _computeTotalArea();
    _computeBBox();
    mesh.addPatch(model::patch_id(patch), this);
}

////////////////////////////////////////////////////////////////////////////////

DistPatch::DistPatch(const mesh::patch_name& patch,
                     DistMesh& mesh,
                     const std::vector<mesh::triangle_global_id_t>& global_indices,
                     DistComp* icomp,
                     DistComp* ocomp)
    : wm::Patch(patch, &mesh, icomp, ocomp, 0.0)
    , meshRef(mesh)
    , ownedArea(0.0) {
    for (auto& tri_global_index: global_indices) {
        const auto tri_local_index = mesh.getLocalIndex(tri_global_index, false);
        if (tri_local_index.valid()) {
            _addTri(tri_local_index);
        }
    }
    _computeTotalArea();
    _computeBBox();
    mesh.addPatch(model::patch_id(patch), global_indices, this);
}

////////////////////////////////////////////////////////////////////////////////

DistPatch::DistPatch(const mesh::patch_name& patch,
                     DistMesh& mesh,
                     const std::vector<mesh::triangle_local_id_t>& local_indices,
                     DistComp* icomp,
                     DistComp* ocomp)
    : wm::Patch(patch, &mesh, icomp, ocomp, 0.0)
    , meshRef(mesh)
    , ownedArea(0.0) {
    for (auto& tri_local_index: local_indices) {
        _addTri(tri_local_index);
    }
    _computeTotalArea();
    _computeBBox();
    mesh.addPatch(model::patch_id(patch), local_indices, this);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<mesh::triangle_global_id_t> DistPatch::getAllTriIndices() const {
    int local_size = ownedTriLocalIndices.size();
    std::vector<int> sizes(steps::util::mpi_comm_size(meshRef.comm_impl()));

    auto err = MPI_Allgather(&local_size, 1, MPI_INT, sizes.data(), 1, MPI_INT,
                             meshRef.comm_impl());
    if (err != MPI_SUCCESS) {
        MPI_Abort(meshRef.comm_impl(), err);
    }

    std::vector<int> offsets(sizes.size() + 1);
    std::partial_sum(sizes.begin(), sizes.end(), offsets.begin() + 1);

    std::vector<osh::GO> global_indices;
    global_indices.reserve(ownedTriLocalIndices.size());
    for (auto &ind : ownedTriLocalIndices) {
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

const std::vector<mesh::triangle_local_id_t> &
DistPatch::getLocalTriIndices(bool owned) const {
    if (owned) {
        return ownedTriLocalIndices;
    } else {
        return triLocalIndices;
    }
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> DistPatch::getBoundMin(bool local) const {
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

std::vector<double> DistPatch::getBoundMax(bool local) const {
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

void DistPatch::_addTri(mesh::triangle_local_id_t local_index) {
  if (meshRef.getTriPatch(local_index) != nullptr) {
    ArgErrLog("Triangle with local index " + std::to_string(local_index) +
              " already belongs to a patch.");
  }
  meshRef.setTriPatch(local_index, this);
  triLocalIndices.push_back(local_index);
  if (meshRef.isOwned(local_index)) {
    ownedArea += meshRef.getTriInfo()[static_cast<size_t>(local_index.get())].area;
    ownedTriLocalIndices.push_back(local_index);
  }
}

////////////////////////////////////////////////////////////////////////////////

void DistPatch::_computeTotalArea() {
  MPI_Allreduce(&ownedArea, &pArea, 1, MPI_DOUBLE, MPI_SUM, meshRef.comm_impl());
}

////////////////////////////////////////////////////////////////////////////////

void DistPatch::_computeBBox() {
    ownedBBoxMin.fill(std::numeric_limits<osh::Real>::max());
    ownedBBoxMax.fill(std::numeric_limits<osh::Real>::lowest());
    const auto& tris2verts = meshRef.ask_verts_of(Omega_h::FACE);
    for (const auto tri: ownedTriLocalIndices) {
        const auto tri2verts = osh::gather_verts<3>(tris2verts, tri.get());
        const auto tri2x = osh::gather_vectors<3, mesh_dimensions()>(
            meshRef.coords(), tri2verts);
        for (const auto &p : tri2x) {
            for (int i = 0; i < mesh_dimensions(); ++i) {
                ownedBBoxMin[i] = std::min(ownedBBoxMin[i], p[i]);
                ownedBBoxMax[i] = std::max(ownedBBoxMax[i], p[i]);
            }
        }
    }
}

}  // namespace dist
}  // namespace steps
