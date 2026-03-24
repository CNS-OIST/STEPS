#include "distpatch.hpp"

#include <Omega_h_for.hpp>
#include <sstream>
#include <vector>

#include "geom/dist/distcomp.hpp"
#include "geom/dist/distmesh.hpp"
#include "util/error.hpp"
#include "util/mpitools.hpp"
#include "util/vocabulary.hpp"

namespace steps::dist {

DistPatch::DistPatch(const mesh::patch_name& patch,
                     DistMesh& mesh,
                     DistComp& icomp,
                     DistComp* ocomp)
    : DistPatch(patch, mesh, patch, icomp, ocomp) {}

////////////////////////////////////////////////////////////////////////////////

DistPatch::DistPatch(const mesh::patch_name& patch,
                     DistMesh& mesh,
                     mesh::patch_physical_tag physical_tag,
                     DistComp& icomp,
                     DistComp* ocomp)
    : DistPatch(patch, mesh, std::to_string(physical_tag), icomp, ocomp) {}

////////////////////////////////////////////////////////////////////////////////

DistPatch::DistPatch(const mesh::patch_name& patch,
                     DistMesh& mesh,
                     std::string tag,
                     DistComp& icomp,
                     DistComp* ocomp)
    : wm::Patch(patch, mesh, icomp, ocomp, 0.0)
    , meshRef(mesh) {
    init(mesh.getEntities(model::patch_id(tag)));
    mesh.addPatch(model::patch_id(patch), this);
}

////////////////////////////////////////////////////////////////////////////////

DistPatch::DistPatch(const mesh::patch_name& patch,
                     DistMesh& mesh,
                     const std::vector<mesh::triangle_global_id_t>& global_indices,
                     DistComp& icomp,
                     DistComp* ocomp)
    : wm::Patch(patch, mesh, icomp, ocomp, 0.0)
    , meshRef(mesh) {
    std::vector<osh::LO> local_tris;
    for (auto& tri_global_index: global_indices) {
        const auto tri = mesh.getLocalIndex(tri_global_index, false);
        if (tri.valid()) {
            local_tris.emplace_back(tri.get());
        }
    }
    osh::Write<osh::LO> local_tris_w(local_tris.size());
    std::copy(local_tris.begin(), local_tris.end(), local_tris_w.data());
    init(local_tris_w);
    mesh.addPatch(model::patch_id(patch), this);
}

////////////////////////////////////////////////////////////////////////////////

DistPatch::DistPatch(const mesh::patch_name& patch,
                     DistMesh& mesh,
                     const std::vector<mesh::triangle_local_id_t>& local_indices,
                     DistComp& icomp,
                     DistComp* ocomp)
    : wm::Patch(patch, mesh, icomp, ocomp, 0.0)
    , meshRef(mesh) {
    osh::Write<osh::LO> local_tris_w(local_indices.size());
    osh::parallel_for(
        local_indices.size(),
        OMEGA_H_LAMBDA(osh::LO i) { local_tris_w[i] = local_indices[i].get(); });
    init(local_tris_w);
    mesh.addPatch(model::patch_id(patch), this);
}

////////////////////////////////////////////////////////////////////////////////

void DistPatch::init(const mesh::triangle_local_ids& localInds) {
    // Synchronize local indices so that all ranks have non-owned indices
    triLocalIndices = syncLocalInds(localInds);

    ownedArea = 0;
    ownedBBoxMin.fill(std::numeric_limits<osh::Real>::max());
    ownedBBoxMax.fill(std::numeric_limits<osh::Real>::lowest());

    const auto& tris2verts = meshRef.ask_verts_of(Omega_h::FACE);

    std::vector<osh::LO> ownedInds;
    ownedInds.reserve(triLocalIndices.size());

    osh::LO cont_id = 0;
    for (const auto tri: triLocalIndices) {
        if (meshRef.getTriPatch(tri) != nullptr) {
            ArgErrLog("Triangle with local index " + std::to_string(tri) +
                      " already belongs to a patch.");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (meshRef.isOwned(tri)) {
            meshRef.setTriPatch(tri, this, container::triangle_id(cont_id++));

            ownedInds.emplace_back(tri.get());
            ownedArea += meshRef.getTriInfo()[tri].area;

            // Update bbox limits
            const auto tri2verts = osh::gather_verts<3>(tris2verts, tri.get());
            const auto tri2x = osh::gather_vectors<3, mesh_dimensions()>(meshRef.coords(),
                                                                         tri2verts);
            for (const auto& p: tri2x) {
                for (int i = 0; i < mesh_dimensions(); ++i) {
                    ownedBBoxMin[i] = std::min(ownedBBoxMin[i], p[i]);
                    ownedBBoxMax[i] = std::max(ownedBBoxMax[i], p[i]);
                }
            }
        } else {
            meshRef.setTriPatch(tri, this, {});
        }
    }

    // Set owned triangle indices
    osh::Write<osh::LO> ownedInds_w(ownedInds.size());
    std::copy(ownedInds.begin(), ownedInds.end(), ownedInds_w.data());
    ownedTriLocalIndices = ownedInds_w;

    // Compute total area
    MPI_Allreduce(&ownedArea, &pArea, 1, MPI_DOUBLE, MPI_SUM, meshRef.comm_impl());
}

////////////////////////////////////////////////////////////////////////////////

mesh::triangle_local_ids DistPatch::syncLocalInds(const mesh::triangle_local_ids& localInds) const {
    osh::Write<osh::LO> patch_mask(meshRef.owned_bounds_mask().size(), 0);
    osh::parallel_for(
        localInds.size(), OMEGA_H_LAMBDA(osh::LO i) { patch_mask[localInds[i].get()] = 1; });
    auto sync_patch_mask = meshRef.sync_array(meshRef.dim() - 1, osh::Read(patch_mask), 1);

    std::vector<mesh::triangle_local_id_t> sync_localInds;
    sync_localInds.reserve(localInds.size());
    const auto fillInds = [&sync_patch_mask, &sync_localInds](osh::LO i) {
        if (sync_patch_mask[i] > 0) {
            sync_localInds.emplace_back(i);
        }
    };
    osh::parallel_for(sync_patch_mask.size(), fillInds);

    osh::Write<osh::LO> finalInds(sync_localInds.size());
    osh::parallel_for(
        sync_localInds.size(),
        OMEGA_H_LAMBDA(osh::LO i) { finalInds[i] = sync_localInds[i].get(); });
    return finalInds;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<mesh::triangle_global_id_t> DistPatch::getAllTriIndices() const {
    auto local_size = ownedTriLocalIndices.size();
    std::vector<int> sizes(util::mpi_comm_size(meshRef.comm_impl()));

    auto err =
        MPI_Allgather(&local_size, 1, MPI_INT, sizes.data(), 1, MPI_INT, meshRef.comm_impl());
    if (static_cast<int>(err) != MPI_SUCCESS) {
        MPI_Abort(meshRef.comm_impl(), static_cast<int>(err));
    }

    std::vector<int> offsets(sizes.size() + 1);
    std::partial_sum(sizes.begin(), sizes.end(), offsets.begin() + 1);

    std::vector<mesh::triangle_global_id_t> global_indices;
    global_indices.reserve(ownedTriLocalIndices.size());
    for (const auto ind: ownedTriLocalIndices) {
        global_indices.emplace_back(meshRef.getGlobalIndex(ind));
    }

    std::vector<osh::GO> all_indices(offsets.back());

    err = MPI_Allgatherv(global_indices.data(),
                         global_indices.size(),
                         MPI_INT64_T,
                         all_indices.data(),
                         sizes.data(),
                         offsets.data(),
                         MPI_INT64_T,
                         meshRef.comm_impl());
    if (static_cast<int>(err) != MPI_SUCCESS) {
        MPI_Abort(meshRef.comm_impl(), static_cast<int>(err));
    }

    return {all_indices.begin(), all_indices.end()};
}

////////////////////////////////////////////////////////////////////////////////

std::vector<mesh::triangle_local_id_t> DistPatch::getLocalTriIndices(bool owned) const {
    std::vector<mesh::triangle_local_id_t> ret;
    if (owned) {
        ret.reserve(ownedTriLocalIndices.size());
        for (auto v: ownedTriLocalIndices) {
            ret.emplace_back(v);
        }
    } else {
        ret.reserve(triLocalIndices.size());
        for (auto v: triLocalIndices) {
            ret.emplace_back(v);
        }
    }
    return ret;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> DistPatch::getBoundMin(bool local) const {
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

std::vector<double> DistPatch::getBoundMax(bool local) const {
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

}  // namespace steps::dist
