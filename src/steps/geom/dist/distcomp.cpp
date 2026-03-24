#include "distcomp.hpp"

#include <Omega_h_array.hpp>
#include <Omega_h_for.hpp>
#include <sstream>

#include "util/error.hpp"
#include "util/mpitools.hpp"
#include "util/vocabulary.hpp"

namespace steps::dist {

DistComp::DistComp(const mesh::compartment_name& compartment, DistMesh& mesh, double cond)
    : DistComp(compartment, mesh, compartment, cond) {}

////////////////////////////////////////////////////////////////////////////////

DistComp::DistComp(const mesh::compartment_name& compartment,
                   DistMesh& mesh,
                   mesh::compartment_physical_tag physical_tag,
                   double cond)
    : DistComp(compartment, mesh, std::to_string(physical_tag), cond) {}

////////////////////////////////////////////////////////////////////////////////

DistComp::DistComp(const mesh::compartment_name& compartment,
                   DistMesh& mesh,
                   std::string tag,
                   double cond)
    : wm::Comp(compartment, mesh, 0.0)
    , meshRef(mesh)
    , ownedVol(0.0)
    , pConductivity(cond) {
    init(mesh.getEntities(model::compartment_id(tag)));
    meshRef.addComp(model::compartment_id(compartment), this);
}

////////////////////////////////////////////////////////////////////////////////

DistComp::DistComp(const mesh::compartment_name& compartment,
                   DistMesh& mesh,
                   const std::vector<mesh::tetrahedron_global_id_t>& global_indices,
                   double cond)
    : wm::Comp(compartment, mesh, 0.0)
    , meshRef(mesh)
    , ownedVol(0.0)
    , pConductivity(cond) {
    std::vector<osh::LO> local_tets;
    for (auto& tet_global_index: global_indices) {
        const auto tet = mesh.getLocalIndex(tet_global_index, false);
        if (tet.valid()) {
            local_tets.emplace_back(tet.get());
        }
    }
    osh::Write<osh::LO> local_tets_w(local_tets.size());
    std::copy(local_tets.begin(), local_tets.end(), local_tets_w.data());
    init(local_tets_w);
    meshRef.addComp(model::compartment_id(compartment), this);
}

////////////////////////////////////////////////////////////////////////////////

DistComp::DistComp(const mesh::compartment_name& compartment,
                   DistMesh& mesh,
                   const std::vector<mesh::tetrahedron_local_id_t>& local_indices,
                   double cond)
    : wm::Comp(compartment, mesh, 0.0)
    , meshRef(mesh)
    , ownedVol(0.0)
    , pConductivity(cond) {
    osh::Write<osh::LO> local_tets_w(local_indices.size());
    osh::parallel_for(
        local_indices.size(),
        OMEGA_H_LAMBDA(osh::LO i) { local_tets_w[i] = local_indices[i].get(); });
    init(local_tets_w);
    meshRef.addComp(model::compartment_id(compartment), this);
}

////////////////////////////////////////////////////////////////////////////////

void DistComp::init(const mesh::tetrahedron_local_ids& localInds) {
    // Synchronize local indices so that all ranks have non-owned indices
    tetLocalIndices = syncLocalInds(localInds);

    ownedVol = 0;
    ownedBBoxMin.fill(std::numeric_limits<osh::Real>::max());
    ownedBBoxMax.fill(std::numeric_limits<osh::Real>::lowest());

    const auto& tets2verts = meshRef.ask_elem_verts();

    std::vector<osh::LO> ownedInds;
    ownedInds.reserve(tetLocalIndices.size());

    osh::LO cont_id = 0;
    for (const auto tet: tetLocalIndices) {
        if (meshRef.getTetComp(tet) != nullptr) {
            ArgErrLog("Tetrahedron with local index " + std::to_string(tet) +
                      " already belongs to a compartment.");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (meshRef.isOwned(tet)) {
            meshRef.setTetComp(tet, this, container::tetrahedron_id(cont_id++));

            ownedInds.emplace_back(tet.get());
            ownedVol += meshRef.getTetInfo()[tet].vol;

            const auto tet2verts = osh::gather_verts<4>(tets2verts, tet.get());
            const auto tet2x = osh::gather_vectors<4, mesh_dimensions()>(meshRef.coords(),
                                                                         tet2verts);
            for (const auto& p: tet2x) {
                for (int i = 0; i < mesh_dimensions(); ++i) {
                    ownedBBoxMin[i] = std::min(ownedBBoxMin[i], p[i]);
                    ownedBBoxMax[i] = std::max(ownedBBoxMax[i], p[i]);
                }
            }
        } else {
            meshRef.setTetComp(tet, this, {});
        }
    }

    // Set owned tetrahedron indices
    osh::Write<osh::LO> ownedInds_w(ownedInds.size());
    std::copy(ownedInds.begin(), ownedInds.end(), ownedInds_w.data());
    ownedTetLocalIndices = ownedInds_w;

    // Compute total area
    MPI_Allreduce(&ownedVol, &pVol, 1, MPI_DOUBLE, MPI_SUM, meshRef.comm_impl());
}

////////////////////////////////////////////////////////////////////////////////

mesh::tetrahedron_local_ids DistComp::syncLocalInds(
    const mesh::tetrahedron_local_ids& localInds) const {
    osh::Write<osh::LO> comp_mask(meshRef.owned_elems_mask().size(), 0);
    osh::parallel_for(
        localInds.size(), OMEGA_H_LAMBDA(osh::LO i) { comp_mask[localInds[i].get()] = 1; });
    auto sync_comp_mask = meshRef.sync_array(meshRef.dim(), osh::Read(comp_mask), 1);

    std::vector<mesh::tetrahedron_local_id_t> sync_localInds;
    sync_localInds.reserve(localInds.size());
    const auto fillInds = [&sync_comp_mask, &sync_localInds](osh::LO i) {
        if (sync_comp_mask[i] > 0) {
            sync_localInds.emplace_back(i);
        }
    };
    osh::parallel_for(sync_comp_mask.size(), fillInds);

    osh::Write<osh::LO> finalInds(sync_localInds.size());
    osh::parallel_for(
        sync_localInds.size(),
        OMEGA_H_LAMBDA(osh::LO i) { finalInds[i] = sync_localInds[i].get(); });
    return finalInds;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<mesh::tetrahedron_global_id_t> DistComp::getAllTetIndices() const {
    std::vector<mesh::tetrahedron_global_id_t> tets_global;
    tets_global.reserve(ownedTetLocalIndices.size());
    for (const auto tet: ownedTetLocalIndices) {
        tets_global.push_back(meshRef.getGlobalIndex(tet));
    }
    return meshRef.allGatherEntities(tets_global, MPI_INT64_T);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<mesh::tetrahedron_local_id_t> DistComp::getLocalTetIndices(bool owned) const {
    std::vector<mesh::tetrahedron_local_id_t> ret;
    if (owned) {
        ret.reserve(ownedTetLocalIndices.size());
        for (auto v: ownedTetLocalIndices) {
            ret.emplace_back(v);
        }
    } else {
        ret.reserve(tetLocalIndices.size());
        for (auto v: tetLocalIndices) {
            ret.emplace_back(v);
        }
    }
    return ret;
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

}  // namespace steps::dist
