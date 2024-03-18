#include "distmemb.hpp"

namespace steps::dist {

DistMemb::DistMemb(const model::membrane_id& membrane,
                   DistMesh& mesh,
                   const std::set<model::patch_id>& patches,
                   double capac)
    : pID(membrane)
    , meshRef(mesh)
    , pPatches(patches)
    , pCapacitance(capac) {
    mesh.addMembrane(model::membrane_id(membrane), this);
}

}  // namespace steps::dist
