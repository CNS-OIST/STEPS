/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   

 */

#include "tmpatch.hpp"

#include <unordered_set>

#include "tetmesh.hpp"
#include "tmcomp.hpp"

#include "math/point.hpp"
#include "util/collections.hpp"
#include "util/error.hpp"

namespace steps::tetmesh {

using math::point3d;

////////////////////////////////////////////////////////////////////////////////

TmPatch::TmPatch(std::string const& id,
                 Tetmesh& container,
                 std::vector<index_t> const& tris,
                 wm::Comp& wmicomp,
                 wm::Comp* wmocomp)
    : wm::Patch(id, container, wmicomp, wmocomp, 0.0)
    , pTetmesh(container) {
    // upcast the compartment pointers for this overloaded constructor

    // Note that for well-mixed compartment this will return nullptr
    auto& icomp = dynamic_cast<TmComp&>(wmicomp);
    auto ocomp = dynamic_cast<TmComp*>(wmocomp);

    // The maximum triangle index in tetrahedral mesh
    auto maxidx = pTetmesh.countTris() - 1;

    // The patch's area - contributed to from all triangles
    double area = 0.0;

    std::unordered_set<triangle_global_id> visited_tris(tris.size());
    for (auto tri_idx: tris) {
        triangle_global_id tri(tri_idx);
        if (visited_tris.count(tri) > 0) {
            continue;
        }
        visited_tris.insert(tri);

        ArgErrLogIf(tri > maxidx, "Invalid triangle index " + std::to_string(tri) + ".");

        ArgErrLogIf(pTetmesh.getTriPatch(tri) != nullptr,
                    "Triangle with index " + std::to_string(tri) + " already belongs to a patch.");

        ArgErrLogIf(pTetmesh.getTriDiffBoundary(tri) != nullptr,
                    "Triangle with index " + std::to_string(tri) +
                        " belongs to a diffusion boundary.");

        // Add triangle if compartments match those of patch, flipping
        // triangle neighbours if required.

        auto tri_tets = pTetmesh._getTriTetNeighb(tri);

        // Weiliang: bug fix for trianglar patch with pure well-mixed compartments

        // this implementation guarantee the following:
        // 1. if a tri's tet neighbor is in inner comp, its id is in tri_tet[0]
        // 2. if a tri's tet neighbor is in outer comp, its id is in tri_tet[1]

        // note that if the tet neighbor is not in either comp (well-mixed
        // compartment), there is no guarantee that tri_tets[0] != -1

        // if tri_tet[0] doesn't exist but tri_tet[1] exists
        if (tri_tets[0].unknown() && tri_tets[1].valid()) {
            auto tri_comp = pTetmesh.getTetComp(tri_tets[1]);
            // if tri_comp is inner
            // we flip it so that tri_tets[0] is now in inner comp
            if (tri_comp != nullptr && tri_comp == &icomp) {
                pTetmesh._flipTriTetNeighb(tri);
            }
        } else if (tri_tets[0].valid() && tri_tets[1].unknown()) {
            auto tri_comp = pTetmesh.getTetComp(tri_tets[0]);
            // if it is outcomp or nullptr
            // we flip it so that now tri_tets[1] is outer comp
            if (tri_comp != nullptr && tri_comp == ocomp) {
                pTetmesh._flipTriTetNeighb(tri);
            }
        }
        // both sides have tetrahedrons
        else if (tri_tets[0].valid() && tri_tets[1].valid()) {
            auto tri_comp0 = pTetmesh.getTetComp(tri_tets[0]);
            auto tri_comp1 = pTetmesh.getTetComp(tri_tets[1]);
            // if tri_comp0 is not nullptr (no well-mixed) and is outcomp,
            // we flip so that tri_tets[1] is now ocomp
            if (tri_comp0 != nullptr && tri_comp0 == ocomp) {
                pTetmesh._flipTriTetNeighb(tri);
            }
            // if tri_comp1 is not nullptr (no well-mixed) and is innercomp,
            // we flip so that tri_tets[0] is now icomp
            else if (tri_comp1 != nullptr && tri_comp1 == &icomp) {
                pTetmesh._flipTriTetNeighb(tri);
            }
            // in other cases (including both comps are well-mixed), we don't need to
            // flip

        }
        // both sides have no tet, this shouldn't happen
        else {
            ArgErrLog("Triangle with index " + std::to_string(tri) +
                      " has no neighboring tetrahedron.");
        }

        pTri_indices.push_back(tri);

        // from above we are sure that if updated_tri_tets are used in the
        // simulation, updated_tri_tets[0] will be inner tet, updated_tri_tets[1]
        // will be outer tet. note that they may not be used in the simulation, in
        // which case flipping or not doesn't matter

        auto updated_tri_tets = pTetmesh._getTriTetNeighb(tri);

        // tri_tets[0] exists
        if (updated_tri_tets[0].valid()) {
            point3d b_to_b = pTetmesh._getTriBarycenter(tri) -
                             pTetmesh._getTetBarycenter(updated_tri_tets[0]);
            if (b_to_b.dot(pTetmesh._getTriNorm(tri)) < 0) {
                pTetmesh._flipTriVerts(tri);
            }
        }
        // tri_tets[1] exists
        else if (updated_tri_tets[1].valid()) {
            const point3d b_to_b = pTetmesh._getTriBarycenter(tri) -
                                   pTetmesh._getTetBarycenter(updated_tri_tets[1]);
            if (b_to_b.dot(pTetmesh._getTriNorm(tri)) >= 0) {
                pTetmesh._flipTriVerts(tri);
            }
        }

        // Update area, patch, bounding box.
        area += pTetmesh.getTriArea(tri);
        pTetmesh.setTriPatch(tri, this);

        const auto tri_verts = pTetmesh._getTri(tri);
        for (uint j = 0; j < 3; ++j) {
            pBBox.insert(pTetmesh._getVertex(tri_verts[j]));
        }
    }  // end of loop over all tris (argument to constructor)

    pTrisN = pTri_indices.size();
    setArea(area);
    // needed for autopartition
    container._addTmPatch(this);
}

////////////////////////////////////////////////////////////////////////////////

TmPatch::~TmPatch() {
    pTetmesh._delTmPatch(this->getID());
}

////////////////////////////////////////////////////////////////////////////////

std::vector<bool> TmPatch::isTriInside(const std::vector<index_t>& tris) const {
    return util::map_membership(tris, pTri_indices);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TmPatch::getBoundMin() const {
    return util::as_vector(pBBox.min());
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TmPatch::getBoundMax() const {
    return util::as_vector(pBBox.max());
}

////////////////////////////////////////////////////////////////////////////////

void TmPatch::_addEndocyticZone(EndocyticZone* endoZone) {
    pEndoZones.push_back(endoZone);
}

}  // namespace steps::tetmesh
