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

#include "tmcomp.hpp"

#include <sstream>
#include <unordered_set>
#include <vector>

#include "tetmesh.hpp"

#include "util/error.hpp"

namespace steps::tetmesh {

using math::point3d;

TmComp::TmComp(std::string const& id, Tetmesh& container, std::vector<index_t> const& tets)
    : wm::Comp(id, container, 0.0)
    , pTetmesh(container) {
    ArgErrLogIf(tets.empty(), "No tetrahedrons provided to TmComp initializer function.");

    // The maximum tetrahedron index in tetrahedral mesh
    auto maxidx = static_cast<index_t>(pTetmesh.countTets() - 1);

    std::unordered_set<tetrahedron_global_id> visited_tets(tets.size());
    pTet_indices.reserve(tets.size());
    for (const auto tet_idx: tets) {
        tetrahedron_global_id tet(tet_idx);
        if (visited_tets.count(tet) != 0) {
            continue;
        }
        visited_tets.insert(tet);

        ArgErrLogIf(tet > maxidx, "Invalid tetrahedron index " + std::to_string(tet) + ".");

        ArgErrLogIf(pTetmesh.getTetComp(tet) != nullptr,
                    "Tetrahedron with index " + std::to_string(tet) +
                        " already belongs to a compartment.");

        // Add the tetrahedron to this compartment.
        pTet_indices.push_back(tet);
        pTetmesh.setTetComp(tet, this);

        // Add this tetrahedron volume to the total.
        pVol += pTetmesh.getTetVol(tet);

        // Grow bounding box.
        const auto tet_verts = pTetmesh._getTet(tet);
        for (uint j = 0; j < 4; ++j) {
            pBBox.insert(pTetmesh._getVertex(tet_verts[j]));
        }
    }

    pTetsN = pTet_indices.size();

    // needed for autopartition
    pTetmesh._addTmComp(this);
}

////////////////////////////////////////////////////////////////////////////////

TmComp::~TmComp() {
    pTetmesh._delTmComp(this->getID());
}

////////////////////////////////////////////////////////////////////////////////

void TmComp::setVol(double /* vol */) {
    NotImplErrLog("Cannot set volume of Tetmesh comp object; vol calculated internally.");
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TmComp::getBoundMin() const {
    const point3d& p = pBBox.min();
    return {p.begin(), p.end()};
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TmComp::getBoundMax() const {
    const point3d& p = pBBox.max();
    return {p.begin(), p.end()};
}

////////////////////////////////////////////////////////////////////////////////

std::vector<bool> TmComp::isTetInside(const std::vector<index_t>& tets) const {
    return util::map_membership(tets, pTet_indices);
}

}  // namespace steps::tetmesh
