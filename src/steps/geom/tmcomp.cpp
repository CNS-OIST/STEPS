/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
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

// STL headers.
#include <cassert>
#include <algorithm>
#include <unordered_set>
#include <vector>
#include <sstream>

// STEPS headers.
#include "steps/common.h"
#include "steps/math/point.hpp"
#include "steps/geom/tetmesh.hpp"
#include "steps/geom/tmcomp.hpp"
#include "steps/error.hpp"
#include "steps/util/collections.hpp"
// logging
#include "easylogging++.h"
namespace stetmesh = steps::tetmesh;

using steps::math::point3d;

////////////////////////////////////////////////////////////////////////////////

stetmesh::TmComp::TmComp(std::string const & id, Tetmesh * container,
                         std::vector<uint> const & tets)
: steps::wm::Comp(id, container, 0.0)
, pTetmesh(container)
, pTetsN(0)
, pTet_indices(0)
{
    if (pTetmesh == 0)
        ArgErrLog("No mesh provided to TmComp initializer function.");

    if (tets.size() == 0)
        ArgErrLog("No tetrahedrons provided to TmComp initializer function.");

    // The maximum tetrahedron index in tetrahedral mesh
    uint maxidx = (pTetmesh-> countTets())-1;

    std::unordered_set<uint> visited_tets(tets.size());
    for (uint tet: tets) {
        if (visited_tets.count(tet)) continue;
        visited_tets.insert(tet);

        if (tet > maxidx)
            ArgErrLog("Invalid tetrahedron index "+std::to_string(tet)+".");

        if (pTetmesh->getTetComp(tet) != nullptr)
            ArgErrLog("Tetrahedron with index "+std::to_string(tet)+" already belongs to a compartment.");

        // Add the tetrahedron to this compartment.
        pTet_indices.push_back(tet);
        pTetmesh->setTetComp(tet, this);

        // Add this tetrahedron volume to the total.
        pVol += pTetmesh->getTetVol(tet);

        // Grow bounding box.
        const uint *tet_verts = pTetmesh->_getTet(tet);
        for (uint j = 0; j < 4; ++j) {
            pBBox.insert(pTetmesh->_getVertex(tet_verts[j]));
        }
    }

    pTetsN = pTet_indices.size();
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::TmComp::setVol(double vol)
{
    NotImplErrLog("""Cannot set volume of Tetmesh comp object; vol calculated internally.");
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::TmComp::getBoundMin(void) const
{
    const point3d &p = pBBox.min();
    return std::vector<double>(p.begin(),p.end());
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::TmComp::getBoundMax(void) const
{
    const point3d &p = pBBox.max();
    return std::vector<double>(p.begin(),p.end());
}

////////////////////////////////////////////////////////////////////////////////

std::vector<bool> stetmesh::TmComp::isTetInside(const std::vector<uint> &tets) const
{
    return steps::util::map_membership(tets, pTet_indices);
}

////////////////////////////////////////////////////////////////////////////////

// END

