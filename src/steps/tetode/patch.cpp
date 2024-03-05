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

// STEPS headers.
#include "patch.hpp"
#include "solver/compdef.hpp"
#include "tri.hpp"
// logging
#include "util/error.hpp"

#include "util/checkpointing.hpp"

namespace steps::tetode {


Patch::Patch(solver::Patchdef* patchdef)
    : pPatchdef(patchdef)
    , pArea(0)


{
    AssertLog(pPatchdef != nullptr);
}

////////////////////////////////////////////////////////////////////////////////

Patch::~Patch() = default;


////////////////////////////////////////////////////////////////////////////////

void Patch::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pArea);
}

////////////////////////////////////////////////////////////////////////////////

void Patch::restore(std::fstream& cp_file) {
    util::compare(cp_file, pArea);
}

////////////////////////////////////////////////////////////////////////////////

void Patch::addTri(Tri* tri) {
    AssertLog(tri->patchdef() == &def());
    triangle_local_id lidx(static_cast<index_t>(pTris.size()));
    pTris.push_back(tri);
    pTris_GtoL.emplace(tri->idx(), lidx);

    pArea += tri->area();
}

////////////////////////////////////////////////////////////////////////////////

Tri* Patch::getTri(triangle_local_id lidx) {
    AssertLog(lidx < static_cast<index_t>(pTris.size()));

    return pTris[lidx.get()];
}

////////////////////////////////////////////////////////////////////////////////

steps::triangle_local_id Patch::getTri_GtoL(triangle_global_id gidx) {
    const auto lidx_it = pTris_GtoL.find(gidx);
    AssertLog(lidx_it != pTris_GtoL.end());
    return lidx_it->second;
}

}  // namespace steps::tetode
