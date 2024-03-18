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
#include "comp.hpp"
// logging
#include "util/error.hpp"

#include "util/checkpointing.hpp"

namespace steps::tetode {

////////////////////////////////////////////////////////////////////////////////

Comp::Comp(solver::Compdef* compdef)
    : pCompdef(compdef) {
    AssertLog(pCompdef != nullptr);
}

////////////////////////////////////////////////////////////////////////////////

Comp::~Comp() = default;

////////////////////////////////////////////////////////////////////////////////

void Comp::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pVol);
}

////////////////////////////////////////////////////////////////////////////////

void Comp::restore(std::fstream& cp_file) {
    util::compare(cp_file, pVol);
}

////////////////////////////////////////////////////////////////////////////////

void Comp::addTet(Tet* tet) {
    AssertLog(tet->compdef() == &def());
    tetrahedron_local_id lidx(static_cast<index_t>(pTets.size()));
    pTets.push_back(tet);
    pTets_GtoL.emplace(tet->idx(), lidx);
    pVol += tet->vol();
}

////////////////////////////////////////////////////////////////////////////////

Tet* Comp::getTet(tetrahedron_local_id lidx) {
    AssertLog(lidx < static_cast<index_t>(pTets.size()));
    return pTets[lidx.get()];
}

////////////////////////////////////////////////////////////////////////////////

steps::tetrahedron_local_id Comp::getTet_GtoL(tetrahedron_global_id gidx) {
    auto lidx_it = pTets_GtoL.find(gidx);
    AssertLog(lidx_it != pTets_GtoL.end());
    return lidx_it->second;
}

}  // namespace steps::tetode
