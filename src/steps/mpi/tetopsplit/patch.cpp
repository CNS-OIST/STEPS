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

#include "patch.hpp"

// Standard library & STL headers.
#include <vector>

// STEPS headers.
#include "model/reac.hpp"
#include "solver/compdef.hpp"
#include "util/checkpointing.hpp"
#include "util/error.hpp"

#include <easylogging++.h>

namespace steps::mpi::tetopsplit {

Patch::Patch(solver::Patchdef* patchdef)
    : pPatchdef(patchdef) {
    AssertLog(pPatchdef != nullptr);
}

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
    AssertLog(tri->patchdef() == def());
    pTris.push_back(tri);
    pArea += tri->area();
}

////////////////////////////////////////////////////////////////////////////////

Tri* Patch::pickTriByArea(double rand01) const {
    if (countTris() == 0) {
        return nullptr;
    }
    if (countTris() == 1) {
        return pTris[0];
    }
    double accum = 0.0;
    double selector = rand01 * area();
    for (auto const& t: pTris) {
        accum += t->area();
        if (selector <= accum) {
            return t;
        }
    }

    return tris().back();
}

}  // namespace steps::mpi::tetopsplit
