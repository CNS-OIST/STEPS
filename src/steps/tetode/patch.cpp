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


// Standard library & STL headers.
#include <vector>

// STEPS headers.
#include "patch.hpp"
#include "tri.hpp"
#include "solver/compdef.hpp"
// logging
#include "util/error.hpp"
#include <easylogging++.h>
////////////////////////////////////////////////////////////////////////////////

namespace stode = steps::tetode;
namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

stode::Patch::Patch(ssolver::Patchdef * patchdef)
: pPatchdef(patchdef)
, pTris()
, pArea(0)
, pTris_GtoL()

{
    AssertLog(pPatchdef != 0);
}

////////////////////////////////////////////////////////////////////////////////

stode::Patch::~Patch()
= default;


////////////////////////////////////////////////////////////////////////////////

void stode::Patch::checkpoint(std::fstream & cp_file)
{
    cp_file.write(reinterpret_cast<char*> (&pArea), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void stode::Patch::restore(std::fstream & cp_file)
{
    cp_file.read(reinterpret_cast<char*>(&pArea), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void stode::Patch::addTri(stode::Tri * tri)
{
    AssertLog(tri->patchdef() == def());
    uint lidx = pTris.size();
    pTris.push_back(tri);
    pTris_GtoL.emplace(tri->idx(), lidx);

    pArea+=tri->area();

}

////////////////////////////////////////////////////////////////////////////////

stode::Tri * stode::Patch::getTri(uint lidx)
{
    AssertLog(lidx < pTris.size());
    return pTris[lidx];
}

////////////////////////////////////////////////////////////////////////////////

steps::triangle_id_t stode::Patch::getTri_GtoL(triangle_id_t gidx)
{
    const auto lidx_it = pTris_GtoL.find(gidx);
    AssertLog(lidx_it != pTris_GtoL.end());
    return lidx_it->second;
}
////////////////////////////////////////////////////////////////////////////////


// END
