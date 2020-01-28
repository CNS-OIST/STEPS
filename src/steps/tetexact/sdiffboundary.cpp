/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2020 Okinawa Institute of Science and Technology, Japan.
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
 
// Standard library & STL headers.
#include <vector>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/solver/sdiffboundarydef.hpp"
#include "steps/tetexact/patch.hpp"
#include "steps/tetexact/sdiffboundary.hpp"


// logging
#include "easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

namespace stex = steps::tetexact;
namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

stex::SDiffBoundary::SDiffBoundary(steps::solver::SDiffBoundarydef * sdbdef)
: pSDiffBoundarydef(sdbdef)
, pSetPatches(false)
, pPatchA(nullptr)
, pPatchB(nullptr)
, pTris()
, pTriDirection()
{
    AssertLog(sdbdef != 0);
}

////////////////////////////////////////////////////////////////////////////////

stex::SDiffBoundary::~SDiffBoundary()
= default;

////////////////////////////////////////////////////////////////////////////////

void stex::SDiffBoundary::checkpoint(std::fstream & /*cp_file*/)
{
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void stex::SDiffBoundary::restore(std::fstream & /*cp_file*/)
{
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void stex::SDiffBoundary::setPatches(stex::Patch * patcha, stex::Patch * patchb)
{
    AssertLog(pSetPatches == false);
    AssertLog(patcha != nullptr);
    AssertLog(patchb != nullptr);
    AssertLog(patcha != patchb);

    pPatchA = patcha;
    pPatchB = patchb;
    pSetPatches = true;
}

////////////////////////////////////////////////////////////////////////////////

stex::Patch * stex::SDiffBoundary::patchA()
{
    AssertLog(pSetPatches == true);
    return pPatchA;
}

////////////////////////////////////////////////////////////////////////////////

stex::Patch * stex::SDiffBoundary::patchB()
{
    AssertLog(pSetPatches == true);
    return pPatchB;
}

////////////////////////////////////////////////////////////////////////////////

void stex::SDiffBoundary::setTriDirection(triangle_id_t tri, uint direction)
{
    AssertLog(direction < 3);

    pTris.push_back(tri);
    pTriDirection.push_back(direction);
}

////////////////////////////////////////////////////////////////////////////////

// END

