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
#include <string>
#include <cassert>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/types.hpp"
#include "steps/error.hpp"
#include "steps/solver/statedef.hpp"
#include "steps/solver/sdiffboundarydef.hpp"
#include "steps/geom/patch.hpp"

// logging
#include "easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace ssolver = steps::solver;
namespace stetmesh = steps::tetmesh;

////////////////////////////////////////////////////////////////////////////////

ssolver::SDiffBoundarydef::SDiffBoundarydef(Statedef * sd, uint idx, stetmesh::SDiffBoundary * sdb)
: pStatedef(sd)
, pIdx(idx)
, pName()
, pBars()
, pPatchA_temp(0)
, pPatchB_temp(0)
, pPatchA(0)
, pPatchB(0)
, pSetupdone(false)
{
    AssertLog(pStatedef != 0);
    AssertLog(sdb != 0);

    pName = sdb->getID();
    pBars = sdb->_getAllBarIndices();
    std::vector<steps::wm::Patch *> patches = sdb->getPatches();
    pPatchA_temp = patches[0];
    pPatchB_temp = patches[1];
    AssertLog(pPatchA_temp != 0);
    AssertLog(pPatchB_temp != 0);

}

////////////////////////////////////////////////////////////////////////////////

ssolver::SDiffBoundarydef::~SDiffBoundarydef(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::SDiffBoundarydef::checkpoint(std::fstream & cp_file)
{
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::SDiffBoundarydef::restore(std::fstream & cp_file)
{
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::SDiffBoundarydef::setup(void)
{
    AssertLog(pSetupdone == false);

    pPatchA = pStatedef->getPatchIdx(pPatchA_temp);
    pPatchB = pStatedef->getPatchIdx(pPatchB_temp);
    AssertLog(pPatchA >= 0);
    AssertLog(pPatchB >= 0);
    pSetupdone = true;

}

////////////////////////////////////////////////////////////////////////////////

std::string const ssolver::SDiffBoundarydef::name(void) const
{
    return pName;
}

////////////////////////////////////////////////////////////////////////////////


