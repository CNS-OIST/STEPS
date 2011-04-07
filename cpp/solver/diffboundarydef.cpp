////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2011ÊOkinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006ÊUniversity of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPSÊis free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPSÊis distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.ÊIf not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////

// STL headers.
#include <string>
#include <cassert>

// STEPS headers.
#include "../common.h"
#include "types.hpp"
#include "../error.hpp"
#include "statedef.hpp"
#include "diffboundarydef.hpp"
#include "../geom/comp.hpp"

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::solver, ssolver);
NAMESPACE_ALIAS(steps::tetmesh, stetmesh);

////////////////////////////////////////////////////////////////////////////////

ssolver::DiffBoundarydef::DiffBoundarydef(Statedef * sd, uint idx, stetmesh::DiffBoundary * db)
: pStatedef(sd)
, pIdx(idx)
, pName()
, pTris()
, pCompA_temp(0)
, pCompB_temp(0)
, pCompA(0)
, pCompB(0)
, pSetupdone(false)
{
    assert(pStatedef != 0);
    assert(db != 0);

    pName = db->getID();
    pTris = db->_getAllTriIndices();
    std::vector<steps::wm::Comp *> comps = db->getComps();
    pCompA_temp = comps[0];
    pCompB_temp = comps[1];
    assert (pCompA_temp != 0);
    assert (pCompB_temp != 0);

}

////////////////////////////////////////////////////////////////////////////////

ssolver::DiffBoundarydef::~DiffBoundarydef(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::DiffBoundarydef::checkpoint(std::fstream & cp_file)
{
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::DiffBoundarydef::restore(std::fstream & cp_file)
{
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::DiffBoundarydef::setup(void)
{
	assert (pSetupdone == false);

	pCompA = pStatedef->getCompIdx(pCompA_temp);
	pCompB = pStatedef->getCompIdx(pCompB_temp);
	assert(pCompA >= 0);
	assert(pCompB >= 0);
	pSetupdone = true;

}

////////////////////////////////////////////////////////////////////////////////

std::string const ssolver::DiffBoundarydef::name(void) const
{
	return pName;
}

////////////////////////////////////////////////////////////////////////////////


