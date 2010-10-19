////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2009ÊOkinawa Institute of Science and Technology, Japan.
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

// Standard library & STL headers.
#include <vector>

// STEPS headers.
#include "../common.h"
#include "../solver/diffboundarydef.hpp"
#include "diffboundary.hpp"
#include "comp.hpp"

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::tetexact, stex);
NAMESPACE_ALIAS(steps::solver, ssolver);

////////////////////////////////////////////////////////////////////////////////

stex::DiffBoundary::DiffBoundary(steps::solver::DiffBoundarydef * dbdef)
: pDiffBoundarydef(dbdef)
, pCompA(0)
, pCompB(0)
, pTets()
, pTetDirection()
, pSetComps(false)
{
	assert(dbdef != 0);
}

////////////////////////////////////////////////////////////////////////////////

stex::DiffBoundary::~DiffBoundary(void)
{

}

////////////////////////////////////////////////////////////////////////////////

void stex::DiffBoundary::setComps(stex::Comp * compa, stex::Comp * compb)
{
	assert (pSetComps == false);
    assert(compa != 0);
    assert(compb != 0);
    assert(compa != compb);

    pCompA = compa;
    pCompB = compb;
    pSetComps = true;
}

////////////////////////////////////////////////////////////////////////////////

stex::Comp * stex::DiffBoundary::compA(void)
{
	assert(pSetComps == true);
	return pCompA;
}

////////////////////////////////////////////////////////////////////////////////

stex::Comp * stex::DiffBoundary::compB(void)
{
	assert(pSetComps == true);
	return pCompB;
}

////////////////////////////////////////////////////////////////////////////////

void stex::DiffBoundary::setTetDirection(uint tet, uint direction)
{
	assert(direction < 4);

	pTets.push_back(tet);
	pTetDirection.push_back(direction);
}

////////////////////////////////////////////////////////////////////////////////

// END

