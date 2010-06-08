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

/*
 *  Last Changed Rev:  $Rev: 307 $
 *  Last Changed Date: $Date: 2010-03-24 09:25:37 +0900 (Wed, 24 Mar 2010) $
 *  Last Changed By:   $Author: wchen $
 */


// Standard library & STL headers.
#include <vector>

// STEPS headers.
#include "../common.h"
#include "../solver/compdef.hpp"
#include "comp.hpp"
#include "kproc.hpp"
#include "reac.hpp"
#include "tet.hpp"

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::tetexact, stex);
NAMESPACE_ALIAS(steps::solver, ssolver);

////////////////////////////////////////////////////////////////////////////////

stex::Comp::Comp(steps::solver::Compdef * compdef)
: pCompdef(compdef)
, pVol(0.0)
, pTets()
{
	assert(pCompdef != 0);
}

////////////////////////////////////////////////////////////////////////////////

stex::Comp::~Comp(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void stex::Comp::addTet(stex::Tet * tet)
{
	assert (tet->compdef() == def());
	pTets.push_back(tet);
	pVol += tet->vol();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Comp::modCount(uint slidx, double count)
{
	assert (slidx < def()->countSpecs());
	double newcount = (def()->pools()[slidx] + count);
	assert (newcount >= 0.0);
    def()->setCount(slidx, newcount);
}

////////////////////////////////////////////////////////////////////////////////

stex::Tet * stex::Comp::pickTetByVol(double rand01) const
{
	if (countTets() == 0) return 0;
	if (countTets() == 1) return pTets[0];

	double accum = 0.0;
	double selector = rand01 * vol();
	TetPVecCI t_end = endTet();
	for (TetPVecCI t = bgnTet(); t != t_end; ++t)
	{
		accum += (*t)->vol();
		if (selector <= accum) return (*t);
	}

	assert(false);
	return 0;
}

////////////////////////////////////////////////////////////////////////////////

// END
