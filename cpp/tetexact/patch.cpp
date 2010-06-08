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
#include "patch.hpp"
#include "kproc.hpp"
#include "reac.hpp"
#include "tri.hpp"

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::tetexact, stex);
NAMESPACE_ALIAS(steps::solver, ssolver);

////////////////////////////////////////////////////////////////////////////////

stex::Patch::Patch(ssolver::Patchdef * patchdef)
: pPatchdef(patchdef)
, pTris()
{
    assert(pPatchdef != 0);
}

////////////////////////////////////////////////////////////////////////////////

stex::Patch::~Patch(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void stex::Patch::addTri(stex::Tri * tri)
{
    assert(tri->patchdef() == def());
    pTris.push_back(tri);
    pArea += tri->area();
}

////////////////////////////////////////////////////////////////////////////////

void stex::Patch::modCount(uint slidx, double count)
{
	assert (slidx < def()->countSpecs());
	double newcount = (def()->pools()[slidx] + count);
	assert (newcount >= 0.0);
    def()->setCount(slidx, newcount);
}

////////////////////////////////////////////////////////////////////////////////

stex::Tri * stex::Patch::pickTriByArea(double rand01) const
{
    if (countTris() == 0) return 0;
    if (countTris() == 1) return pTris[0];

    double accum = 0.0;
    double selector = rand01 * area();
    TriPVecCI t_end = endTri();
    for (TriPVecCI t = bgnTri(); t != t_end; ++t)
    {
        accum += (*t)->area();
        if (selector <= accum) return *t;
    }

    assert(false);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

// END
