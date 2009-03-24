////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2008 Stefan Wils. All rights reserved.
//
// This file is part of STEPS.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA
//
//
////////////////////////////////////////////////////////////////////////////////

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// Standard library & STL headers.
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/solver/compdef.hpp>
#include <steps/tetexact/patch.hpp>
#include <steps/tetexact/kproc.hpp>
#include <steps/tetexact/reac.hpp>
#include <steps/tetexact/tri.hpp>

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
        if (selector < accum) return *t;
    }

    assert(false);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

// END
