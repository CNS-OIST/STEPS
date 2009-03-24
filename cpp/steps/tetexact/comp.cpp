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
#include <steps/tetexact/comp.hpp>
#include <steps/tetexact/kproc.hpp>
#include <steps/tetexact/reac.hpp>
#include <steps/tetexact/tet.hpp>

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
		if (selector < accum) return (*t);
	}
	assert(false);
	return 0;
}

////////////////////////////////////////////////////////////////////////////////

// END
