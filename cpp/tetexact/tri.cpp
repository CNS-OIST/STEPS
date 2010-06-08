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
#include <cassert>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iostream>

// STEPS headers.
#include "../common.h"
#include "../solver/patchdef.hpp"
#include "../solver/sreacdef.hpp"
#include "sreac.hpp"
#include "tet.hpp"
#include "tri.hpp"
#include "kproc.hpp"
#include "tetexact.hpp"

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::tetexact, stex);
NAMESPACE_ALIAS(steps::solver, ssolver);

////////////////////////////////////////////////////////////////////////////////

stex::Tri::Tri(steps::solver::Patchdef * patchdef, double area,
			   int tetinner, int tetouter)
: pPatchdef(patchdef)
, pArea(area)
, pInnerTet(0)
, pOuterTet(0)
, pTets()
, pPoolCount(0)
, pPoolFlags(0)
, pKProcs()
{
	assert(pPatchdef != 0);
	assert (pArea > 0.0);

	pTets[0] = tetinner;
	pTets[1] = tetouter;

	uint nspecs = pPatchdef->countSpecs();
    pPoolCount = new uint[nspecs];
    pPoolFlags = new uint[nspecs];
    std::fill_n(pPoolCount, nspecs, 0);
    std::fill_n(pPoolFlags, nspecs, 0);
    pKProcs.resize(pPatchdef->countSReacs());
}

////////////////////////////////////////////////////////////////////////////////

stex::Tri::~Tri(void)
{
    delete[] pPoolCount;
    delete[] pPoolFlags;
    KProcPVecCI e = pKProcs.end();
    for (std::vector<stex::KProc *>::const_iterator i = pKProcs.begin();
         i != e; ++i) delete *i;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tri::setInnerTet(stex::Tet * t)
{
	pInnerTet = t;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tri::setOuterTet(stex::Tet * t)
{
	pOuterTet = t;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tri::setupKProcs(stex::Tetexact * tex)
{
	uint j = 0;
	// Create surface reaction kprocs
	uint nsreacs = patchdef()->countSReacs();
	for (uint i=0; i < nsreacs; ++i)
	{
		ssolver::SReacdef * srdef = patchdef()->sreacdef(i);
		stex::SReac * sr = new SReac(srdef, this);
		assert(sr != 0);
		pKProcs[j++] = sr;
		tex->addKProc(sr);
	}
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tri::reset(void)
{
    uint nspecs = patchdef()->countSpecs();
    std::fill_n(pPoolCount, nspecs, 0);
    std::fill_n(pPoolFlags, nspecs, 0);
    std::for_each(pKProcs.begin(), pKProcs.end(),
        std::mem_fun(&stex::KProc::reset));
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tri::setCount(uint lidx, uint count)
{
	assert (lidx < patchdef()->countSpecs());
	double oldcount = pPoolCount[lidx];
	double c = static_cast<double>(count);
	pPoolCount[lidx] = c;

	/* 16/01/10 IH: Counts no longer stored in patch object.
	// Now update the count in this tri's patch
	double diff = c - oldcount;
	double newcount = (patchdef()->pools()[lidx]) + diff;
	// Patchdef method will do the checking on the double argument (should be positive!)
	patchdef()->setCount(lidx, newcount);
	*/
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tri::setClamped(uint lidx, bool clamp)
{
	if (clamp == true) pPoolFlags[lidx] |= CLAMPED;
    else pPoolFlags[lidx] &= ~CLAMPED;
}

////////////////////////////////////////////////////////////////////////////////

stex::SReac * stex::Tri::sreac(uint lidx) const
{
    assert(lidx < patchdef()->countSReacs());
    return dynamic_cast<stex::SReac*>(pKProcs[lidx]);
}

////////////////////////////////////////////////////////////////////////////////

//END
