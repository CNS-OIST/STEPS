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
// #include <vector>
#include <algorithm>

// STEPS headers.
#include "../common.h"
#include "../solver/compdef.hpp"
#include "comp.hpp"
#include "kproc.hpp"
#include "reac.hpp"
#include "wmdirect.hpp"

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::wmdirect, swmd);
NAMESPACE_ALIAS(steps::solver, ssolver);

////////////////////////////////////////////////////////////////////////////////

swmd::Comp::Comp(steps::solver::Compdef * compdef)
: pCompdef(compdef)
, pKProcs()
, pIPatches()
, pOPatches()
{
	assert (pCompdef != 0);
}

////////////////////////////////////////////////////////////////////////////////

swmd::Comp::~Comp(void)
{
    for (KProcPVecCI k = pKProcs.begin(); k != pKProcs.end(); ++k)
    {
        delete (*k);
    }
}

////////////////////////////////////////////////////////////////////////////////

void swmd::Comp::reset(void)
{
    std::for_each(pKProcs.begin(), pKProcs.end(), std::mem_fun(&swmd::KProc::reset));
}

////////////////////////////////////////////////////////////////////////////////

void swmd::Comp::setupKProcs(swmd::Wmdirect * wmd)
{
    // Create reaction kproc's.
    uint nreacs = def()->countReacs();
    pKProcs.resize(nreacs);
    for (uint i = 0; i < nreacs; ++i)
    {
        ssolver::Reacdef * rdef = def()->reacdef(i);
        swmd::Reac * r = new swmd::Reac(rdef, this);
        pKProcs[i] = r;
        wmd->addKProc(r);
    }
}

////////////////////////////////////////////////////////////////////////////////

steps::wmdirect::KProc * swmd::Comp::reac(uint lridx) const
{
	assert (lridx < pKProcs.size());
	return pKProcs[lridx];
}

////////////////////////////////////////////////////////////////////////////////

void swmd::Comp::setupDeps(void)
{
    std::for_each(pKProcs.begin(), pKProcs.end(),
        std::mem_fun(&KProc::setupDeps));
}

////////////////////////////////////////////////////////////////////////////////


void swmd::Comp::addIPatch(swmd::Patch * p)
{
    assert(std::find(pIPatches.begin(), pIPatches.end(), p) == pIPatches.end());
    pIPatches.push_back(p);
}

////////////////////////////////////////////////////////////////////////////////

void swmd::Comp::addOPatch(swmd::Patch * p)
{
    assert(std::find(pOPatches.begin(), pOPatches.end(), p) == pOPatches.end());
    pOPatches.push_back(p);
}
////////////////////////////////////////////////////////////////////////////////

// END
