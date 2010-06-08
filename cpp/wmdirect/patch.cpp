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
#include "patch.hpp"
#include "comp.hpp"
#include "kproc.hpp"
#include "sreac.hpp"
#include "wmdirect.hpp"
#include "../solver/statedef.hpp"
#include "../solver/patchdef.hpp"
#include "../solver/types.hpp"

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::wmdirect, swmd);
NAMESPACE_ALIAS(steps::solver, ssolver);

////////////////////////////////////////////////////////////////////////////////

swmd::Patch::Patch(steps::solver::Patchdef * patchdef, swmd::Comp * icomp, swmd::Comp * ocomp)
: pPatchdef(patchdef)
, pKProcs()
, pIComp(icomp)
, pOComp(ocomp)
{
	assert(pPatchdef != 0);
    if (iComp() != 0) iComp()->addIPatch(this);
    if (oComp() != 0) oComp()->addOPatch(this);
}

////////////////////////////////////////////////////////////////////////////////

swmd::Patch::~Patch(void)
{
    for (KProcPVecCI k = pKProcs.begin(); k != pKProcs.end(); ++k)
    {
        delete (*k);
    }
}

////////////////////////////////////////////////////////////////////////////////

void swmd::Patch::setupKProcs(swmd::Wmdirect * wmd)
{
    // Create surface reaction kproc's.
    uint nsreacs = def()->countSReacs();
    pKProcs.resize(nsreacs);
    for (uint i = 0; i < nsreacs; ++i)
    {
        ssolver::SReacdef * srdef = def()->sreacdef(i);
        swmd::SReac * sr = new swmd::SReac(srdef, this);
        pKProcs[i] = sr;
        wmd->addKProc(sr);
    }
}

////////////////////////////////////////////////////////////////////////////////

void swmd::Patch::setupDeps(void)
{
    std::for_each(pKProcs.begin(), pKProcs.end(),
        std::mem_fun(&KProc::setupDeps));
}

////////////////////////////////////////////////////////////////////////////////

swmd::KProc * swmd::Patch::sreac(uint lsridx) const
{
	assert (lsridx < pKProcs.size());
	return pKProcs[lsridx];
}

////////////////////////////////////////////////////////////////////////////////

void swmd::Patch::reset(void)
{
    std::for_each(pKProcs.begin(), pKProcs.end(), std::mem_fun(&KProc::reset));
}

////////////////////////////////////////////////////////////////////////////////

// END
