/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
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


// Standard library & STL headers.
// #include <vector>
#include <algorithm>

// STEPS headers.
#include "steps/common.h"
#include "steps/wmdirect/patch.hpp"
#include "steps/wmdirect/comp.hpp"
#include "steps/wmdirect/kproc.hpp"
#include "steps/wmdirect/sreac.hpp"
#include "steps/wmdirect/wmdirect.hpp"
#include "steps/solver/statedef.hpp"
#include "steps/solver/patchdef.hpp"
#include "steps/solver/types.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace swmd = steps::wmdirect;
namespace ssolver = steps::solver;

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

void swmd::Patch::checkpoint(std::fstream & cp_file)
{
    for (KProcPVecCI k = pKProcs.begin(); k != pKProcs.end(); ++k)
    {
        (*k)->checkpoint(cp_file);
    }
}

////////////////////////////////////////////////////////////////////////////////

void swmd::Patch::restore(std::fstream & cp_file)
{
    for (KProcPVecCI k = pKProcs.begin(); k != pKProcs.end(); ++k)
    {
        (*k)->restore(cp_file);
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
