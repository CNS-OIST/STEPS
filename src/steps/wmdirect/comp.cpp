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
#include "steps/solver/compdef.hpp"
#include "steps/wmdirect/comp.hpp"
#include "steps/wmdirect/kproc.hpp"
#include "steps/wmdirect/reac.hpp"
#include "steps/wmdirect/wmdirect.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace swmd = steps::wmdirect;
namespace ssolver = steps::solver;

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

void swmd::Comp::checkpoint(std::fstream & cp_file)
{
    for (KProcPVecCI k = pKProcs.begin(); k != pKProcs.end(); ++k)
    {
        (*k)->checkpoint(cp_file);
    }
}

////////////////////////////////////////////////////////////////////////////////

void swmd::Comp::restore(std::fstream & cp_file)
{
    for (KProcPVecCI k = pKProcs.begin(); k != pKProcs.end(); ++k)
    {
        (*k)->restore(cp_file);
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
