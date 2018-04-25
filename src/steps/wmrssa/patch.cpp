/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
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
#include "steps/wmrssa/patch.hpp"
#include "steps/wmrssa/comp.hpp"
#include "steps/wmrssa/kproc.hpp"
#include "steps/wmrssa/sreac.hpp"
#include "steps/wmrssa/wmrssa.hpp"
#include "steps/solver/statedef.hpp"
#include "steps/solver/patchdef.hpp"
#include "steps/solver/types.hpp"

// logging
#include "steps/error.hpp"
#include "easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

namespace swmrssa = steps::wmrssa;
namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

swmrssa::Patch::Patch(steps::solver::Patchdef * patchdef, swmrssa::Comp * icomp, swmrssa::Comp * ocomp)
: pPatchdef(patchdef)
, pKProcs()
, pIComp(icomp)
, pOComp(ocomp)
{
    AssertLog(pPatchdef != 0);
    if (iComp() != 0) iComp()->addIPatch(this);
    if (oComp() != 0) oComp()->addOPatch(this);
    uint nspecs = patchdef->countSpecs();
    pPoolLB = new double[nspecs];
    pPoolUB = new double[nspecs];
}

////////////////////////////////////////////////////////////////////////////////

swmrssa::Patch::~Patch(void)
{
    for (KProcPVecCI k = pKProcs.begin(); k != pKProcs.end(); ++k)
    {
        delete (*k);
    }
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::Patch::checkpoint(std::fstream & cp_file)
{
    for (KProcPVecCI k = pKProcs.begin(); k != pKProcs.end(); ++k)
    {
        (*k)->checkpoint(cp_file);
    }
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::Patch::restore(std::fstream & cp_file)
{
    for (KProcPVecCI k = pKProcs.begin(); k != pKProcs.end(); ++k)
    {
        (*k)->restore(cp_file);
    }
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::Patch::setupKProcs(swmrssa::Wmrssa * wmd)
{
    // Create surface reaction kproc's.
    uint nsreacs = def()->countSReacs();
    pKProcs.resize(nsreacs);
    for (uint i = 0; i < nsreacs; ++i)
    {
        ssolver::SReacdef * srdef = def()->sreacdef(i);
        swmrssa::SReac * sr = new swmrssa::SReac(srdef, this);
        pKProcs[i] = sr;
        wmd->addKProc(sr);
    }
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::Patch::setupDeps(void)
{
    std::for_each(pKProcs.begin(), pKProcs.end(),
        std::mem_fun(&KProc::setupDeps));
}

////////////////////////////////////////////////////////////////////////////////

swmrssa::KProc * swmrssa::Patch::sreac(uint lsridx) const
{
    assert (lsridx < pKProcs.size());
    return pKProcs[lsridx];
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::Patch::reset(void)
{
    std::for_each(pKProcs.begin(), pKProcs.end(), std::mem_fun(&KProc::reset));
}


////////////////////////////////////////////////////////////////////////////////

void swmrssa::Patch::setBounds(uint i, int nc)
{
    const double delta = .05;
    if (nc > 3/delta)
    {
        pPoolLB[i] = nc*(1 - delta);
        pPoolUB[i] = nc*(1 + delta);
    }
    else if (nc > 3)
    {
        pPoolLB[i] = nc - 3;
        pPoolUB[i] = nc + 3;
    }
    else if (nc > 0)
    {
        pPoolLB[i] = 1;
        pPoolUB[i] = 2*nc;
    }
    else{
        pPoolLB[i] = 0;
        pPoolUB[i] = 0;
    }
    pPoolLB[i] -= delta;
    pPoolUB[i] += delta;
    /*pPoolLB[i] = std::max(std::min(nc*(1 - delta), nc - 3.), 0.); //nc/(3 - delta - 2*(1 - delta)/(1 + 2./(nc+1))); //
    pPoolUB[i] = std::max(nc*(1 + delta), nc + 3.); //nc*(3 - delta - 2*(1 - delta)/(1 + 2./(nc+1))); /*/
}

////////////////////////////////////////////////////////////////////////////////

bool swmrssa::Patch::isOutOfBound(uint i, int nc)
{
    AssertLog(i < def()->countSpecs());
    if (nc > pPoolLB[i] && nc < pPoolUB[i])
        return false;
    setBounds(i, nc);
    return true;
}

////////////////////////////////////////////////////////////////////////////////

double* swmrssa::Patch::pools(steps::wmrssa::PropensityRSSA prssa) const
{
    switch(prssa)
    {
        case steps::wmrssa::CURRENT:
            return def()->pools();
        case steps::wmrssa::LOWERBOUND:
            return pPoolLB;
        case steps::wmrssa::BOUNDS:
            return pPoolUB;
    }
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::Patch::setupSpecDeps(void)
{
    uint nspecs = def()->countSpecs();
    localSpecUpdKProcs.resize(nspecs);
    for (uint slidx = 0; slidx < nspecs; slidx++) {
        uint sgidx = def()->specL2G(slidx);
        for (KProcPVecCI k = pKProcs.begin(); k != pKProcs.end(); ++k)
        {
            if ((*k)->depSpecPatch(sgidx, this) == true) {
                localSpecUpdKProcs[slidx].push_back(*k);
            }
        }
        if (pIComp)
        {
            KProcPVecCI kprocend = pIComp->kprocEnd();
            for (KProcPVecCI k = pIComp->kprocBegin(); k != kprocend; ++k) {
                if ((*k)->depSpecPatch(sgidx, this) == true) {
                    localSpecUpdKProcs[slidx].push_back(*k);
                }
            }
        }
        if (pOComp)
        {
            KProcPVecCI kprocend = pOComp->kprocEnd();
            for (KProcPVecCI k = pOComp->kprocBegin(); k != kprocend; ++k) {
                if ((*k)->depSpecPatch(sgidx, this) == true) {
                    localSpecUpdKProcs[slidx].push_back(*k);
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

std::vector<swmrssa::KProc*> const & swmrssa::Patch::getSpecUpdKProcs(uint slidx)
{
    return localSpecUpdKProcs[slidx];
}

// END
