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
#include "steps/solver/compdef.hpp"
#include "steps/wmrssa/comp.hpp"
#include "steps/wmrssa/kproc.hpp"
#include "steps/wmrssa/reac.hpp"
#include "steps/wmrssa/wmrssa.hpp"

// logging
#include "steps/error.hpp"
#include "easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

namespace swmrssa = steps::wmrssa;
namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

swmrssa::Comp::Comp(steps::solver::Compdef * compdef)
: pCompdef(compdef)
, pKProcs()
, pIPatches()
, pOPatches()
{
    assert (pCompdef != 0);
    uint nspecs = compdef->countSpecs();
    pPoolLB = new double[nspecs];
    pPoolUB = new double[nspecs];
}

////////////////////////////////////////////////////////////////////////////////

swmrssa::Comp::~Comp(void)
{
    for (KProcPVecCI k = pKProcs.begin(); k != pKProcs.end(); ++k)
    {
        delete (*k);
    }
    delete[] pPoolLB;
    delete[] pPoolUB;
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::Comp::checkpoint(std::fstream & cp_file)
{
    for (KProcPVecCI k = pKProcs.begin(); k != pKProcs.end(); ++k)
    {
        (*k)->checkpoint(cp_file);
    }
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::Comp::restore(std::fstream & cp_file)
{
    for (KProcPVecCI k = pKProcs.begin(); k != pKProcs.end(); ++k)
    {
        (*k)->restore(cp_file);
    }
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::Comp::reset(void)
{
    std::for_each(pKProcs.begin(), pKProcs.end(), std::mem_fun(&swmrssa::KProc::reset));
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::Comp::setupKProcs(swmrssa::Wmrssa * wmd)
{
    // Create reaction kproc's.
    uint nreacs = def()->countReacs();
    pKProcs.resize(nreacs);
    for (uint i = 0; i < nreacs; ++i)
    {
        ssolver::Reacdef * rdef = def()->reacdef(i);
        swmrssa::Reac * r = new swmrssa::Reac(rdef, this);
        pKProcs[i] = r;
        wmd->addKProc(r);
    }
}

////////////////////////////////////////////////////////////////////////////////

steps::wmrssa::KProc * swmrssa::Comp::reac(uint lridx) const
{
    assert (lridx < pKProcs.size());
    return pKProcs[lridx];
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::Comp::setupDeps(void)
{
    std::for_each(pKProcs.begin(), pKProcs.end(),
        std::mem_fun(&KProc::setupDeps));
}

////////////////////////////////////////////////////////////////////////////////


void swmrssa::Comp::addIPatch(swmrssa::Patch * p)
{
    AssertLog(std::find(pIPatches.begin(), pIPatches.end(), p) == pIPatches.end());
    pIPatches.push_back(p);
}

////////////////////////////////////////////////////////////////////////////////

void swmrssa::Comp::addOPatch(swmrssa::Patch * p)
{
    AssertLog(std::find(pOPatches.begin(), pOPatches.end(), p) == pOPatches.end());
    pOPatches.push_back(p);
}
////////////////////////////////////////////////////////////////////////////////

void swmrssa::Comp::setBounds(uint i, int nc)
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
    /*pPoolLB[i] = std::max(std::min(nc*(1 - delta), nc - 3.), 0.); // nc/(7 - delta - 2*(3 - delta)/(1 + 2./(nc+1))); //*(1 - delta);//
    pPoolUB[i] = nc > 0 ? nc*(7 - delta - 2*(3 - delta)/(1 + 2./(nc+1))) : 3; //std::max(nc*(1 + delta), nc + 3.); // (1 + delta);/*/
}
////////////////////////////////////////////////////////////////////////////////

bool swmrssa::Comp::isOutOfBound(uint i, int nc)
{
    AssertLog(i < def()->countSpecs());
    if (nc > pPoolLB[i] && nc < pPoolUB[i])
        return false;
    setBounds(i, nc);
    return true;
}

double* swmrssa::Comp::pools(steps::wmrssa::PropensityRSSA prssa) const
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

void swmrssa::Comp::setupSpecDeps(void)
{
    uint nspecs = def()->countSpecs();
    localSpecUpdKProcs.resize(nspecs);
    for (uint slidx = 0; slidx < nspecs; slidx++) {
        uint sgidx = def()->specL2G(slidx);
        for (KProcPVecCI k = pKProcs.begin(); k != pKProcs.end(); ++k)
        {
            if ((*k)->depSpecComp(sgidx, this) == true) {
                localSpecUpdKProcs[slidx].push_back(*k);
            }
        }
        PatchPVecCI ip_bgn = beginIPatches();
        PatchPVecCI ip_end = endIPatches();
        KProcPVecCI kprocend;
        for (PatchPVecCI ip = ip_bgn; ip != ip_end; ++ip)
        {
            kprocend = (*ip)->kprocEnd();
            for (KProcPVecCI k = (*ip)->kprocBegin(); k != kprocend; ++k)
            {
                if ((*k)->depSpecComp(sgidx, this) == true)
                    localSpecUpdKProcs[slidx].push_back(*k);
            }
        }
        ip_bgn = beginOPatches();
        ip_end = endOPatches();
        for (PatchPVecCI ip = ip_bgn; ip != ip_end; ++ip)
        {
            kprocend = (*ip)->kprocEnd();
            for (KProcPVecCI k = (*ip)->kprocBegin(); k != kprocend; ++k)
            {
                if ((*k)->depSpecComp(sgidx, this) == true)
                    localSpecUpdKProcs[slidx].push_back(*k);
            }
        }
    }
}

std::vector<swmrssa::KProc*> const & swmrssa::Comp::getSpecUpdKProcs(uint slidx)
{
    return localSpecUpdKProcs[slidx];
}

// END
